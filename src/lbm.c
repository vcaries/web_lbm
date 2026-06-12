/*
 * lbm.c — D2Q9 Lattice Boltzmann wind-tunnel engine.
 *
 * Compiled to WebAssembly with Emscripten (see build.sh). All field arrays
 * are allocated here with malloc and exposed to JavaScript as raw pointers;
 * JS builds TypedArray views directly over WebAssembly.Memory — zero copies.
 *
 * Collision operator: BGK (single relaxation time) by default.
 * Compile with -DUSE_MRT for the multiple-relaxation-time operator
 * (Lallemand & Luo 2000), which is markedly more stable at high Reynolds.
 *
 * Lattice (D2Q9), velocity set e_i and weights w_i:
 *
 *        6   2   5          i : (ex, ey)  w
 *         \  |  /           0 : ( 0, 0)  4/9
 *      3 --- 0 --- 1        1 : ( 1, 0)  1/9   2 : (0, 1)  1/9
 *         /  |  \           3 : (-1, 0)  1/9   4 : (0,-1)  1/9
 *        7   4   8          5 : ( 1, 1)  1/36  6 : (-1, 1) 1/36
 *                           7 : (-1,-1)  1/36  8 : ( 1,-1) 1/36
 *
 * Grid indexing: idx = y * Nx + x, x in [0, Nx), y in [0, Ny).
 * y = 0 is the bottom wall. Flow enters at x = 0 (inlet), exits at x = Nx-1.
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ------------------------------------------------------------------ */
/* Lattice constants                                                   */
/* ------------------------------------------------------------------ */

static const int   ex[9]  = { 0, 1, 0,-1, 0, 1,-1,-1, 1 };
static const int   ey[9]  = { 0, 0, 1, 0,-1, 1, 1,-1,-1 };
static const float w9[9]  = { 4.f/9.f,
                              1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f,
                              1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f };
/* opp[i]: direction opposite to i (used by bounce-back). */
static const int   opp[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
/* mirY[i]: direction i with its y-component negated (specular wall). */
static const int   mirY[9]= { 0, 1, 4, 3, 2, 8, 7, 6, 5 };

/* Obstacle mask values. Painted and preset cells are both solid for the
 * physics, but lbm_resize() treats them differently: painted cells are
 * resampled nearest-neighbour, preset cells are re-rasterised analytically
 * from the shape registry so they stay crisp at any resolution. */
#define OBS_FLUID   0
#define OBS_PAINTED 1
#define OBS_PRESET  2

/* ------------------------------------------------------------------ */
/* Simulation state (Structure of Arrays)                              */
/* ------------------------------------------------------------------ */

static int Nx = 0, Ny = 0;

/* For each of the 9 velocity directions, a flat array of length Nx*Ny.   */
static float  *f[9];      /* current distributions                        */
static float  *f_tmp[9];  /* post-stream buffer (ping-pong, swapped)      */
static uint8_t *obstacle; /* 0 = fluid, 1 = painted solid, 2 = preset     */

/* Macroscopic fields, recomputed on demand by lbm_compute_fields().
 * schl_f = |grad rho| (numerical schlieren); dil_f = div u (dilatation,
 * isolates the acoustic radiation: it vanishes for incompressible
 * vortical motion). */
static float *rho_f, *ux_f, *uy_f, *speed_f, *vort_f, *schl_f, *dil_f;

/* Wall velocity of solid cells (moving-wall bounce-back). Zero everywhere
 * except on the footprint of rotating cylinder presets, where it holds the
 * local solid-body rotation velocity. */
static float *uwx, *uwy;

/* Physical parameters (lattice units). */
static float u0_mag   = 0.10f;  /* inlet speed                  */
static float angle_rad = 0.0f;  /* inlet flow angle             */
static float tau      = 0.55f;  /* BGK relaxation time          */
static int   wall_noslip = 0;   /* 0 = free-slip, 1 = no-slip   */

/* Smagorinsky subgrid model (LES). At low tau the bare collision operator
 * develops checkerboard pressure oscillations; the eddy viscosity
 * nu_t = (Cs*dx)^2 |S| grows exactly where gradients steepen and damps
 * them. les_on toggles it at runtime; smag_c2 = Cs^2. */
static int   les_on  = 1;
static float smag_c2 = 0.18f * 0.18f;

/* Hermite-regularized collision (Latt & Chopard 2006). Before relaxing,
 * the non-equilibrium part of each population is replaced by its
 * projection onto the second-order Hermite (stress) basis:
 *
 *   f_i^neq -> f_i^(1) = (w_i / 2 c_s^4) Q_i : Pi^neq,
 *   Q_i = e_i e_i - c_s^2 I,   Pi^neq = sum_i e_i e_i (f_i - f_i^eq).
 *
 * Everything orthogonal to the hydrodynamic moments (the "ghost" modes
 * that destabilise BGK at low tau) is discarded every step. Costs ~20
 * flops/cell, keeps the kernel branchless, and combines with the LES. */
static int regularized = 1;

/* Bulk-viscosity boost, applied inside the regularized reconstruction:
 * the trace of Pi^neq (the dilatational / acoustic stress) is scaled by
 * BULK_DAMP < 1, which relaxes the acoustic modes faster than the shear
 * modes (an enhanced bulk viscosity). Damping grows like k^2, so the
 * under-resolved short-wavelength "wrinkles" the lattice supports die
 * fastest while vortical motion — carried by the traceless part — is
 * untouched. */
#define BULK_DAMP 0.4f

/* Soft start: the inlet/initial speed ramps up over the first RAMP_STEPS
 * steps after a (re)initialisation — an impulsive start against an
 * obstacle is otherwise the most fragile moment at high Re. */
#define RAMP_STEPS 400
static int ramp_t = RAMP_STEPS;

static inline float ramp_factor(void)
{
    return ramp_t >= RAMP_STEPS ? 1.f
                                : ((float)ramp_t + 1.f) / (float)RAMP_STEPS;
}

/* Which derivative fields lbm_compute_fields() fills, as a bitmask:
 * 1 = vorticity, 2 = schlieren, 4 = dilatation. The moments (rho, u,
 * speed) are always computed; the UI displays one field at a time, so
 * it narrows the mask to skip the unused finite-difference passes. */
#define FIELD_VORT 1
#define FIELD_SCHL 2
#define FIELD_DIL  4
static int field_mask = FIELD_VORT | FIELD_SCHL | FIELD_DIL;

/* Inlet flow profile: 0 = uniform, 1 = shear layer, 2 = jet.
 * Applied both to the initial condition and to the Zou-He inlet, so the
 * selected profile is sustained, not just a transient. */
static int flow_profile = 0;

/* Streamwise inlet speed at row y (fraction of the full grid height). */
static inline float profile_mag(int y)
{
    float yn = ((float)y + 0.5f) / (float)Ny;
    switch (flow_profile) {
        case 1: /* shear layer: tanh ramp between 0.1*u0 and u0 */
            return u0_mag * (0.55f + 0.45f * tanhf(8.f * (yn - 0.5f)));
        case 2: { /* jet: gaussian core (width 0.12 Ny) + 15% co-flow */
            float g = (yn - 0.5f) / 0.12f;
            return u0_mag * (0.15f + 0.85f * expf(-g * g));
        }
        default:
            return u0_mag;
    }
}

/* Aerodynamic force accumulators (momentum-exchange method).
 * Fx_step/Fy_step are re-zeroed every step; Fx_s/Fy_s are exponential
 * moving averages so the HUD readout is steady. */
static float Fx_step = 0.f, Fy_step = 0.f;
static float Fx_s = 0.f, Fy_s = 0.f;

/* Characteristic length (streamwise obstacle extent, in cells); cached. */
static float char_len = 0.f;
static int   obstacle_dirty = 1;

/* ------------------------------------------------------------------ */
/* Preset shape registry — replayed by lbm_resize() so analytic shapes  */
/* survive grid-resolution changes (re-rasterised, never resampled).    */
/* All coordinates are stored normalised: cx in [0,1] of Nx, cy in      */
/* [0,1] of Ny, scale as a fraction of Ny.                              */
/* ------------------------------------------------------------------ */

#define MAX_SHAPES 128
/* omega: surface speed of the shape boundary (lattice units, CCW positive).
 * Stored as a speed rather than an angular rate so the spin survives grid
 * resizes unchanged (the cylinder radius in cells varies with resolution). */
typedef struct { int id; float cx, cy, scale, angle, omega; } Shape;
static Shape shapes[MAX_SHAPES];
static int   n_shapes = 0;

/* ------------------------------------------------------------------ */
/* Helpers                                                              */
/* ------------------------------------------------------------------ */

static inline int idx(int x, int y) { return y * Nx + x; }

/*
 * (c) Equilibrium distribution / macro-variable relation.
 * The D2Q9 second-order equilibrium is
 *
 *   f_i^eq = w_i * rho * [ 1 + 3 (e_i . u) + 9/2 (e_i . u)^2 - 3/2 |u|^2 ]
 *
 * which reproduces the Navier–Stokes equations to O(u^3) with sound
 * speed c_s^2 = 1/3 and kinematic viscosity nu = (tau - 1/2) / 3.
 */
static inline float feq(int i, float rho, float ux, float uy)
{
    float eu = (float)ex[i] * ux + (float)ey[i] * uy; /* e_i . u        */
    float u2 = ux * ux + uy * uy;                     /* |u|^2          */
    return w9[i] * rho * (1.f + 3.f * eu + 4.5f * eu * eu - 1.5f * u2);
}

static void free_all(void)
{
    for (int i = 0; i < 9; i++) { free(f[i]); f[i] = NULL;
                                  free(f_tmp[i]); f_tmp[i] = NULL; }
    free(obstacle); obstacle = NULL;
    free(rho_f);  rho_f  = NULL;
    free(ux_f);   ux_f   = NULL;
    free(uy_f);   uy_f   = NULL;
    free(speed_f);speed_f= NULL;
    free(vort_f); vort_f = NULL;
    free(schl_f); schl_f = NULL;
    free(dil_f);  dil_f  = NULL;
    free(uwx);    uwx    = NULL;
    free(uwy);    uwy    = NULL;
}

/* Allocate every array for an nx-by-ny grid. Returns 0 on success.
 * On failure everything already allocated is freed again, so a failed
 * (re)allocation can never leak or leave dangling pointers. */
static int alloc_all(int nx, int ny)
{
    size_t n = (size_t)nx * (size_t)ny;
    int ok = 1;
    for (int i = 0; i < 9; i++) {
        f[i]     = (float *)malloc(n * sizeof(float));
        f_tmp[i] = (float *)malloc(n * sizeof(float));
        if (!f[i] || !f_tmp[i]) ok = 0;
    }
    obstacle = (uint8_t *)calloc(n, 1);
    rho_f    = (float *)malloc(n * sizeof(float));
    ux_f     = (float *)malloc(n * sizeof(float));
    uy_f     = (float *)malloc(n * sizeof(float));
    speed_f  = (float *)malloc(n * sizeof(float));
    vort_f   = (float *)malloc(n * sizeof(float));
    schl_f   = (float *)malloc(n * sizeof(float));
    dil_f    = (float *)malloc(n * sizeof(float));
    uwx      = (float *)calloc(n, sizeof(float));
    uwy      = (float *)calloc(n, sizeof(float));
    if (!obstacle || !rho_f || !ux_f || !uy_f || !speed_f || !vort_f
        || !schl_f || !dil_f || !uwx || !uwy) ok = 0;
    if (!ok) { free_all(); return -1; }
    return 0;
}

/* Set every fluid cell to equilibrium at the selected inlet profile
 * (uniform / shear / jet), swept downstream unchanged. */
static void init_distributions(void)
{
    float ca = cosf(angle_rad), sa = sinf(angle_rad);
    ramp_t = 0;                          /* soft-start the inflow */
    float rf = ramp_factor();
    for (int y = 0; y < Ny; y++) {
        float mag = profile_mag(y) * rf;
        float ux0 = mag * ca, uy0 = mag * sa;
        for (int x = 0; x < Nx; x++) {
            int k = idx(x, y);
            float r = 1.f, vx = ux0, vy = uy0;
            if (obstacle[k]) { vx = 0.f; vy = 0.f; }
            for (int i = 0; i < 9; i++) {
                float v = feq(i, r, vx, vy);
                f[i][k] = v; f_tmp[i][k] = v;
            }
            rho_f[k] = r; ux_f[k] = vx; uy_f[k] = vy;
            speed_f[k] = obstacle[k] ? 0.f : mag;
            vort_f[k] = 0.f; schl_f[k] = 0.f; dil_f[k] = 0.f;
        }
    }
    Fx_s = Fy_s = 0.f;
}

/* ------------------------------------------------------------------ */
/* Shape rasterisers (all analytic — no bitmaps)                        */
/* ------------------------------------------------------------------ */

static void raster_cylinder(float cx, float cy, float scale, float us)
{
    /* Circle of diameter scale*Ny centred at (cx*Nx, cy*Ny):
     * solid where (x-xc)^2 + (y-yc)^2 <= r^2.
     * us is the boundary surface speed (CCW positive): every solid cell
     * gets the solid-body rotation velocity u = us/R * (-(y-yc), x-xc),
     * which the moving-wall bounce-back picks up at the surface. */
    float xc = cx * (float)Nx, yc = cy * (float)Ny;
    float r  = 0.5f * scale * (float)Ny;
    float ir = (r > 0.5f) ? us / r : 0.f;
    int x0 = (int)floorf(xc - r) - 1, x1 = (int)ceilf(xc + r) + 1;
    int y0 = (int)floorf(yc - r) - 1, y1 = (int)ceilf(yc + r) + 1;
    for (int y = y0; y <= y1; y++) {
        if (y < 0 || y >= Ny) continue;
        for (int x = x0; x <= x1; x++) {
            if (x < 0 || x >= Nx) continue;
            float dx = (float)x - xc, dy = (float)y - yc;
            if (dx * dx + dy * dy <= r * r) {
                int k = idx(x, y);
                obstacle[k] = OBS_PRESET;
                uwx[k] = -ir * dy;
                uwy[k] =  ir * dx;
            }
        }
    }
}

static void raster_plate(float cx, float cy, float scale, float angle_deg)
{
    /* Thin rectangle: length scale*Ny, thickness max(2, 4% of length),
     * rotated CCW by angle_deg about its centre. A grid cell is solid if
     * its centre, expressed in the plate frame (inverse rotation), lies
     * inside the half-length / half-thickness box. */
    float xc = cx * (float)Nx, yc = cy * (float)Ny;
    float len = scale * (float)Ny;
    float th  = fmaxf(2.f, 0.04f * len);
    float a = angle_deg * (float)M_PI / 180.f;
    float ca = cosf(a), sa = sinf(a);
    float rad = 0.5f * sqrtf(len * len + th * th) + 1.f;
    int x0 = (int)floorf(xc - rad), x1 = (int)ceilf(xc + rad);
    int y0 = (int)floorf(yc - rad), y1 = (int)ceilf(yc + rad);
    for (int y = y0; y <= y1; y++) {
        if (y < 0 || y >= Ny) continue;
        for (int x = x0; x <= x1; x++) {
            if (x < 0 || x >= Nx) continue;
            float dx = (float)x - xc, dy = (float)y - yc;
            float lx =  dx * ca + dy * sa;   /* along plate   */
            float ly = -dx * sa + dy * ca;   /* across plate  */
            if (fabsf(lx) <= 0.5f * len && fabsf(ly) <= 0.5f * th)
                obstacle[idx(x, y)] = OBS_PRESET;
        }
    }
}

/* Standard NACA 4-digit profile, evaluated analytically.
 * s = x/c in [0,1] along the chord. Half-thickness (closed trailing edge):
 *   yt(s) = 5 t (0.2969 sqrt(s) - 0.1260 s - 0.3516 s^2
 *                + 0.2843 s^3 - 0.1036 s^4)
 * Mean camber line:
 *   yc(s) = m/p^2  (2 p s - s^2)            for s <  p
 *   yc(s) = m/(1-p)^2 ((1-2p) + 2 p s - s^2) for s >= p
 * (m = max camber, p = its chordwise position, t = thickness, all /chord). */
static void raster_naca(float m, float p, float t,
                        float cx, float cy, float scale, float angle_deg)
{
    float xc = cx * (float)Nx, yc = cy * (float)Ny;
    float chord = scale * (float)Ny;
    float a = angle_deg * (float)M_PI / 180.f;
    float ca = cosf(a), sa = sinf(a);
    float rad = 0.75f * chord + 2.f;
    int x0 = (int)floorf(xc - rad), x1 = (int)ceilf(xc + rad);
    int y0 = (int)floorf(yc - rad), y1 = (int)ceilf(yc + rad);
    for (int y = y0; y <= y1; y++) {
        if (y < 0 || y >= Ny) continue;
        for (int x = x0; x <= x1; x++) {
            if (x < 0 || x >= Nx) continue;
            float dx = (float)x - xc, dy = (float)y - yc;
            /* Inverse-rotate into the chord-aligned airfoil frame,
             * with the chord spanning lx in [-c/2, +c/2]. */
            float lx =  dx * ca + dy * sa;
            float ly = -dx * sa + dy * ca;
            float s = (lx + 0.5f * chord) / chord;
            if (s < 0.f || s > 1.f) continue;
            float yt = 5.f * t * (0.2969f * sqrtf(s) - 0.1260f * s
                       - 0.3516f * s * s + 0.2843f * s * s * s
                       - 0.1036f * s * s * s * s);
            float ycam = (s < p && p > 0.f)
                       ? (m / (p * p)) * (2.f * p * s - s * s)
                       : (p < 1.f ? (m / ((1.f - p) * (1.f - p)))
                                    * ((1.f - 2.f * p) + 2.f * p * s - s * s)
                                  : 0.f);
            /* Thickness applied vertically about the camber line —
             * the standard small-camber approximation. The 0.6-cell
             * half-thickness floor keeps the thin trailing edge
             * watertight at any grid resolution (cell-centre sampling
             * would otherwise leave gaps where yt < 1 cell); the floor
             * is in grid units, so the rasterised profile converges to
             * the true shape as the resolution rises. */
            if (fabsf(ly - ycam * chord) <= fmaxf(yt * chord, 0.6f))
                obstacle[idx(x, y)] = OBS_PRESET;
        }
    }
}

static void draw_shape(const Shape *sh)
{
    switch (sh->id) {
        case 0: raster_cylinder(sh->cx, sh->cy, sh->scale, sh->omega); break;
        case 1: raster_plate(sh->cx, sh->cy, sh->scale, sh->angle); break;
        case 2: raster_naca(0.00f, 0.0f, 0.12f,
                            sh->cx, sh->cy, sh->scale, sh->angle); break;
        case 3: raster_naca(0.04f, 0.4f, 0.12f,
                            sh->cx, sh->cy, sh->scale, sh->angle); break;
        default: break;
    }
}

/* ------------------------------------------------------------------ */
/* Core LBM kernels                                                     */
/* ------------------------------------------------------------------ */

#ifdef USE_MRT
/* MRT transform matrix and row norms (see the collision comment below).
 * Because M's rows are mutually orthogonal, M^{-1} = M^T D^{-1} with
 * D_k = sum_j M_kj^2 = diag(9,36,36,6,12,6,12,4,4). */
static const float MRT_M[9][9] = {
    { 1, 1, 1, 1, 1, 1, 1, 1, 1},
    {-4,-1,-1,-1,-1, 2, 2, 2, 2},
    { 4,-2,-2,-2,-2, 1, 1, 1, 1},
    { 0, 1, 0,-1, 0, 1,-1,-1, 1},
    { 0,-2, 0, 2, 0, 1,-1,-1, 1},
    { 0, 0, 1, 0,-1, 1, 1,-1,-1},
    { 0, 0,-2, 0, 2, 1, 1,-1,-1},
    { 0, 1,-1, 1,-1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 1,-1, 1,-1}
};
static const float MRT_D[9] = {9, 36, 36, 6, 12, 6, 12, 4, 4};
#endif

/*
 * (a) Per-cell collision. Reads the nine populations of cell k, writes the
 * post-collision values into post[], and returns the local density (the
 * moving-wall bounce-back needs it). The fused collide+stream pass below
 * scatters post[] straight to the neighbours, so unlike the classic
 * two-kernel formulation nothing is ever written back into f[].
 *
 * BGK (default): each population relaxes toward local equilibrium at a
 * single rate omega = 1/tau,
 *
 *   f_i^pc = f_i - (1/tau) * (f_i - f_i^eq(rho, u)),
 *
 * with rho = sum_i f_i, rho*u = sum_i e_i f_i, nu = (tau - 1/2)/3.
 *
 * MRT (-DUSE_MRT, Lallemand & Luo 2000): populations are mapped to moment
 * space, m = M f, each moment is relaxed at its own rate, and the result
 * is mapped back:  f^pc = f - M^{-1} S (m - m^eq).  Moments (rows of M):
 * rho, e (energy), eps (energy^2), jx, qx, jy, qy, pxx, pxy. Only
 * s_pxx = s_pxy = 1/tau set the shear viscosity; the free rates (s_e,
 * s_eps, s_q) damp the non-hydrodynamic "ghost" modes, which is what buys
 * extra stability at high Reynolds compared with BGK.
 */
static inline float collide_cell(int k, float inv_tau, float *restrict post)
{
#ifndef USE_MRT
    float f0=f[0][k], f1=f[1][k], f2=f[2][k], f3=f[3][k], f4=f[4][k],
          f5=f[5][k], f6=f[6][k], f7=f[7][k], f8=f[8][k];
    float rho = f0+f1+f2+f3+f4+f5+f6+f7+f8;
    float inv = 1.f / fmaxf(rho, 0.05f);   /* density floor: NaN guard */
    float ux = (f1 - f3 + f5 - f6 - f7 + f8) * inv;
    float uy = (f2 - f4 + f5 + f6 - f7 - f8) * inv;
    /* Numerical safety: cap |u| (the LBM low-Mach expansion is only
     * valid for u << c_s; runaway u would otherwise produce NaN).
     * Branchless (select), so the surrounding loop stays vectorizable. */
    float u2 = ux*ux + uy*uy;
    float s = (u2 > 0.09f) ? 0.3f / sqrtf(u2) : 1.f;
    ux *= s; uy *= s; u2 = (u2 > 0.09f) ? 0.09f : u2;
    /* Equilibria with shared subexpressions: the axis pair (i, opp(i))
     * differs only in the sign of the linear 3(e.u) term, and the
     * diagonal pairs share (ux+uy) / (ux-uy), so the 9 feq evaluations
     * reduce to ~25 multiplies. */
    float c15 = 1.f - 1.5f * u2;
    float r49 = (4.f/9.f) * rho, r19 = (1.f/9.f) * rho,
          r136 = (1.f/36.f) * rho;
    float sd = ux + uy, dd = ux - uy;
    float feq0 = r49 * c15;
    float axc = r19 * (c15 + 4.5f * ux * ux), axl = r19 * 3.f * ux;
    float ayc = r19 * (c15 + 4.5f * uy * uy), ayl = r19 * 3.f * uy;
    float asc = r136 * (c15 + 4.5f * sd * sd), asl = r136 * 3.f * sd;
    float adc = r136 * (c15 + 4.5f * dd * dd), adl = r136 * 3.f * dd;
    float feq1 = axc + axl, feq3 = axc - axl;
    float feq2 = ayc + ayl, feq4 = ayc - ayl;
    float feq5 = asc + asl, feq7 = asc - asl;
    float feq8 = adc + adl, feq6 = adc - adl;
    /* Smagorinsky closure: the non-equilibrium momentum flux
     *   Pi_ab = sum_i e_ia e_ib (f_i - f_i^eq)
     * is proportional to the strain rate, so the effective
     * relaxation time solves a quadratic with the closed form
     *   tau_eff = 1/2 ( tau + sqrt(tau^2 + 18 Cs^2 |Pi| / rho) ),
     * |Pi| = sqrt(2 Pi:Pi)  (Hou et al. 1996). smag_eff = 0 (LES off)
     * degenerates to tau_eff = tau exactly. */
    float nxx = (f1+f3+f5+f6+f7+f8) - (feq1+feq3+feq5+feq6+feq7+feq8);
    float nyy = (f2+f4+f5+f6+f7+f8) - (feq2+feq4+feq5+feq6+feq7+feq8);
    float nxy = (f5 - f6 + f7 - f8) - (feq5 - feq6 + feq7 - feq8);
    float q = sqrtf(2.f * (nxx*nxx + nyy*nyy + 2.f*nxy*nxy));
    float smag_eff = les_on ? smag_c2 : 0.f;
    float tau_eff = 0.5f * (tau + sqrtf(tau*tau + 18.f*smag_eff*q*inv));
    float omr = 1.f - 1.f / tau_eff;
    (void)inv_tau;
    /* Hermite reconstruction of f^neq (regularized collision) with the
     * trace (acoustic) part damped by BULK_DAMP; see the `regularized`
     * and BULK_DAMP comments at the top of the file. */
    float regf = regularized ? 1.f : 0.f;
    float tdc = BULK_DAMP * (nxx + nyy);
    float dev = 0.5f * (nxx - nyy);
    float g0  = -(2.f/3.f) * tdc;
    float g13 =  0.5f * dev + tdc * (1.f/12.f);
    float g24 = -0.5f * dev + tdc * (1.f/12.f);
    float gp  = tdc * (1.f/12.f) + 0.25f * nxy;
    float gm  = tdc * (1.f/12.f) - 0.25f * nxy;
    float n0 = f0-feq0, n1 = f1-feq1, n2 = f2-feq2, n3 = f3-feq3,
          n4 = f4-feq4, n5 = f5-feq5, n6 = f6-feq6, n7 = f7-feq7,
          n8 = f8-feq8;
    post[0] = feq0 + omr * (n0 + (g0  - n0) * regf);
    post[1] = feq1 + omr * (n1 + (g13 - n1) * regf);
    post[2] = feq2 + omr * (n2 + (g24 - n2) * regf);
    post[3] = feq3 + omr * (n3 + (g13 - n3) * regf);
    post[4] = feq4 + omr * (n4 + (g24 - n4) * regf);
    post[5] = feq5 + omr * (n5 + (gp  - n5) * regf);
    post[6] = feq6 + omr * (n6 + (gm  - n6) * regf);
    post[7] = feq7 + omr * (n7 + (gp  - n7) * regf);
    post[8] = feq8 + omr * (n8 + (gm  - n8) * regf);
    return rho;
#else
    float fv[9], m[9], meq[9], S[9];
    for (int i = 0; i < 9; i++) fv[i] = f[i][k];
    for (int r = 0; r < 9; r++) {
        float acc = 0.f;
        for (int j = 0; j < 9; j++) acc += MRT_M[r][j] * fv[j];
        m[r] = acc;
    }
    float rho = m[0], jx = m[3], jy = m[5];
    float rl = fmaxf(rho, 0.05f);          /* density floor: NaN guard */
    S[0] = 0.f;          /* rho: conserved          */
    S[1] = 1.4f;         /* energy                  */
    S[2] = 1.4f;         /* energy squared          */
    S[3] = 0.f;          /* jx: conserved           */
    S[4] = 1.2f;         /* qx                      */
    S[5] = 0.f;          /* jy: conserved           */
    S[6] = 1.2f;         /* qy                      */
    S[7] = inv_tau;      /* pxx — sets viscosity    */
    S[8] = inv_tau;      /* pxy — sets viscosity    */
    /* Numerical safety, mirroring the BGK branch: cap |u| at 0.3.
     * Momentum moments are normally conserved (S[3] = S[5] = 0), so
     * when the cap engages we also switch those rates to 1 to pull
     * the cell's momentum down to the capped value. */
    float u2 = (jx * jx + jy * jy) / (rl * rl);
    if (u2 > 0.09f) {
        float sc = 0.3f / sqrtf(u2);
        jx *= sc; jy *= sc;
        S[3] = S[5] = 1.f;
    }
    float j2  = (jx * jx + jy * jy) / rl;
    meq[0] = rho;
    meq[1] = -2.f * rho + 3.f * j2;
    meq[2] =  rho       - 3.f * j2;
    meq[3] =  jx;
    meq[4] = -jx;
    meq[5] =  jy;
    meq[6] = -jy;
    meq[7] = (jx * jx - jy * jy) / rl;
    meq[8] = (jx * jy) / rl;
    if (les_on) {
        /* Smagorinsky in moment space: m7, m8 are the traceless
         * stress moments, so their non-equilibrium part gives |Pi|
         * directly (Pi_xx = -Pi_yy = (m7-m7eq)/2, Pi_xy = m8-m8eq).
         * Same effective-tau closed form as the BGK branch. */
        float pxx = m[7] - meq[7], pxy = m[8] - meq[8];
        float q = sqrtf(2.f * (0.5f * pxx * pxx + 2.f * pxy * pxy));
        float tau_eff = 0.5f * (tau + sqrtf(tau*tau + 18.f*smag_c2*q/rl));
        S[7] = S[8] = 1.f / tau_eff;
    }
    for (int j = 0; j < 9; j++) {
        float acc = 0.f;
        for (int r = 0; r < 9; r++)
            acc += (MRT_M[r][j] / MRT_D[r]) * S[r] * (m[r] - meq[r]);
        post[j] = fv[j] - acc;
    }
    return rho;
#endif
}

/*
 * (e) Moving-wall halfway bounce-back of population i pushed from fluid
 * cell k into solid cell kt, plus the momentum-exchange force (Ladd 1994).
 * For a wall moving with velocity u_w the reflected population is
 *
 *   f_opp(x, t+1) = f_i^pc(x) - 6 w_i rho (e_i . u_w),
 *
 * and the momentum the wall absorbs over the link is e_i (f_i^pc + f_opp).
 * A resting wall (u_w = 0) reduces to the classic f_opp = f_i^pc with
 * dF = 2 e_i f_i^pc. u_w is nonzero only on rotating cylinder cells.
 */
static inline void bounce_into_solid(int k, int kt, int i, float rho,
                                     float fi, float *fx, float *fy)
{
    float du = (float)ex[i] * uwx[kt] + (float)ey[i] * uwy[kt];
    float fr = fi - 6.f * w9[i] * rho * du;
    f_tmp[opp[i]][k] = fr;
    *fx += (fi + fr) * (float)ex[i];
    *fy += (fi + fr) * (float)ey[i];
}

/* Fixup over rows [y0, y1): convert the populations pushed into solid
 * cells into bounce-back + momentum-exchange force (accumulated through
 * fx and fy so concurrent slices don't share an accumulator). See the
 * call site in collide_stream for the full story. */
static void fixup_rows(int y0, int y1, float *fx, float *fy)
{
    for (int y = y0; y < y1; y++) {
        for (int x = 0; x < Nx; x++) {
            int kt = idx(x, y);
            if (!obstacle[kt]) continue;
            for (int j = 1; j < 9; j++) {
                int nx2 = x + ex[j], ny2 = y + ey[j];
                if (nx2 < 0 || nx2 >= Nx || ny2 < 0 || ny2 >= Ny) continue;
                int n = idx(nx2, ny2);
                if (obstacle[n]) continue;
                int i = opp[j];               /* direction n -> solid */
                /* rho of n: collision conserves it, so the pre-collision
                 * sum equals the value collide_cell saw. */
                float rho = f[0][n]+f[1][n]+f[2][n]+f[3][n]+f[4][n]
                          + f[5][n]+f[6][n]+f[7][n]+f[8][n];
                bounce_into_solid(n, kt, i, rho, f_tmp[i][kt], fx, fy);
            }
        }
    }
}

/* Push the post-collision populations of a domain-border fluid cell.
 * Only ~2(Nx+Ny) cells take this path, so full per-link branching is
 * cheap here; the interior uses the unrolled fast path in collide_stream.
 */
static void push_border_cell(int x, int y, const float *post)
{
    int k = idx(x, y);
    for (int i = 1; i < 9; i++) {
        int tx = x + ex[i], ty = y + ey[i];
        if (ty < 0 || ty >= Ny) {
            /* Top/bottom wall. No-slip: halfway bounce-back, the
             * population returns reversed into its source cell.
             * Free-slip: specular reflection — tangential momentum is
             * preserved, normal momentum flips, so the reflected
             * population lands one tangential step over (the push-scheme
             * mirror image of the old pull-scheme wall handling). */
            if (!wall_noslip && tx >= 0 && tx < Nx && !obstacle[idx(tx, y)])
                f_tmp[mirY[i]][idx(tx, y)] = post[i];
            else
                f_tmp[opp[i]][k] = post[i];
        } else if (tx < 0 || tx >= Nx) {
            /* Leaves the domain through the inlet/outlet column; the
             * Zou-He / zero-gradient pass overwrites these slots. */
        } else {
            /* Plain push — also into solid cells: the fixup pass in
             * collide_stream reads the pushed value back out of the
             * solid and converts it into the bounce-back + force. */
            f_tmp[i][idx(tx, ty)] = post[i];
        }
    }
}

static void sweep_interior(int y0, int y1);

/*
 * Worker pool. The interior sweep is row-partitioned across threads:
 * every f_tmp slot has exactly one producer cell, so concurrent slices
 * write disjoint locations and read only the immutable f arrays — no
 * locks needed, just a start barrier (releases the workers) and an end
 * barrier (everyone finished) per step. Workers block on the barrier
 * between steps (futex wait, no spinning). Threads require the page to
 * be crossOriginIsolated (SharedArrayBuffer); the single-thread build
 * compiles the no-op fallback below.
 */
#ifdef __EMSCRIPTEN_PTHREADS__
#include <pthread.h>
#define MAX_THREADS 16
static int n_threads = 1;            /* total participants incl. main  */
static int pool_started = 0;
static pthread_barrier_t bar_start, bar_mid, bar_end;

/* Per-thread force accumulators, padded to separate cache lines. The
 * partials are reduced on the main thread in tid order, so the result
 * is deterministic for a given thread count. */
typedef struct { float fx, fy; char pad[56]; } ForcePart;
static ForcePart fpart[MAX_THREADS];

static void sweep_slice(int tid)
{
    int rows = Ny - 2;
    if (rows <= 0) return;
    int chunk = (rows + n_threads - 1) / n_threads;
    int y0 = 1 + tid * chunk;
    int y1 = y0 + chunk;
    if (y1 > Ny - 1) y1 = Ny - 1;
    if (y0 < y1) sweep_interior(y0, y1);
}

static void fixup_slice(int tid)
{
    int chunk = (Ny + n_threads - 1) / n_threads;
    int y0 = tid * chunk;
    int y1 = y0 + chunk;
    if (y1 > Ny) y1 = Ny;
    fpart[tid].fx = 0.f; fpart[tid].fy = 0.f;
    if (y0 < y1) fixup_rows(y0, y1, &fpart[tid].fx, &fpart[tid].fy);
}

static volatile int pool_go = 0;

static void *worker_main(void *arg)
{
    int tid = (int)(intptr_t)arg;
    /* Spin until lbm_set_threads has counted the spawned workers and
     * sized the barriers (a one-time, microsecond window). */
    while (!__atomic_load_n(&pool_go, __ATOMIC_ACQUIRE)) {}
    for (;;) {
        pthread_barrier_wait(&bar_start);
        sweep_slice(tid);
        /* The fixup reads values the interior sweep pushed into solid
         * cells from *other* slices, so it needs the full sweep done. */
        pthread_barrier_wait(&bar_mid);
        fixup_slice(tid);
        pthread_barrier_wait(&bar_end);
    }
    return NULL;
}

/* Start the worker pool (one-shot: call once, before heavy stepping).
 * Threads must come from the preallocated Emscripten pool
 * (PTHREAD_POOL_SIZE), otherwise they cannot start until the main
 * thread yields to the event loop. The barriers are sized by the count
 * that actually spawned, and the workers are only released after, so a
 * partial spawn can never deadlock the step barrier. */
EXPORT void lbm_set_threads(int n)
{
    if (pool_started || n < 2) return;
    if (n > MAX_THREADS) n = MAX_THREADS;
    int spawned = 1;                 /* main thread is participant 0 */
    for (int t = 1; t < n; t++) {
        pthread_t th;
        if (pthread_create(&th, NULL, worker_main,
                           (void *)(intptr_t)spawned) == 0)
            spawned++;
    }
    n_threads = spawned;
    pool_started = 1;
    if (spawned < 2) return;
    pthread_barrier_init(&bar_start, NULL, spawned);
    pthread_barrier_init(&bar_mid, NULL, spawned);
    pthread_barrier_init(&bar_end, NULL, spawned);
    __atomic_store_n(&pool_go, 1, __ATOMIC_RELEASE);
}

EXPORT int lbm_get_threads(void) { return n_threads; }

/* Run interior sweep + solid-link fixup across the pool; the per-thread
 * force partials are reduced here in tid order (deterministic). */
static void sweep_dispatch(void)
{
    if (n_threads > 1) {
        pthread_barrier_wait(&bar_start);   /* release the workers */
        sweep_slice(0);                     /* main takes slice 0  */
        pthread_barrier_wait(&bar_mid);
        fixup_slice(0);
        pthread_barrier_wait(&bar_end);
        for (int t = 0; t < n_threads; t++) {
            Fx_step += fpart[t].fx;
            Fy_step += fpart[t].fy;
        }
    } else {
        sweep_interior(1, Ny - 1);
        fixup_rows(0, Ny, &Fx_step, &Fy_step);
    }
}
#else
EXPORT void lbm_set_threads(int n) { (void)n; }
EXPORT int  lbm_get_threads(void)  { return 1; }
static void sweep_dispatch(void)
{
    sweep_interior(1, Ny - 1);
    fixup_rows(0, Ny, &Fx_step, &Fy_step);
}
#endif

static void collide_stream(void)
{
    /*
     * (a)+(b) fused: collide each fluid cell, then push ("scatter") the
     * post-collision populations one lattice link outward in the same
     * pass:
     *
     *   f_i(x + e_i, t+1) = f_i^pc(x, t)
     *
     * Fusing the two kernels halves the memory traffic per step — one
     * read and one write of all nine arrays instead of two of each —
     * and removes the old full-grid obstacle repair pass entirely: a
     * push that lands in a solid cell is reflected back into its source
     * on the spot (bounce_into_solid), which is also where the momentum-
     * exchange force and the moving-wall (rotation) term live. For each
     * direction the write stream f_tmp[i][k + shift_i] advances
     * sequentially with k, so all nine stores stay cache-friendly.
     * Solid cells are skipped entirely: nothing ever reads their
     * populations (their f_tmp entries simply keep stale finite values).
     */
    Fx_step = 0.f; Fy_step = 0.f;
    float inv_tau = 1.f / tau;
    float bpost[9];

    /* Domain border first (rows y = 0 / Ny-1, columns x = 0 / Nx-1), so
     * the interior sweep below is a single self-contained block whose
     * restrict-qualified pointers never overlap other writers. */
    for (int x = 0; x < Nx; x++) {
        int kb = x, kt2 = (Ny - 1) * Nx + x;
        if (!obstacle[kb]) {
            collide_cell(kb, inv_tau, bpost);
            f_tmp[0][kb] = bpost[0];
            push_border_cell(x, 0, bpost);
        }
        if (!obstacle[kt2]) {
            collide_cell(kt2, inv_tau, bpost);
            f_tmp[0][kt2] = bpost[0];
            push_border_cell(x, Ny - 1, bpost);
        }
    }
    for (int y = 1; y < Ny - 1; y++) {
        int k = y * Nx, ke = k + Nx - 1;
        if (!obstacle[k]) {
            collide_cell(k, inv_tau, bpost);
            f_tmp[0][k] = bpost[0];
            push_border_cell(0, y, bpost);
        }
        if (!obstacle[ke]) {
            collide_cell(ke, inv_tau, bpost);
            f_tmp[0][ke] = bpost[0];
            push_border_cell(Nx - 1, y, bpost);
        }
    }

    /*
     * Interior sweep + fixup pass, row-partitioned across the thread
     * pool when one is running. The fixup converts every population
     * pushed into a solid cell during the sweep (each fluid neighbour n
     * pushed f_i^pc into the solid's f_tmp slot, i = direction n->solid)
     * into halfway bounce-back with the moving-wall term, accumulating
     * the momentum-exchange force, and overwrites the stray values solid
     * cells pushed into their fluid neighbours. For fluid-dominated
     * grids the scan reads only the byte mask — same cost as the repair
     * pass the old two-kernel formulation needed anyway.
     */
    sweep_dispatch();

    /* Exponential moving average for a readable HUD (~40-step memory). */
    Fx_s += 0.05f * (Fx_step - Fx_s);
    Fy_s += 0.05f * (Fy_step - Fy_s);
}

/*
 * Interior sweep over rows [y0, y1). Completely branchless: every push
 * target is in-grid, and solid cells are collided and pushed like fluid
 * (their populations sit at the w9 equilibrium fixed point or hold
 * stale-but-finite values, so the arithmetic is harmless and the values
 * they leak into neighbours are overwritten by the fixup pass).
 * The BGK body is inlined here with restrict-qualified pointer copies:
 * without them LLVM cannot prove the nine source and nine destination
 * arrays are disjoint and refuses to vectorize. The LES toggle is folded
 * into the arithmetic (smag_eff = 0 makes tau_eff = tau exactly), so the
 * loop body has no branches at all.
 *
 * Thread-safety: every f_tmp slot has exactly one producer cell, so
 * row-partitioned concurrent sweeps write disjoint locations and read
 * only the immutable f arrays — race-free with no locks.
 */
static void sweep_interior(int y0, int y1)
{
    float inv_tau = 1.f / tau;
#ifndef USE_MRT
    {
        const float *restrict p0 = f[0], *restrict p1 = f[1],
                    *restrict p2 = f[2], *restrict p3 = f[3],
                    *restrict p4 = f[4], *restrict p5 = f[5],
                    *restrict p6 = f[6], *restrict p7 = f[7],
                    *restrict p8 = f[8];
        float *restrict t0 = f_tmp[0], *restrict t1 = f_tmp[1],
              *restrict t2 = f_tmp[2], *restrict t3 = f_tmp[3],
              *restrict t4 = f_tmp[4], *restrict t5 = f_tmp[5],
              *restrict t6 = f_tmp[6], *restrict t7 = f_tmp[7],
              *restrict t8 = f_tmp[8];
        float smag_eff = les_on ? smag_c2 : 0.f;
        float regf = regularized ? 1.f : 0.f;
        float tau_l = tau;
        (void)inv_tau;
        for (int y = y0; y < y1; y++) {
            int k = y * Nx + 1;
            for (int x = 1; x < Nx - 1; x++, k++) {
                float f0=p0[k], f1=p1[k], f2=p2[k], f3=p3[k], f4=p4[k],
                      f5=p5[k], f6=p6[k], f7=p7[k], f8=p8[k];
                float rho = f0+f1+f2+f3+f4+f5+f6+f7+f8;
                /* Density floor: only engages in already-broken states,
                 * where it stops 1/rho from seeding a NaN cascade. */
                float inv = 1.f / fmaxf(rho, 0.05f);
                float ux = (f1 - f3 + f5 - f6 - f7 + f8) * inv;
                float uy = (f2 - f4 + f5 + f6 - f7 - f8) * inv;
                float u2 = ux*ux + uy*uy;
                float s = (u2 > 0.09f) ? 0.3f / sqrtf(u2) : 1.f;
                ux *= s; uy *= s; u2 = (u2 > 0.09f) ? 0.09f : u2;
                float c15 = 1.f - 1.5f * u2;
                float r49 = (4.f/9.f) * rho, r19 = (1.f/9.f) * rho,
                      r136 = (1.f/36.f) * rho;
                float sd = ux + uy, dd = ux - uy;
                float feq0 = r49 * c15;
                float axc = r19 * (c15 + 4.5f * ux * ux), axl = r19 * 3.f * ux;
                float ayc = r19 * (c15 + 4.5f * uy * uy), ayl = r19 * 3.f * uy;
                float asc = r136 * (c15 + 4.5f * sd * sd), asl = r136 * 3.f * sd;
                float adc = r136 * (c15 + 4.5f * dd * dd), adl = r136 * 3.f * dd;
                float feq1 = axc + axl, feq3 = axc - axl;
                float feq2 = ayc + ayl, feq4 = ayc - ayl;
                float feq5 = asc + asl, feq7 = asc - asl;
                float feq8 = adc + adl, feq6 = adc - adl;
                float nxx = (f1+f3+f5+f6+f7+f8) - (feq1+feq3+feq5+feq6+feq7+feq8);
                float nyy = (f2+f4+f5+f6+f7+f8) - (feq2+feq4+feq5+feq6+feq7+feq8);
                float nxy = (f5 - f6 + f7 - f8) - (feq5 - feq6 + feq7 - feq8);
                float q = sqrtf(2.f * (nxx*nxx + nyy*nyy + 2.f*nxy*nxy));
                float tau_eff = 0.5f * (tau_l
                              + sqrtf(tau_l*tau_l + 18.f*smag_eff*q*inv));
                float omr = 1.f - 1.f / tau_eff;
                /* Hermite reconstruction of f^neq from its stress moments
                 * (regularized collision); regf selects it branchlessly.
                 * The trace (acoustic) part is damped by BULK_DAMP. */
                float tdc = BULK_DAMP * (nxx + nyy);
                float dev = 0.5f * (nxx - nyy);
                float g0  = -(2.f/3.f) * tdc;
                float g13 =  0.5f * dev + tdc * (1.f/12.f);
                float g24 = -0.5f * dev + tdc * (1.f/12.f);
                float gp  = tdc * (1.f/12.f) + 0.25f * nxy;
                float gm  = tdc * (1.f/12.f) - 0.25f * nxy;
                float n0 = f0-feq0, n1 = f1-feq1, n2 = f2-feq2, n3 = f3-feq3,
                      n4 = f4-feq4, n5 = f5-feq5, n6 = f6-feq6, n7 = f7-feq7,
                      n8 = f8-feq8;
                t0[k]          = feq0 + omr * (n0 + (g0  - n0) * regf);
                t1[k + 1]      = feq1 + omr * (n1 + (g13 - n1) * regf);
                t2[k + Nx]     = feq2 + omr * (n2 + (g24 - n2) * regf);
                t3[k - 1]      = feq3 + omr * (n3 + (g13 - n3) * regf);
                t4[k - Nx]     = feq4 + omr * (n4 + (g24 - n4) * regf);
                t5[k + Nx + 1] = feq5 + omr * (n5 + (gp  - n5) * regf);
                t6[k + Nx - 1] = feq6 + omr * (n6 + (gm  - n6) * regf);
                t7[k - Nx - 1] = feq7 + omr * (n7 + (gp  - n7) * regf);
                t8[k - Nx + 1] = feq8 + omr * (n8 + (gm  - n8) * regf);
            }
        }
    }
#else
    for (int y = y0; y < y1; y++) {
        int k = y * Nx + 1;
        for (int x = 1; x < Nx - 1; x++, k++) {
            float post[9];
            collide_cell(k, inv_tau, post);
            int kn = k + Nx, ks = k - Nx;
            f_tmp[0][k]      = post[0];
            f_tmp[1][k + 1]  = post[1];
            f_tmp[2][kn]     = post[2];
            f_tmp[3][k - 1]  = post[3];
            f_tmp[4][ks]     = post[4];
            f_tmp[5][kn + 1] = post[5];
            f_tmp[6][kn - 1] = post[6];
            f_tmp[7][ks - 1] = post[7];
            f_tmp[8][ks + 1] = post[8];
        }
    }
#endif
}

static void apply_boundaries(void)
{
    /* Inlet (x = 0): Zou-He velocity boundary. After streaming, the
     * populations pointing into the domain (f1, f5, f8) are unknown.
     * Imposing rho*u = sum e_i f_i with the prescribed u and using the
     * bounce-back of the non-equilibrium part normal to the boundary
     * gives the classic closed-form solution (Zou & He 1997):
     *
     *   rho = (f0 + f2 + f4 + 2 (f3 + f6 + f7)) / (1 - ux)
     *   f1  = f3 + 2/3 rho ux
     *   f5  = f7 - 1/2 (f2 - f4) + 1/6 rho ux + 1/2 rho uy
     *   f8  = f6 + 1/2 (f2 - f4) + 1/6 rho ux - 1/2 rho uy
     */
    float ca = cosf(angle_rad), sa = sinf(angle_rad);
    float rf = ramp_factor();                /* soft start after (re)init */
    for (int y = 0; y < Ny; y++) {
        int k = idx(0, y);
        if (obstacle[k]) continue;
        float mag = profile_mag(y) * rf;     /* uniform / shear / jet */
        float ux0 = mag * ca, uy0 = mag * sa;
        float f0=f_tmp[0][k], f2=f_tmp[2][k], f3=f_tmp[3][k],
              f4=f_tmp[4][k], f6=f_tmp[6][k], f7=f_tmp[7][k];
        float rho = (f0 + f2 + f4 + 2.f * (f3 + f6 + f7)) / (1.f - ux0);
        f_tmp[1][k] = f3 + (2.f/3.f) * rho * ux0;
        f_tmp[5][k] = f7 - 0.5f * (f2 - f4) + (1.f/6.f) * rho * ux0
                         + 0.5f * rho * uy0;
        f_tmp[8][k] = f6 + 0.5f * (f2 - f4) + (1.f/6.f) * rho * ux0
                         - 0.5f * rho * uy0;
    }
    /* Outlet (x = Nx-1): zero-gradient extrapolation for the populations
     * pointing back into the domain (f3, f6, f7). First-order, robust,
     * and non-reflecting enough for a wind-tunnel demo. */
    for (int y = 0; y < Ny; y++) {
        int k = idx(Nx - 1, y), km = idx(Nx - 2, y);
        if (obstacle[k]) continue;
        f_tmp[3][k] = f_tmp[3][km];
        f_tmp[6][k] = f_tmp[6][km];
        f_tmp[7][k] = f_tmp[7][km];
    }
}

/*
 * Absorbing (sponge) layers. Over the columns nearest the inlet and the
 * outlet, relax f_tmp toward the freestream equilibrium with a strength
 * that grows quadratically toward the boundary. Outgoing vortices and —
 * more importantly — the spurious short-wavelength acoustic waves the
 * lattice supports are absorbed before the open boundaries can reflect
 * them back into the domain (the Zou-He inlet is acoustically rigid).
 * The relaxation target follows the selected inlet profile, so shear
 * layers and jets pass through the inlet sponge unharmed.
 */
static void apply_sponge(void)
{
    int W = Nx / 25 + 6;
    if (W > 32) W = 32;          /* ~30 cells absorbs everything relevant */
    if (W > Nx / 3) W = Nx / 3;
    float ca = cosf(angle_rad), sa = sinf(angle_rad);
    float rf = ramp_factor();
    /* Row-major sweep (cache-friendly), with the 9-value equilibrium
     * target hoisted per row — it only depends on y via the profile. */
    for (int y = 0; y < Ny; y++) {
        float mag = profile_mag(y) * rf;
        float ux0 = mag * ca, uy0 = mag * sa;
        float ft[9];
        for (int i = 0; i < 9; i++) ft[i] = feq(i, 1.f, ux0, uy0);
        int row = y * Nx;
        for (int d = 0; d < W; d++) {
            float t = (float)(W - d) / (float)W;
            float sig = 0.25f * t * t;
            int kl = row + d, kr = row + Nx - 1 - d;
            if (!obstacle[kl])
                for (int i = 0; i < 9; i++)
                    f_tmp[i][kl] += sig * (ft[i] - f_tmp[i][kl]);
            if (!obstacle[kr])
                for (int i = 0; i < 9; i++)
                    f_tmp[i][kr] += sig * (ft[i] - f_tmp[i][kr]);
        }
    }
}

/* ------------------------------------------------------------------ */
/* Public (exported) API                                                */
/* ------------------------------------------------------------------ */

EXPORT int lbm_init(int nx, int ny)
{
    if (nx < 8 || ny < 8) return -1;
    free_all();
    n_shapes = 0;
    Nx = nx; Ny = ny;
    if (alloc_all(nx, ny) != 0) { Nx = Ny = 0; return -1; }
    obstacle_dirty = 1;
    init_distributions();
    return 0;
}

EXPORT void lbm_set_params(float u0, float angle_deg, float tau_in)
{
    u0_mag = u0;
    if (u0_mag < 0.f)    u0_mag = 0.f;
    if (u0_mag > 0.3f)   u0_mag = 0.3f;   /* low-Mach validity bound  */
    angle_rad = angle_deg * (float)M_PI / 180.f;
    tau = tau_in;
    /* nu > 0 stability floor. 0.501 gives nu ~ 3.3e-4: bare BGK is
     * unstable this close to 0.5, but the Smagorinsky eddy viscosity
     * raises the effective tau exactly where gradients steepen, which
     * is what makes Re ~ 10^4 wind-tunnel runs possible. */
    if (tau < 0.501f) tau = 0.501f;
}

EXPORT void lbm_set_wall_mode(int noslip) { wall_noslip = noslip ? 1 : 0; }

/* Toggle the Smagorinsky subgrid model; cs <= 0 keeps the current
 * constant (default Cs = 0.18). */
EXPORT void lbm_set_les(int on, float cs)
{
    les_on = on ? 1 : 0;
    if (cs > 0.f) smag_c2 = cs * cs;
}

/* Toggle the Hermite-regularized collision (BGK only; MRT damps its
 * ghost modes through its own relaxation rates). */
EXPORT void lbm_set_regularized(int on) { regularized = on ? 1 : 0; }

/* Select which derivative fields lbm_compute_fields() fills (bitmask:
 * 1 vorticity, 2 schlieren, 4 dilatation); skipping unused ones saves
 * a finite-difference pass per frame. */
EXPORT void lbm_set_active_field(int mask) { field_mask = mask & 7; }

/* Select the inlet/initial flow profile: 0 uniform, 1 shear layer, 2 jet. */
EXPORT void lbm_set_flow_profile(int p)
{
    flow_profile = (p < 0) ? 0 : (p > 2 ? 2 : p);
}

/* Re-initialize the flow to the current profile, keeping obstacles. */
EXPORT void lbm_reset_flow(void)
{
    if (!Nx) return;
    init_distributions();
}

EXPORT void lbm_step(void)
{
    if (!Nx) return;
    if (ramp_t < RAMP_STEPS) ramp_t++;
    collide_stream();
    apply_boundaries();
    apply_sponge();
    for (int i = 0; i < 9; i++) {        /* ping-pong swap */
        float *t = f[i]; f[i] = f_tmp[i]; f_tmp[i] = t;
    }
}

EXPORT void lbm_compute_fields(void)
{
    if (!Nx) return;
    size_t n = (size_t)Nx * (size_t)Ny;
    /*
     * (c) Macro-variable extraction. The hydrodynamic fields are the
     * low-order velocity moments of the distributions:
     *
     *   rho   = sum_i f_i
     *   rho u = sum_i e_i f_i
     *
     * Pressure follows from the ideal lattice gas: p = c_s^2 rho = rho/3.
     */
    for (size_t k = 0; k < n; k++) {
        if (obstacle[k]) {
            rho_f[k] = 1.f; ux_f[k] = 0.f; uy_f[k] = 0.f; speed_f[k] = 0.f;
            continue;
        }
        float f0=f[0][k], f1=f[1][k], f2=f[2][k], f3=f[3][k], f4=f[4][k],
              f5=f[5][k], f6=f[6][k], f7=f[7][k], f8=f[8][k];
        float rho = f0+f1+f2+f3+f4+f5+f6+f7+f8;
        float inv = 1.f / rho;
        float ux = (f1 - f3 + f5 - f6 - f7 + f8) * inv;
        float uy = (f2 - f4 + f5 + f6 - f7 - f8) * inv;
        rho_f[k] = rho; ux_f[k] = ux; uy_f[k] = uy;
        speed_f[k] = sqrtf(ux * ux + uy * uy);
    }
    /*
     * (d) Derivative fields by finite differences: second-order central
     * stencil in the interior, falling back to first-order one-sided
     * differences on the domain border (denominator 1 instead of 2).
     *
     *   vorticity   omega = d(uy)/dx - d(ux)/dy   (rotation of the flow)
     *   schlieren   |grad rho|                    (numerical schlieren —
     *               wave fronts and shear layers, like the experimental
     *               knife-edge technique)
     *   dilatation  div u = d(ux)/dx + d(uy)/dy   (zero for incompressible
     *               motion, so it isolates the acoustic radiation the
     *               weakly-compressible LBM actually carries)
     */
    int wv = field_mask & FIELD_VORT, ws = field_mask & FIELD_SCHL,
        wd = field_mask & FIELD_DIL;
    if (!field_mask) return;
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            int k = idx(x, y);
            if (obstacle[k]) {
                vort_f[k] = 0.f; schl_f[k] = 0.f; dil_f[k] = 0.f;
                continue;
            }
            int xm = x > 0 ? x - 1 : x, xp = x < Nx - 1 ? x + 1 : x;
            int ym = y > 0 ? y - 1 : y, yp = y < Ny - 1 ? y + 1 : y;
            float ddx = 1.f / (float)(xp - xm), ddy = 1.f / (float)(yp - ym);
            if (wv) {
                float duy_dx = (uy_f[idx(xp, y)] - uy_f[idx(xm, y)]) * ddx;
                float dux_dy = (ux_f[idx(x, yp)] - ux_f[idx(x, ym)]) * ddy;
                vort_f[k] = duy_dx - dux_dy;
            }
            if (ws) {
                float drho_dx = (rho_f[idx(xp, y)] - rho_f[idx(xm, y)]) * ddx;
                float drho_dy = (rho_f[idx(x, yp)] - rho_f[idx(x, ym)]) * ddy;
                schl_f[k] = sqrtf(drho_dx * drho_dx + drho_dy * drho_dy);
            }
            if (wd) {
                float dux_dx = (ux_f[idx(xp, y)] - ux_f[idx(xm, y)]) * ddx;
                float duy_dy = (uy_f[idx(x, yp)] - uy_f[idx(x, ym)]) * ddy;
                dil_f[k]  = dux_dx + duy_dy;
            }
        }
    }
}

EXPORT void lbm_set_obstacle(int x, int y, int solid)
{
    if (!Nx || x < 0 || x >= Nx || y < 0 || y >= Ny) return;
    int k = idx(x, y);
    uint8_t v = solid <= 0 ? OBS_FLUID
              : (solid >= OBS_PRESET ? OBS_PRESET : OBS_PAINTED);
    if (obstacle[k] && v == OBS_FLUID) {
        /* Cell re-opened to fluid: refill with rest-state equilibrium so
         * the surrounding flow relaxes into it without a pressure shock. */
        for (int i = 0; i < 9; i++) { f[i][k] = w9[i]; f_tmp[i][k] = w9[i]; }
    }
    obstacle[k] = v;
    uwx[k] = uwy[k] = 0.f;   /* painted / erased cells are resting walls */
    obstacle_dirty = 1;
}

EXPORT void lbm_add_preset(int shape_id, float cx, float cy,
                           float scale, float angle)
{
    if (!Nx || shape_id < 0 || shape_id > 3) return;
    Shape sh = { shape_id, cx, cy, scale, angle, 0.f };
    if (n_shapes < MAX_SHAPES) shapes[n_shapes++] = sh;
    draw_shape(&sh);
    obstacle_dirty = 1;
}

/* Set the boundary surface speed of an already-placed preset (lattice
 * units, CCW positive). Only cylinders rotate; the index is the placement
 * order, i.e. the same order lbm_add_preset was called in. Re-rasterising
 * the shape rewrites its wall-velocity footprint at the current grid. */
EXPORT void lbm_set_shape_omega(int index, float us)
{
    if (!Nx || index < 0 || index >= n_shapes) return;
    if (us >  0.3f) us =  0.3f;     /* same low-Mach bound as the flow */
    if (us < -0.3f) us = -0.3f;
    shapes[index].omega = us;
    if (shapes[index].id == 0) draw_shape(&shapes[index]);
}

EXPORT void lbm_clear_obstacles(void)
{
    if (!Nx) return;
    size_t n = (size_t)Nx * (size_t)Ny;
    for (size_t k = 0; k < n; k++) {
        if (obstacle[k]) {
            obstacle[k] = OBS_FLUID;
            uwx[k] = uwy[k] = 0.f;
            for (int i = 0; i < 9; i++) {
                f[i][k] = w9[i]; f_tmp[i][k] = w9[i];
            }
        }
    }
    n_shapes = 0;
    obstacle_dirty = 1;
}

EXPORT float lbm_get_char_length(void)
{
    /* Characteristic length = streamwise (x) extent of the obstacle's
     * bounding box: the chord for airfoils/plates, the diameter for a
     * cylinder. Falls back to the channel height when the tunnel is
     * empty. Cached until the obstacle mask changes. */
    if (!Nx) return 0.f;
    if (obstacle_dirty) {
        int minx = Nx, maxx = -1;
        for (int y = 0; y < Ny; y++)
            for (int x = 0; x < Nx; x++)
                if (obstacle[idx(x, y)]) {
                    if (x < minx) minx = x;
                    if (x > maxx) maxx = x;
                }
        char_len = (maxx >= minx) ? (float)(maxx - minx + 1) : (float)Ny;
        obstacle_dirty = 0;
    }
    return char_len;
}

EXPORT void lbm_get_forces(float *Cl, float *Cd)
{
    /* Normalise the momentum-exchange force by the dynamic pressure
     * q = 1/2 rho0 u0^2 L (rho0 = 1 in lattice units), then project on
     * the freestream direction d = (cos th, sin th) for drag and its
     * perpendicular l = (-sin th, cos th) for lift. */
    float L = lbm_get_char_length();
    float q = 0.5f * u0_mag * u0_mag * L;
    if (q < 1e-12f) { if (Cl) *Cl = 0.f; if (Cd) *Cd = 0.f; return; }
    float ca = cosf(angle_rad), sa = sinf(angle_rad);
    float drag = ( Fx_s * ca + Fy_s * sa) / q;
    float lift = (-Fx_s * sa + Fy_s * ca) / q;
    if (Cd) *Cd = drag;
    if (Cl) *Cl = lift;
}

EXPORT int lbm_resize(int nx, int ny)
{
    if (nx < 8 || ny < 8 || !Nx) return -1;
    int oldNx = Nx, oldNy = Ny;
    uint8_t *old_obs = obstacle;   /* keep the mask alive for resampling */
    obstacle = NULL;

    /* Free every flow array, then allocate at the new resolution.
     * alloc_all() cleans up after itself on failure, so the only
     * remaining allocation to release on error is the saved mask. */
    free_all();
    Nx = nx; Ny = ny;
    if (alloc_all(nx, ny) != 0) { free(old_obs); Nx = Ny = 0; return -1; }

    /* Persistent obstacles, part 1: hand-painted cells are resampled
     * nearest-neighbour from the old mask. */
    if (old_obs) {
        for (int y = 0; y < ny; y++) {
            int sy = (int)(((long long)y * oldNy) / ny);
            for (int x = 0; x < nx; x++) {
                int sx = (int)(((long long)x * oldNx) / nx);
                if (old_obs[sy * oldNx + sx] == OBS_PAINTED)
                    obstacle[idx(x, y)] = OBS_PAINTED;
            }
        }
        free(old_obs);
    }
    /* Persistent obstacles, part 2: preset shapes are re-rasterised
     * analytically from the registry — crisp at any resolution. */
    for (int s = 0; s < n_shapes; s++) draw_shape(&shapes[s]);

    obstacle_dirty = 1;
    init_distributions();
    return 0;
}

/* Zero-copy pointer accessors: JS wraps these addresses in TypedArray
 * views over WebAssembly.Memory and uploads them straight to WebGL. */
EXPORT float   *lbm_get_ux_ptr(void)        { return ux_f; }
EXPORT float   *lbm_get_uy_ptr(void)        { return uy_f; }
EXPORT float   *lbm_get_rho_ptr(void)       { return rho_f; }
EXPORT float   *lbm_get_speed_ptr(void)     { return speed_f; }
EXPORT float   *lbm_get_vorticity_ptr(void) { return vort_f; }
EXPORT float   *lbm_get_schlieren_ptr(void) { return schl_f; }
EXPORT float   *lbm_get_dilatation_ptr(void){ return dil_f; }
EXPORT uint8_t *lbm_get_obstacle_ptr(void)  { return obstacle; }

#ifdef LBM_TEST_MAIN
/* Native/Node smoke-test harness — not part of the shipped wasm build.
 * Compile with -DLBM_TEST_MAIN; see README. */
#include <stdio.h>
#ifdef __EMSCRIPTEN__
#include <emscripten/heap.h>
#endif
static int field_finite(void)
{
    size_t n = (size_t)Nx * (size_t)Ny;
    for (size_t k = 0; k < n; k++)
        if (!isfinite(rho_f[k]) || !isfinite(ux_f[k]) || !isfinite(uy_f[k])
            || !isfinite(vort_f[k]) || !isfinite(schl_f[k])
            || !isfinite(dil_f[k]))
            return 0;
    return 1;
}
int main(void)
{
    int fails = 0;
#ifdef __EMSCRIPTEN_PTHREADS__
    /* Threaded builds run the whole suite on the worker pool; results
     * must be bitwise identical to the single-thread run (the row
     * partition is deterministic and forces accumulate serially). */
    lbm_set_threads(4);
    printf("-- running with %d threads --\n", lbm_get_threads());
#endif
    /* 1. Vortex-street regime: cylinder, Re ~ 220, 3000 steps. */
    lbm_init(300, 150);
    lbm_set_params(0.10f, 0.f, 0.55f);
    lbm_add_preset(0, 0.25f, 0.5f, 0.25f, 0.f);
    for (int s = 0; s < 3000; s++) lbm_step();
    lbm_compute_fields();
    if (!field_finite()) { printf("FAIL: NaN at Re~220\n"); fails++; }
    float Cl, Cd;
    lbm_get_forces(&Cl, &Cd);
    if (!isfinite(Cl) || !isfinite(Cd) || Cd <= 0.f) {
        printf("FAIL: bad forces Cl=%f Cd=%f\n", Cl, Cd); fails++;
    } else printf("ok   Re~220: Cl=%+.3f Cd=%.3f L=%.0f\n",
                  Cl, Cd, lbm_get_char_length());

    /* 2. Steady regime, Re ~ 40: symmetric wake, 10000 steps. */
    lbm_init(300, 150);
    lbm_set_params(0.08f, 0.f, 0.72f);   /* nu=0.0733, D=38 -> Re~41 */
    lbm_add_preset(0, 0.25f, 0.5f, 0.25f, 0.f);
    for (int s = 0; s < 10000; s++) lbm_step();
    lbm_compute_fields();
    if (!field_finite()) { printf("FAIL: NaN at Re~40\n"); fails++; }
    float asym = 0.f, ref = 0.f;
    for (int d = 1; d < 70; d++)
        for (int x = 100; x < 290; x++) {
            float a = uy_f[idx(x, 75 + d)], b = uy_f[idx(x, 75 - d)];
            asym += fabsf(a + b); ref += fabsf(a) + fabsf(b);
        }
    float rel = ref > 0.f ? asym / ref : 0.f;
    if (rel > 0.15f) { printf("FAIL: wake asymmetry %.3f\n", rel); fails++; }
    else printf("ok   Re~40 wake symmetry: rel asym %.4f\n", rel);

    /* 3. Rapid resize x10 with a NACA 0012 present. */
    lbm_add_preset(2, 0.4f, 0.5f, 0.45f, 0.f);
    static const int sizes[10][2] = {
        {100,50},{600,300},{150,80},{512,256},{200,100},
        {450,225},{120,60},{600,300},{300,150},{300,150}};
    for (int r = 0; r < 10; r++) {
        if (lbm_resize(sizes[r][0], sizes[r][1]) != 0) {
            printf("FAIL: resize %d\n", r); fails++; break;
        }
        for (int s = 0; s < 50; s++) lbm_step();
    }
    lbm_compute_fields();
    if (!field_finite()) { printf("FAIL: NaN after resizes\n"); fails++; }
    else printf("ok   10 rapid resizes, fields finite\n");
#ifdef __EMSCRIPTEN__
    printf("heap after resizes: %zu bytes\n",
           (size_t)emscripten_get_heap_size());
#endif

    /* 4. Painting while running. */
    for (int s = 0; s < 500; s++) {
        lbm_step();
        int x = 150 + (s % 50), y = 40 + (s % 30);
        lbm_set_obstacle(x, y, 1);
        if (s % 3 == 0) lbm_set_obstacle(x, y, 0);
    }
    lbm_compute_fields();
    if (!field_finite()) { printf("FAIL: NaN after painting\n"); fails++; }
    else printf("ok   painting during run\n");

    /* 5. High-Re stress test: tau at the floor, fast inflow, LES on.
     * Without the subgrid model this blows up in checkerboard pressure
     * oscillations within a few hundred steps. */
    lbm_init(400, 200);
    lbm_set_params(0.20f, 0.f, 0.51f);   /* nu=0.0033, D=50 -> Re~3000 */
    lbm_set_les(1, 0.18f);
    lbm_add_preset(0, 0.20f, 0.5f, 0.25f, 0.f);
    for (int s = 0; s < 6000; s++) lbm_step();
    lbm_compute_fields();
    if (!field_finite()) { printf("FAIL: NaN at Re~3000 with LES\n"); fails++; }
    else {
        float Cl5, Cd5;
        lbm_get_forces(&Cl5, &Cd5);
        printf("ok   Re~3000 LES stable: Cl=%+.3f Cd=%.3f\n", Cl5, Cd5);
    }

    /* 6. Shear-layer and jet inlet profiles. */
    for (int prof = 1; prof <= 2; prof++) {
        lbm_clear_obstacles();
        lbm_set_flow_profile(prof);
        lbm_reset_flow();
        for (int s = 0; s < 1500; s++) lbm_step();
        lbm_compute_fields();
        if (!field_finite()) {
            printf("FAIL: NaN with flow profile %d\n", prof); fails++;
        } else printf("ok   flow profile %d stable\n", prof);
    }
    lbm_set_flow_profile(0);

    /* 7. Magnus effect: a spinning cylinder must generate lift whose sign
     * follows the spin direction. Run the same configuration with the
     * spin reversed and require a clear, opposite-signed Cl both times. */
    float cl_spin[2] = { 0.f, 0.f };
    for (int dir = 0; dir < 2; dir++) {
        lbm_init(300, 150);
        lbm_set_params(0.10f, 0.f, 0.60f);
        lbm_add_preset(0, 0.30f, 0.5f, 0.25f, 0.f);
        lbm_set_shape_omega(0, dir == 0 ? 0.10f : -0.10f);
        for (int s = 0; s < 4000; s++) lbm_step();
        lbm_compute_fields();
        if (!field_finite()) { printf("FAIL: NaN with spin %d\n", dir); fails++; }
        float Cd7;
        lbm_get_forces(&cl_spin[dir], &Cd7);
    }
    if (!isfinite(cl_spin[0]) || !isfinite(cl_spin[1])
        || cl_spin[0] * cl_spin[1] >= 0.f
        || fabsf(cl_spin[0]) < 0.05f || fabsf(cl_spin[1]) < 0.05f) {
        printf("FAIL: Magnus lift Cl(+)=%f Cl(-)=%f\n",
               cl_spin[0], cl_spin[1]);
        fails++;
    } else printf("ok   Magnus: Cl(+spin)=%+.3f Cl(-spin)=%+.3f\n",
                  cl_spin[0], cl_spin[1]);

    /* 8. Spin survives a resize (registry replays omega analytically). */
    if (lbm_resize(200, 100) != 0) { printf("FAIL: resize w/ spin\n"); fails++; }
    for (int s = 0; s < 800; s++) lbm_step();
    lbm_compute_fields();
    float cl8, cd8;
    lbm_get_forces(&cl8, &cd8);
    if (!field_finite() || !isfinite(cl8)) {
        printf("FAIL: NaN after resize with spinning cylinder\n"); fails++;
    } else printf("ok   spin after resize: Cl=%+.3f Cd=%.3f\n", cl8, cd8);

    /* 9. Very high Re: tau at the new floor, fast inflow, LES carrying
     * all the stabilisation. Re = u D / nu ~ 0.25*50/3.3e-4 ~ 37000. */
    lbm_init(400, 200);
    lbm_set_params(0.25f, 0.f, 0.501f);
    lbm_set_les(1, 0.18f);
    lbm_add_preset(0, 0.20f, 0.5f, 0.25f, 0.f);
    for (int s = 0; s < 8000; s++) lbm_step();
    lbm_compute_fields();
    if (!field_finite()) { printf("FAIL: NaN at Re~4e4 with LES\n"); fails++; }
    else {
        float cl9, cd9;
        lbm_get_forces(&cl9, &cd9);
        if (!isfinite(cl9) || !isfinite(cd9) || cd9 <= 0.f) {
            printf("FAIL: bad forces at Re~4e4 Cl=%f Cd=%f\n", cl9, cd9);
            fails++;
        } else printf("ok   Re~4e4 LES stable: Cl=%+.3f Cd=%.3f\n", cl9, cd9);
    }

#ifndef USE_MRT
    /* 10. Regularized collision alone (LES off) at low tau: bare BGK
     * blows up within a few hundred steps in this configuration; the
     * Hermite ghost-mode projection must carry it on its own. */
    lbm_init(400, 200);
    lbm_set_params(0.20f, 0.f, 0.505f);   /* nu=0.00167, D=50 -> Re~6000 */
    lbm_set_les(0, 0.f);
    lbm_set_regularized(1);
    lbm_add_preset(0, 0.20f, 0.5f, 0.25f, 0.f);
    for (int s = 0; s < 5000; s++) lbm_step();
    lbm_compute_fields();
    if (!field_finite()) { printf("FAIL: NaN regularized-only Re~6e3\n"); fails++; }
    else printf("ok   regularized BGK (LES off) stable at Re~6e3\n");
    lbm_set_les(1, 0.18f);
#endif

    printf(fails ? "== %d FAILURES ==\n" : "== all tests passed ==\n", fails);
    return fails ? 1 : 0;
}
#endif /* LBM_TEST_MAIN */
