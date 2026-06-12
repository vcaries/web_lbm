# web_lbm — 2D Interactive Wind Tunnel

A self-contained, framework-free web application: a real-time 2D wind
tunnel driven by a Lattice Boltzmann engine written in C and compiled to
WebAssembly, rendered on the GPU with WebGL2. Designed for static hosting
(GitHub Pages) — the deployable artifacts are plain `.html`, `.js`, `.css`
and `.wasm` files with relative paths and no build step on the host.

![architecture] C (LBM, D2Q9) → wasm memory → zero-copy Float32Array views
→ WebGL2 `R32F` texture → colormap LUT in the fragment shader.

## Build

Prerequisite: the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html)
(`emcc` on PATH) — **or** Docker (both build scripts fall back to the
official `emscripten/emsdk` image automatically).

```sh
bash build.sh          # Linux / macOS / Git Bash — BGK collision
bash build.sh --mrt    # MRT collision operator (more stable at high Re)
```

```bat
build.bat              :: Windows cmd, same flags / same fallback
```

Output: two engine variants, each ES-module glue + wasm:

- `dist/lbm_engine.js` — single-threaded, works everywhere;
- `dist/lbm_engine_mt.js` — built with `-pthread` (WebAssembly threads on
  `SharedArrayBuffer`). Threads need the page to be `crossOriginIsolated`
  (COOP/COEP headers); static hosts like GitHub Pages cannot send those,
  so `app/coi-serviceworker.min.js` injects them via a service worker
  (one automatic reload on first visit). The app feature-detects at boot
  and falls back to the single-threaded engine seamlessly — the HUD shows
  the active thread count.
Every emcc flag is documented inline in `build.sh`; the two that matter
most:

- `ALLOW_MEMORY_GROWTH=1` — required because `lbm_resize()` frees and
  re-mallocs every field array at runtime; the wasm heap must be able to
  grow beyond its initial size when the user picks a larger grid.
- `MODULARIZE=1` + `EXPORT_ES6=1` — the runtime is a factory function
  default-exported as an ES module: clean namespacing, multiple engine
  instances per page, `import createLBM from '../dist/lbm_engine.js'`
  with no bundler.

## Run

Wasm cannot be loaded from `file://`; serve the repo root over HTTP:

```sh
python -m http.server 8000
# then open http://localhost:8000/app/
```

To deploy on GitHub Pages, publish the repo (or copy `app/` + `dist/`
preserving their relative layout) — all asset paths are relative.
Note `dist/` is gitignored for development; either build in CI, commit it
on a deployment branch, or remove the ignore entry when publishing.

## Layout

```
src/lbm.c        entire C engine (D2Q9 LBM, BGK / MRT, LES, Zou-He BCs, forces)
dist/            emcc output (gitignored during dev)
app/index.html   standalone entry point
app/main.js      app init, RAF loop, UI wiring, zero-copy views
app/renderer.js  WebGL2 setup, shaders, render()
app/colormaps.js 256x1 RGBA LUTs (Jet, Inferno, Coolwarm, Viridis)
app/obstacles.js analytic preset shapes + NACA 4-digit generator
app/style.css    layout & theme
build.sh / .bat  Emscripten build (flags documented inline)
```

## Physics notes

**Method.** Lattice Boltzmann, D2Q9. Each grid node carries nine particle
populations `f_i` along the discrete velocities `e_i`. One time step is
*collision* (relax toward the local Maxwell equilibrium) followed by
*streaming* (each population hops one link). Macroscopic fields are
moments: `ρ = Σ f_i`, `ρu = Σ e_i f_i`, pressure `p = ρ/3` (lattice
units, `c_s² = 1/3`).

**Collision.** BGK single-relaxation-time by default:
`f_i ← f_i − (f_i − f_i^eq)/τ`, with kinematic viscosity
`ν = (τ − ½)/3`. Compile with `--mrt` for the multiple-relaxation-time
operator (Lallemand & Luo 2000), which relaxes non-hydrodynamic moments
at separate rates and stays stable at higher Reynolds numbers.

**Turbulence model.** A Smagorinsky subgrid closure (LES) is built in and
on by default (`Cs = 0.18`, toggleable in the UI). The non-equilibrium
momentum flux `Π` gives the local strain rate directly, so the effective
relaxation time has the closed form
`τ_eff = ½ (τ + √(τ² + 18 Cs² |Π| / ρ))` — eddy viscosity grows exactly
where gradients steepen, suppressing the checkerboard oscillations BGK
develops at low `τ`. Under MRT the same closure is applied to the
traceless stress moments.

**Flow profiles.** The inlet/initial condition is selectable at runtime
(`lbm_set_flow_profile`): uniform stream, shear layer, or jet.
"Reset flow" (`lbm_reset_flow`) re-initializes the populations to the
current profile while keeping all obstacles.

**Boundaries.**

| Boundary       | Condition                          | Implementation              |
|----------------|------------------------------------|-----------------------------|
| Inlet (left)   | uniform velocity at angle θ        | Zou-He velocity BC          |
| Outlet (right) | non-reflecting outflow             | zero-gradient extrapolation |
| Top/bottom     | selectable no-slip / free-slip     | bounce-back / specular      |
| Obstacles      | no-slip (wall velocity if rotating)| (moving-wall) halfway b-b   |

**Rotating cylinders (Magnus effect).** Each placed cylinder can spin: its
cells carry the solid-body rotation velocity `u_w = u_s/R · (-(y-y_c), x-x_c)`
(`u_s` = surface speed, set per cylinder in the UI), and the bounce-back
picks it up through the moving-wall term (Ladd 1994)

```
f_opp(x, t+1) = f_i^pc(x) - 6 w_i ρ (e_i · u_w)
```

so a spinning cylinder drags the fluid around itself and deflects the wake —
the lift this generates flips sign with the spin direction (verified in the
test suite). The spin is stored as a surface speed, not an angular rate, so
it is resolution-independent and survives grid resizes.

**Forces.** Drag and lift are integrated in C with the momentum-exchange
method (Ladd 1994): every boundary link cut by the obstacle surface
contributes `2 e_i f_i^post-collision` during streaming. Coefficients are
normalized by `½ ρ₀ u₀² L` with `L` the streamwise extent of the obstacle
(chord for airfoils, diameter for the cylinder), and projected onto the
freestream direction (drag) and its perpendicular (lift).

**Acoustics views.** LBM is weakly compressible — it genuinely carries
sound waves at `c_s = 1/√3` (lattice units), so two derived fields make
the acoustics visible. **Schlieren** renders `|∇ρ|` in grayscale like the
experimental knife-edge technique: wave fronts, shear layers, and wake
structure stand out as bright filaments. **Dilatation** renders `∇·u`,
which vanishes for incompressible vortical motion and therefore isolates
the acoustic radiation — with the diverging colormap you can see the
pressure wavefronts a shedding cylinder emits (aeolian tones) and the
startup transient racing through the domain at the speed of sound. Both
are central-difference fields computed alongside vorticity.

**Virtual microphones.** Up to four pressure probes can be placed in the
tunnel (Microphones panel). Each samples `p = ρ/3` once per frame into a
4096-sample ring buffer; the bottom-left panel shows the Hann-windowed
FFT of the pressure fluctuation per probe, with the **Strouhal number**
`St = f L / u₀` of the dominant peak — for a cylinder von Kármán street
expect St ≈ 0.2. Probe positions are stored normalized, so they survive
grid resizes (recordings restart).

**Reynolds number.** `Re = u₀ L / ν` is displayed live; with the default
cylinder (`D ≈ 0.25 Ny`), `τ = 0.55` and `u₀ = 0.1` gives `Re ≈ 220` — a
clear von Kármán vortex street in the vorticity view. Raise `τ` toward
`0.7–0.9` for the steady symmetric-wake regime (`Re ≈ 40`).

**Acoustic noise control.** Being weakly compressible, the lattice also
supports spurious short-wavelength sound waves — visible as "wrinkles" in
the pressure/schlieren views and, worse, able to reflect off the open
boundaries and pollute the flow. Two mechanisms remove them: (1) inside
the regularized reconstruction the *trace* of `Π^neq` (the dilatational
stress) is scaled by 0.4, an enhanced bulk viscosity that damps acoustic
modes like k² — grid-scale waves die fastest, vortices (traceless part)
are untouched; (2) **sponge layers** over the ~25 columns nearest the
inlet and outlet relax the populations toward the freestream equilibrium
with quadratically increasing strength, absorbing outgoing waves and the
wake before the boundaries can reflect them. These are what make the
high-Re force readouts meaningful (cylinder Cd ≈ 1.9–2.6 from Re 3·10³
to 4·10⁴ in the test suite, instead of the unphysical values reflected
acoustics used to drive).

**Stability envelope.** Four mechanisms stack, all verified in the test
suite (stable to Re ≈ 4·10⁴):

1. **Hermite-regularized collision** (Latt & Chopard 2006, on by
   default): before relaxing, the non-equilibrium part of each
   population is replaced by its projection onto the stress basis,
   `f_i^neq → (w_i/2c_s⁴) Q_i : Π^neq` — the non-hydrodynamic "ghost"
   modes that destabilise BGK at low `τ` are discarded every step.
   This alone holds a Re ≈ 6000 cylinder with the LES *off*.
2. **Smagorinsky LES** (on by default): eddy viscosity grows exactly
   where gradients steepen.
3. **Soft start**: the inlet speed ramps over the first 400 steps after
   any (re)initialisation, removing the impulsive-start transient that
   is otherwise the most fragile moment at high Re.
4. **Hard guards**: the UI clamps `u₀ ≤ 0.25` (low-Mach validity) and
   `τ ≥ 0.502` (engine floor 0.501, `ν ≈ 3.3·10⁻⁴`); the kernel caps
   `|u|` at 0.3 and floors `ρ` at 0.05 so even a pathological
   configuration degrades gracefully instead of producing NaNs.

For the best-defined von Kármán streets pair a low `τ` with a higher
grid resolution.

## Performance

Collision and streaming are fused into a single push-scheme pass: each
cell is collided once and its post-collision populations are scattered
straight to the neighbours (`f_i(x+e_i, t+1) = f_i^pc(x,t)`), halving the
memory traffic of the classic two-kernel formulation. The interior sweep
is completely branchless — solid cells are collided like fluid (harmless:
their populations sit at an equilibrium fixed point) and fluid/solid
links are repaired afterwards in a cheap fixup pass that also accumulates
the drag/lift force and applies the moving-wall term. With
`restrict`-qualified pointers LLVM auto-vectorizes the whole loop to
WebAssembly SIMD; measured against the previous kernel (Node, same
flags):

| Grid      | before      | after       | speed-up |
|-----------|-------------|-------------|----------|
| 300×150   | 1.96 ms/step| 0.73 ms/step| 2.7×     |
| 600×300   | 9.40 ms/step| 3.14 ms/step| 3.0×     |
| 1000×500  | 26.0 ms/step| 8.6 ms/step | 3.0×     |

At 4 steps/frame this keeps even a 1000-cell-long tunnel within a 60 fps
frame budget (~34 ms → use 1–2 steps/frame), and the default 300×150 grid
costs under 3 ms of simulation per frame.

On top of that, the multithreaded engine row-partitions the interior
sweep and the solid-link fixup across a pthread pool (race-free by
construction: every destination slot has exactly one producer cell;
forces reduce in thread order, so results are deterministic for a given
thread count). Measured on a 4-core laptop under Node: 1.3–1.8× over the
vectorized single-thread kernel — LBM is memory-bandwidth-bound, so
machines with more memory channels scale further. Finally, only the
derivative field being displayed (vorticity / schlieren / dilatation) is
computed each frame; the others are skipped.

## Zero-copy data path

The C engine mallocs all field arrays inside `WebAssembly.Memory`.
JavaScript wraps the raw pointers (`lbm_get_*_ptr()`) in
`Float32Array`/`Uint8Array` views and passes those views directly to
`gl.texSubImage2D` — the scalar field travels C → GPU without a single
intermediate copy or serialization. Views are rebuilt only when wasm
memory growth swaps the backing `ArrayBuffer` (detected by comparing
buffer identity each frame) or after `lbm_resize()` reallocates.

## Resizing & persistent obstacles

`lbm_resize(nx, ny)` frees and reallocates every array, then restores the
obstacle layout: hand-painted cells are resampled nearest-neighbour from
the old mask, while preset shapes (cylinder, plate, NACA airfoils) are
re-rasterized **analytically** from a shape registry, so they stay crisp
at any resolution. The NACA 4-digit profiles are generated from the
standard thickness/camber equations (see `app/obstacles.js` and
`raster_naca` in `src/lbm.c`) — no bitmaps anywhere.

## Testing

A smoke-test harness ships inside `src/lbm.c` behind `-DLBM_TEST_MAIN`
(steady-wake symmetry at Re≈40, vortex-street regime, 10 rapid resizes,
painting while running, Re≈3000 LES stability, shear/jet profiles, Magnus
lift sign reversal for spinning cylinders, spin-after-resize, NaN/Inf
scans). Run it under Node:

```sh
emcc src/lbm.c -O2 -DLBM_TEST_MAIN -s ENVIRONMENT=node -s ALLOW_MEMORY_GROWTH=1 -o tmp_test/test.js
node tmp_test/test.js
```

(`ALLOW_MEMORY_GROWTH=1` is required — the resize test grows the grid past
the default heap.)

Without a local `emcc`, compile through the same Docker image the build
scripts use (Git Bash / WSL syntax):

```sh
mkdir -p tmp_test
MSYS_NO_PATHCONV=1 docker run --rm -v "$(pwd -W 2>/dev/null || pwd):/src" -w /src \
  emscripten/emsdk emcc src/lbm.c -O2 -DLBM_TEST_MAIN -s ENVIRONMENT=node \
  -s ALLOW_MEMORY_GROWTH=1 -o tmp_test/test.js
node tmp_test/test.js
```

Add `-DUSE_MRT` to either command to test the MRT collision operator.
