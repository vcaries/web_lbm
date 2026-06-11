/*
 * obstacles.js — analytic preset shapes rasterized to the current grid.
 *
 * Every generator returns an array of [x, y] grid cells; nothing is a
 * stored bitmap, so shapes can be re-rasterized crisply at any grid
 * resolution (used when the user resizes the grid mid-run).
 *
 * Coordinates: cx in [0,1] of Nx, cy in [0,1] of Ny, scale as a fraction
 * of Ny (diameter for the cylinder, length/chord for plate and airfoils).
 * angleDeg rotates the shape counter-clockwise (grid y points up).
 */

/** Circle: solid where (x-xc)^2 + (y-yc)^2 <= r^2. */
export function rasterizeCylinder(Nx, Ny, cx, cy, scale) {
  const xc = cx * Nx, yc = cy * Ny, r = 0.5 * scale * Ny;
  const cells = [];
  for (let y = Math.max(0, Math.floor(yc - r) - 1);
       y <= Math.min(Ny - 1, Math.ceil(yc + r) + 1); y++) {
    for (let x = Math.max(0, Math.floor(xc - r) - 1);
         x <= Math.min(Nx - 1, Math.ceil(xc + r) + 1); x++) {
      const dx = x - xc, dy = y - yc;
      if (dx * dx + dy * dy <= r * r) cells.push([x, y]);
    }
  }
  return cells;
}

/** Thin rectangle (length scale*Ny, thickness 4% of length, min 2 cells),
 *  rotated by angleDeg: a cell is solid if its center, inverse-rotated
 *  into the plate frame, lies inside the half-length/half-thickness box. */
export function rasterizePlate(Nx, Ny, cx, cy, scale, angleDeg) {
  const xc = cx * Nx, yc = cy * Ny;
  const len = scale * Ny, th = Math.max(2, 0.04 * len);
  const a = (angleDeg * Math.PI) / 180;
  const ca = Math.cos(a), sa = Math.sin(a);
  const rad = 0.5 * Math.hypot(len, th) + 1;
  const cells = [];
  for (let y = Math.max(0, Math.floor(yc - rad));
       y <= Math.min(Ny - 1, Math.ceil(yc + rad)); y++) {
    for (let x = Math.max(0, Math.floor(xc - rad));
         x <= Math.min(Nx - 1, Math.ceil(xc + rad)); x++) {
      const dx = x - xc, dy = y - yc;
      const lx = dx * ca + dy * sa;    // along the plate
      const ly = -dx * sa + dy * ca;   // across the plate
      if (Math.abs(lx) <= len / 2 && Math.abs(ly) <= th / 2)
        cells.push([x, y]);
    }
  }
  return cells;
}

/**
 * NACA 4-digit profile, generated from the standard equations.
 *
 * For a code MPXX (e.g. 4412): m = M/100 is the maximum camber,
 * p = P/10 its chordwise position, t = XX/100 the thickness, all as
 * fractions of the chord c. With s = x/c in [0,1]:
 *
 *   half-thickness  yt(s) = 5t (0.2969 sqrt(s) - 0.1260 s - 0.3516 s^2
 *                             + 0.2843 s^3 - 0.1036 s^4)   (closed TE)
 *   camber line     yc(s) = m/p^2 (2ps - s^2)                    s <  p
 *                   yc(s) = m/(1-p)^2 ((1-2p) + 2ps - s^2)       s >= p
 *
 * A grid cell is solid when, in the (rotated) chord frame, its offset
 * from the camber line is within the local half-thickness:
 * |y'/c - yc(s)| <= yt(s). Thickness is applied vertically about the
 * camber line (standard small-camber approximation).
 */
export function rasterizeNACA(code, Nx, Ny, cx, cy, scale, angleDeg) {
  const m = parseInt(code[0], 10) / 100;
  const p = parseInt(code[1], 10) / 10;
  const t = parseInt(code.slice(2), 10) / 100;
  const xc = cx * Nx, yc = cy * Ny, chord = scale * Ny;
  const a = (angleDeg * Math.PI) / 180;
  const ca = Math.cos(a), sa = Math.sin(a);
  const rad = 0.75 * chord + 2;
  const cells = [];
  for (let y = Math.max(0, Math.floor(yc - rad));
       y <= Math.min(Ny - 1, Math.ceil(yc + rad)); y++) {
    for (let x = Math.max(0, Math.floor(xc - rad));
         x <= Math.min(Nx - 1, Math.ceil(xc + rad)); x++) {
      const dx = x - xc, dy = y - yc;
      const lx = dx * ca + dy * sa;    // chordwise, centered on mid-chord
      const ly = -dx * sa + dy * ca;
      const s = (lx + chord / 2) / chord;
      if (s < 0 || s > 1) continue;
      const yt = 5 * t * (0.2969 * Math.sqrt(s) - 0.1260 * s
                 - 0.3516 * s * s + 0.2843 * s ** 3 - 0.1036 * s ** 4);
      let ycam = 0;
      if (m > 0) {
        ycam = s < p
          ? (m / (p * p)) * (2 * p * s - s * s)
          : (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * s - s * s);
      }
      // 0.6-cell half-thickness floor: keeps the thin trailing edge
      // watertight at any grid resolution (cell-centre sampling would
      // otherwise leave gaps where yt < 1 cell). The floor is in grid
      // units, so the profile converges to the true shape as the
      // resolution rises. Mirrors raster_naca in src/lbm.c.
      if (Math.abs(ly - ycam * chord) <= Math.max(yt * chord, 0.6))
        cells.push([x, y]);
    }
  }
  return cells;
}

/** Preset descriptors: default placement (overridable by click-to-place)
 *  and size. Defaults sit well upstream so the wake has room to develop. */
export const PRESETS = {
  cylinder:  { label: 'Cylinder',   cx: 0.18, cy: 0.50, scale: 0.25 },
  plate:     { label: 'Flat plate', cx: 0.20, cy: 0.50, scale: 0.40 },
  naca0012:  { label: 'NACA 0012',  cx: 0.22, cy: 0.50, scale: 0.45 },
  naca4412:  { label: 'NACA 4412',  cx: 0.22, cy: 0.50, scale: 0.45 },
};

/** Rasterize any preset by name at the given grid size; cx/cy (normalized
 *  [0,1]) override the default placement when given. */
export function rasterizePreset(name, Nx, Ny, angleDeg, cx, cy) {
  const p = PRESETS[name];
  const x = cx ?? p.cx, y = cy ?? p.cy;
  switch (name) {
    case 'cylinder': return rasterizeCylinder(Nx, Ny, x, y, p.scale);
    case 'plate':    return rasterizePlate(Nx, Ny, x, y, p.scale, angleDeg);
    case 'naca0012': return rasterizeNACA('0012', Nx, Ny, x, y, p.scale, angleDeg);
    case 'naca4412': return rasterizeNACA('4412', Nx, Ny, x, y, p.scale, angleDeg);
    default: return [];
  }
}

/**
 * Closed outline polygon for the ghost preview shown while placing a
 * shape. Returns points in grid-cell units, centered on the shape origin
 * (translate by the placement point and flip y to draw on a canvas).
 * Same analytic definitions as the rasterizers, just sampled as a contour.
 */
export function outlinePreset(name, Ny, angleDeg) {
  const p = PRESETS[name];
  const a = (angleDeg * Math.PI) / 180;
  const ca = Math.cos(a), sa = Math.sin(a);
  const rot = ([x, y]) => [x * ca - y * sa, x * sa + y * ca];
  const pts = [];
  if (name === 'cylinder') {
    const r = 0.5 * p.scale * Ny;
    for (let i = 0; i < 48; i++) {
      const t = (i / 48) * 2 * Math.PI;
      pts.push([r * Math.cos(t), r * Math.sin(t)]);
    }
    return pts;
  }
  if (name === 'plate') {
    const len = p.scale * Ny, th = Math.max(2, 0.04 * len);
    return [
      rot([-len / 2, -th / 2]), rot([len / 2, -th / 2]),
      rot([len / 2, th / 2]), rot([-len / 2, th / 2]),
    ];
  }
  // NACA: sample upper surface nose->tail, then lower surface tail->nose.
  const code = name === 'naca0012' ? '0012' : '4412';
  const m = parseInt(code[0], 10) / 100;
  const pp = parseInt(code[1], 10) / 10;
  const t = parseInt(code.slice(2), 10) / 100;
  const chord = p.scale * Ny;
  const surf = (s) => {
    const yt = 5 * t * (0.2969 * Math.sqrt(s) - 0.1260 * s
               - 0.3516 * s * s + 0.2843 * s ** 3 - 0.1036 * s ** 4);
    let yc = 0;
    if (m > 0) {
      yc = s < pp
        ? (m / (pp * pp)) * (2 * pp * s - s * s)
        : (m / ((1 - pp) * (1 - pp))) * ((1 - 2 * pp) + 2 * pp * s - s * s);
    }
    return [yc + yt, yc - yt];
  };
  const N = 36;
  for (let i = 0; i <= N; i++) {
    const s = i / N;
    pts.push(rot([(s - 0.5) * chord, surf(s)[0] * chord]));
  }
  for (let i = N; i >= 0; i--) {
    const s = i / N;
    pts.push(rot([(s - 0.5) * chord, surf(s)[1] * chord]));
  }
  return pts;
}
