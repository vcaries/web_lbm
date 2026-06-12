/*
 * main.js — application entry point.
 *
 * Owns the wasm module instance, the requestAnimationFrame loop, the UI
 * wiring, and the zero-copy TypedArray views over the engine's memory.
 *
 * Zero-copy contract: the C engine mallocs every field array inside
 * WebAssembly.Memory. We wrap the raw pointers in Float32Array/Uint8Array
 * views and pass those views directly to the WebGL2 renderer. Views are
 * rebuilt only when they detach (wasm memory growth swaps the underlying
 * ArrayBuffer) or when lbm_resize() reallocates the arrays.
 */

import { Renderer } from './renderer.js';
import { COLORMAPS } from './colormaps.js';
import { rasterizePreset, outlinePreset, PRESETS } from './obstacles.js';

const $ = (id) => document.getElementById(id);

/* Engine selection: the -pthread build needs crossOriginIsolated
 * (SharedArrayBuffer). coi-serviceworker provides that on GitHub Pages
 * after a one-time reload; anywhere it isn't available we fall back to
 * the single-threaded engine, which is feature-identical. */
let engineMT = false;
let threadCount = 1;
async function loadEngine() {
  if (window.crossOriginIsolated) {
    try {
      const mod = await import('../dist/lbm_engine_mt.js');
      engineMT = true;
      return mod.default;
    } catch (err) {
      console.warn('Multithreaded engine unavailable, using single-threaded:', err);
    }
  }
  return (await import('../dist/lbm_engine.js')).default;
}

let Module;            // Emscripten module instance
let renderer;
let Nx = 300, Ny = 150;
let views = null;      // TypedArray views over wasm memory
let forcesPtr = 0;     // 8-byte scratch for lbm_get_forces out-params

// UI state
let fieldMode = 'velocity';        // velocity | pressure | vorticity
let scaleMode = 'dynamic';         // dynamic | fixed
let stepsPerFrame = 4;
let paintMode = 'draw';            // draw | erase
let stepCount = 0;
let msFrame = 0;                   // EMA of frame-to-frame interval (ms)
let lastFrameT = 0;                // rAF timestamp of the previous frame
let obstacleDirty = true;
let frame = 0;
let placing = null;                // preset name while click-to-place is armed
let ghostPos = null;               // normalized [u, v] cursor pos for the ghost
let overlayMode = 'none';          // none | streamlines | arrows | particles
let octx = null;                   // 2D context of the overlay canvas
let domainAspect = 2.0;            // tunnel length : height (domlen slider)

// Tracer particles: flat [x0,y0, x1,y1, ...] in fractional grid cells,
// advected through the live velocity field each frame.
const PARTICLE_N = 3500;
let particles = null;
let particleSize = 2;              // CSS px (scaled by devicePixelRatio)
let particleColor = '#ffffff';

// C-side preset registry mirror (cylinders & plates, in lbm_add_preset
// call order) so per-shape rotation can be addressed by index.
const cShapes = [];

// Virtual microphones: pressure probes (p = rho/3) sampled once per frame
// into a ring buffer, FFT'd into the spectrum panel. Positions are
// normalized [0,1]^2 so they survive grid resizes.
const MIC_MAX = 4;
const MIC_N = 4096;                // ring length (frames) — sets df
const MIC_COLORS = ['#ffd166', '#06d6a0', '#ef476f', '#3b9eff'];
const mics = [];                   // {u, v, buf, head, count}
let placingMic = false;

// Presets drawn through the analytic JS rasterizers (the NACA airfoils,
// per spec). Replayed after every grid resize. Cylinder/plate go through
// the C engine's own registry (lbm_add_preset), which replays them itself.
const jsShapes = [];

const DEFAULT_CMAP = {
  velocity: 'INFERNO', pressure: 'COOLWARM', vorticity: 'COOLWARM',
  schlieren: 'GRAY', dilatation: 'COOLWARM',
};
// Derivative-field bitmask for lbm_set_active_field (1 vorticity,
// 2 schlieren, 4 dilatation) — skips unused finite-difference passes.
const FIELD_MASKS = {
  velocity: 0, pressure: 0, vorticity: 1, schlieren: 2, dilatation: 4,
};
const FIXED_DEFAULTS = {
  velocity:   { min: 0,      max: 0.2   },
  pressure:   { min: 0.95,   max: 1.05  },
  vorticity:  { min: -0.05,  max: 0.05  },
  schlieren:  { min: 0,      max: 0.005 },
  dilatation: { min: -0.002, max: 0.002 },
};

/* ------------------------------------------------------------------ */
/* Wasm memory views                                                    */
/* ------------------------------------------------------------------ */

function heapBuffer() {
  return Module.wasmMemory ? Module.wasmMemory.buffer : Module.HEAPF32.buffer;
}

/** (Re)build the TypedArray views over the engine's field arrays. */
function refreshViews() {
  const buf = heapBuffer();
  const n = Nx * Ny;
  views = {
    speed:    new Float32Array(buf, Module._lbm_get_speed_ptr(), n),
    rho:      new Float32Array(buf, Module._lbm_get_rho_ptr(), n),
    vort:     new Float32Array(buf, Module._lbm_get_vorticity_ptr(), n),
    schl:     new Float32Array(buf, Module._lbm_get_schlieren_ptr(), n),
    dil:      new Float32Array(buf, Module._lbm_get_dilatation_ptr(), n),
    ux:       new Float32Array(buf, Module._lbm_get_ux_ptr(), n),
    uy:       new Float32Array(buf, Module._lbm_get_uy_ptr(), n),
    obstacle: new Uint8Array(buf, Module._lbm_get_obstacle_ptr(), n),
  };
}

/** Views detach when wasm memory grows; rebuild on demand. */
function ensureViews() {
  if (!views || views.speed.buffer !== heapBuffer()
      || views.speed.length !== Nx * Ny) {
    refreshViews();
  }
}

function currentFieldView() {
  ensureViews();
  if (fieldMode === 'pressure') return views.rho;
  if (fieldMode === 'vorticity') return views.vort;
  if (fieldMode === 'schlieren') return views.schl;
  if (fieldMode === 'dilatation') return views.dil;
  return views.speed;
}

/* ------------------------------------------------------------------ */
/* Simulation parameter plumbing                                        */
/* ------------------------------------------------------------------ */

function pushParams() {
  const u0 = parseFloat($('vel').value);
  const angle = parseFloat($('angle').value);
  const tau = parseFloat($('tau').value);
  Module._lbm_set_params(u0, angle, tau);
  $('vel-val').textContent = u0.toFixed(3);
  $('angle-val').textContent = `${angle.toFixed(1)}°`;
  $('tau-val').textContent = tau.toFixed(3);
}

function derivedRe() {
  const u0 = parseFloat($('vel').value);
  const tau = parseFloat($('tau').value);
  const nu = (tau - 0.5) / 3;
  const L = Module._lbm_get_char_length();
  return (u0 * L) / nu;
}

/* ------------------------------------------------------------------ */
/* Obstacles: presets, painting, resize replay                          */
/* ------------------------------------------------------------------ */

function paintCells(cells, value) {
  for (const [x, y] of cells) Module._lbm_set_obstacle(x, y, value);
  obstacleDirty = true;
}

function addPresetAt(name, cx, cy) {
  const angle = parseFloat($('angle').value);
  const p = PRESETS[name];
  cx = cx ?? p.cx;
  cy = cy ?? p.cy;
  if (name === 'cylinder' || name === 'plate') {
    // C-side path: the engine rasterizes analytically and registers the
    // shape in its own registry, so lbm_resize() replays it natively.
    const id = name === 'cylinder' ? 0 : 1;
    Module._lbm_add_preset(id, cx, cy, p.scale, name === 'plate' ? angle : 0);
    const omega = name === 'cylinder' ? parseFloat($('spin').value) : 0;
    cShapes.push({ name, omega });
    if (omega !== 0) Module._lbm_set_shape_omega(cShapes.length - 1, omega);
    rebuildCylinderList();
  } else {
    // JS-side analytic path (NACA generator, per spec): rasterize at the
    // current grid and remember the shape for replay after resizes.
    paintCells(rasterizePreset(name, Nx, Ny, angle, cx, cy), 2);
    jsShapes.push({ name, angle, cx, cy });
  }
  obstacleDirty = true;
}

/** Arm/disarm click-to-place mode for a preset. */
function setPlacing(name) {
  placing = name;
  ghostPos = null;
  document.querySelectorAll('[data-preset]').forEach((b) =>
    b.classList.toggle('placing', b.dataset.preset === placing));
}

function clearObstacles() {
  Module._lbm_clear_obstacles();
  jsShapes.length = 0;
  cShapes.length = 0;
  rebuildCylinderList();
  obstacleDirty = true;
}

/** One spin slider per placed cylinder (index = C registry order). */
function rebuildCylinderList() {
  const list = $('cyl-list');
  list.textContent = '';
  cShapes.forEach((s, i) => {
    if (s.name !== 'cylinder') return;
    const row = document.createElement('div');
    row.className = 'cyl-row';
    const lab = document.createElement('span');
    lab.textContent = `Cyl ${i + 1}`;
    const slider = document.createElement('input');
    slider.type = 'range';
    slider.min = '-0.15'; slider.max = '0.15'; slider.step = '0.005';
    slider.value = String(s.omega);
    const val = document.createElement('span');
    val.className = 'val';
    val.textContent = s.omega.toFixed(3);
    slider.addEventListener('input', () => {
      s.omega = parseFloat(slider.value);
      val.textContent = s.omega.toFixed(3);
      Module._lbm_set_shape_omega(i, s.omega);
    });
    row.append(lab, slider, val);
    list.append(row);
  });
}

/* Letterbox both canvases to the domain aspect ratio inside the viewport
 * so grid cells always display square, whatever the tunnel length. */
function layoutCanvases() {
  const vw = window.innerWidth, vh = window.innerHeight;
  let w = vw, h = vw / domainAspect;
  if (h > vh) { h = vh; w = vh * domainAspect; }
  const left = (vw - w) / 2, top = (vh - h) / 2;
  for (const c of [$('glcanvas'), $('overlay')]) {
    c.style.left = `${left}px`;
    c.style.top = `${top}px`;
    c.style.width = `${w}px`;
    c.style.height = `${h}px`;
  }
}

function resizeGrid(newNx) {
  const newNy = Math.min(512, Math.max(48, Math.round(newNx / domainAspect)));
  if (Module._lbm_resize(newNx, newNy) !== 0) {
    console.error(`lbm_resize(${newNx}, ${newNy}) failed`);
    return;
  }
  Nx = newNx; Ny = newNy;
  initParticles();
  for (const m of mics) { m.head = 0; m.count = 0; }   // restart recordings
  // Hand-painted cells were resampled and the C registry replayed inside
  // lbm_resize(); now re-rasterize the JS-side (NACA) shapes crisply.
  for (const s of jsShapes) {
    paintCells(rasterizePreset(s.name, Nx, Ny, s.angle, s.cx, s.cy), 2);
  }
  refreshViews();
  renderer.setGrid(Nx, Ny);
  obstacleDirty = true;
  stepCount = 0;
  $('res-val').textContent = `${Nx} × ${Ny}`;
}

/* Painting: pointer drag rasterizes a brush disc along the stroke.
 * When a preset is armed (click-to-place), the same pointer events drive
 * the ghost preview instead, and a click drops the shape. */
function setupPainting(canvas) {
  let painting = false;
  let last = null;

  const toNorm = (e) => {
    const r = canvas.getBoundingClientRect();
    return [(e.clientX - r.left) / r.width,
            1 - (e.clientY - r.top) / r.height];   // grid y=0 is the bottom
  };
  const toGrid = (e) => {
    const [u, v] = toNorm(e);
    return [Math.floor(u * Nx), Math.floor(v * Ny)];
  };

  const dab = (gx, gy) => {
    const rad = Math.max(1, Math.round(Nx / 150));
    const value = paintMode === 'draw' ? 1 : 0;
    for (let dy = -rad; dy <= rad; dy++) {
      for (let dx = -rad; dx <= rad; dx++) {
        if (dx * dx + dy * dy <= rad * rad) {
          Module._lbm_set_obstacle(gx + dx, gy + dy, value);
        }
      }
    }
    obstacleDirty = true;
  };

  const stroke = (e) => {
    const [gx, gy] = toGrid(e);
    if (last) {
      // Interpolate so fast drags leave a continuous wall.
      const steps = Math.max(Math.abs(gx - last[0]), Math.abs(gy - last[1]));
      for (let s = 1; s <= steps; s++) {
        dab(Math.round(last[0] + ((gx - last[0]) * s) / steps),
            Math.round(last[1] + ((gy - last[1]) * s) / steps));
      }
    }
    dab(gx, gy);
    last = [gx, gy];
  };

  canvas.addEventListener('pointerdown', (e) => {
    if (placingMic) {
      const [u, v] = toNorm(e);
      if (mics.length < MIC_MAX) {
        mics.push({ u, v, buf: new Float32Array(MIC_N), head: 0, count: 0 });
      }
      placingMic = false;
      $('mic-add').classList.remove('placing');
      drawSpectrum();
      return;
    }
    if (placing) {
      const [u, v] = toNorm(e);
      addPresetAt(placing, u, v);
      setPlacing(null);
      return;
    }
    painting = true;
    last = null;
    canvas.setPointerCapture(e.pointerId);
    stroke(e);
  });
  canvas.addEventListener('pointermove', (e) => {
    if (placing) { ghostPos = toNorm(e); return; }
    if (painting) stroke(e);
  });
  canvas.addEventListener('pointerleave', () => { if (placing) ghostPos = null; });
  const stop = () => { painting = false; last = null; };
  canvas.addEventListener('pointerup', stop);
  canvas.addEventListener('pointercancel', stop);
  window.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') {
      setPlacing(null);
      placingMic = false;
      $('mic-add').classList.remove('placing');
    }
  });
}

/* ------------------------------------------------------------------ */
/* 2D overlay: ghost preview, streamlines, velocity arrows              */
/* ------------------------------------------------------------------ */

/** Bilinear velocity sample at fractional grid coordinates. */
function sampleVel(gx, gy) {
  const x0 = Math.min(Nx - 2, Math.max(0, Math.floor(gx)));
  const y0 = Math.min(Ny - 2, Math.max(0, Math.floor(gy)));
  const fx = Math.min(1, Math.max(0, gx - x0));
  const fy = Math.min(1, Math.max(0, gy - y0));
  const k00 = y0 * Nx + x0, k10 = k00 + 1, k01 = k00 + Nx, k11 = k01 + 1;
  const ux = (views.ux[k00] * (1 - fx) + views.ux[k10] * fx) * (1 - fy)
           + (views.ux[k01] * (1 - fx) + views.ux[k11] * fx) * fy;
  const uy = (views.uy[k00] * (1 - fx) + views.uy[k10] * fx) * (1 - fy)
           + (views.uy[k01] * (1 - fx) + views.uy[k11] * fx) * fy;
  return [ux, uy];
}

/** Streamlines: midpoint (RK2) integration of the live velocity field,
 *  seeded along the inlet. Stops at obstacles, walls, or stagnation. */
function drawStreamlines(ctx, w, h) {
  const nSeeds = 22;
  const hStep = 0.8;                      // integration step, in cells
  const maxSteps = Math.round(2.5 * Nx);
  ctx.strokeStyle = 'rgba(255, 255, 255, 0.45)';
  ctx.lineWidth = Math.max(1, w / 1500);
  for (let s = 0; s < nSeeds; s++) {
    let gx = 1.5;
    let gy = ((s + 0.5) / nSeeds) * Ny;
    ctx.beginPath();
    ctx.moveTo((gx / Nx) * w, h - (gy / Ny) * h);
    for (let i = 0; i < maxSteps; i++) {
      const [u1, v1] = sampleVel(gx, gy);
      const sp1 = Math.hypot(u1, v1);
      if (sp1 < 1e-6) break;
      const mx = gx + (0.5 * hStep * u1) / sp1;
      const my = gy + (0.5 * hStep * v1) / sp1;
      const [u2, v2] = sampleVel(mx, my);
      const sp2 = Math.hypot(u2, v2);
      if (sp2 < 1e-6) break;
      gx += (hStep * u2) / sp2;
      gy += (hStep * v2) / sp2;
      if (gx < 0 || gx >= Nx - 1 || gy < 0 || gy >= Ny - 1) break;
      if (views.obstacle[Math.round(gy) * Nx + Math.round(gx)]) break;
      ctx.lineTo((gx / Nx) * w, h - (gy / Ny) * h);
    }
    ctx.stroke();
  }
}

/** Velocity glyphs on a coarse grid: direction + length ~ |u| / u0. */
function drawArrows(ctx, w, h) {
  const spacing = Math.max(6, Math.round(Nx / 40));   // in cells
  const u0 = parseFloat($('vel').value);
  const cellW = w / Nx, cellH = h / Ny;
  const maxLen = spacing * 0.85 * cellW;
  ctx.strokeStyle = 'rgba(255, 255, 255, 0.55)';
  ctx.fillStyle = 'rgba(255, 255, 255, 0.55)';
  ctx.lineWidth = Math.max(1, w / 1800);
  for (let gy = spacing >> 1; gy < Ny; gy += spacing) {
    for (let gx = spacing >> 1; gx < Nx; gx += spacing) {
      const k = gy * Nx + gx;
      if (views.obstacle[k]) continue;
      const ux = views.ux[k], uy = views.uy[k];
      const sp = Math.hypot(ux, uy);
      if (sp < 1e-5) continue;
      const len = Math.min(maxLen, (sp / Math.max(u0, 1e-6)) * maxLen * 0.7);
      const px = (gx + 0.5) * cellW, py = h - (gy + 0.5) * cellH;
      const dx = (ux / sp) * len, dy = -(uy / sp) * len;   // canvas y down
      ctx.beginPath();
      ctx.moveTo(px - dx / 2, py - dy / 2);
      ctx.lineTo(px + dx / 2, py + dy / 2);
      ctx.stroke();
      // arrowhead
      const hx = px + dx / 2, hy = py + dy / 2;
      const ah = Math.max(2, len * 0.3);
      const a = Math.atan2(dy, dx);
      ctx.beginPath();
      ctx.moveTo(hx, hy);
      ctx.lineTo(hx - ah * Math.cos(a - 0.45), hy - ah * Math.sin(a - 0.45));
      ctx.lineTo(hx - ah * Math.cos(a + 0.45), hy - ah * Math.sin(a + 0.45));
      ctx.closePath();
      ctx.fill();
    }
  }
}

/* Tracer particles: massless markers advected by the live velocity field
 * (forward Euler with bilinear sampling, inlined to avoid per-particle
 * allocations). Particles that exit the domain, hit an obstacle, or get
 * trapped are re-seeded just behind the inlet. */

function seedParticle(i, anywhere) {
  particles[2 * i]     = anywhere ? Math.random() * (Nx - 2)
                                  : 0.5 + Math.random() * 2;
  particles[2 * i + 1] = 0.5 + Math.random() * (Ny - 1.5);
}

function initParticles() {
  particles = new Float32Array(PARTICLE_N * 2);
  for (let i = 0; i < PARTICLE_N; i++) seedParticle(i, true);
}

function advectParticles() {
  const ux = views.ux, uy = views.uy, obs = views.obstacle;
  const dt = stepsPerFrame;          // the flow advances this many steps/frame
  for (let i = 0; i < PARTICLE_N; i++) {
    let gx = particles[2 * i], gy = particles[2 * i + 1];
    let x0 = gx | 0, y0 = gy | 0;
    if (x0 > Nx - 2) x0 = Nx - 2;
    if (y0 > Ny - 2) y0 = Ny - 2;
    const fx = gx - x0, fy = gy - y0;
    const k = y0 * Nx + x0;
    const vx = (ux[k] * (1 - fx) + ux[k + 1] * fx) * (1 - fy)
             + (ux[k + Nx] * (1 - fx) + ux[k + Nx + 1] * fx) * fy;
    const vy = (uy[k] * (1 - fx) + uy[k + 1] * fx) * (1 - fy)
             + (uy[k + Nx] * (1 - fx) + uy[k + Nx + 1] * fx) * fy;
    gx += vx * dt;
    gy += vy * dt;
    if (gx < 0.5 || gx >= Nx - 1 || gy < 0.5 || gy >= Ny - 1
        || obs[(gy | 0) * Nx + (gx | 0)]) {
      seedParticle(i, false);
    } else {
      particles[2 * i] = gx;
      particles[2 * i + 1] = gy;
    }
  }
}

function drawParticles(ctx, w, h) {
  const sx = w / Nx, sy = h / Ny;
  const r = Math.max(1, Math.round(particleSize * (window.devicePixelRatio || 1)));
  ctx.fillStyle = particleColor;
  ctx.globalAlpha = 0.8;
  ctx.beginPath();
  for (let i = 0; i < PARTICLE_N; i++) {
    ctx.rect(particles[2 * i] * sx, h - particles[2 * i + 1] * sy, r, r);
  }
  ctx.fill();
  ctx.globalAlpha = 1;
}

/* ------------------------------------------------------------------ */
/* Virtual microphones: recording, FFT, spectrum panel                  */
/* ------------------------------------------------------------------ */

function recordMics() {
  for (const m of mics) {
    const gx = Math.min(Nx - 2, Math.max(1, Math.round(m.u * Nx)));
    const gy = Math.min(Ny - 2, Math.max(1, Math.round(m.v * Ny)));
    m.buf[m.head] = views.rho[gy * Nx + gx] / 3;   // p = rho / 3
    m.head = (m.head + 1) % MIC_N;
    m.count++;
  }
}

/** In-place iterative radix-2 complex FFT (n = power of two). */
function fft(re, im) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j |= bit;
    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (-2 * Math.PI) / len;
    const wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let j = 0; j < len / 2; j++) {
        const a = i + j, b = i + j + len / 2;
        const vr = re[b] * cr - im[b] * ci;
        const vi = re[b] * ci + im[b] * cr;
        re[b] = re[a] - vr; im[b] = im[a] - vi;
        re[a] += vr;        im[a] += vi;
        const ncr = cr * wr - ci * wi;
        ci = cr * wi + ci * wr;
        cr = ncr;
      }
    }
  }
}

/** Spectrum of one mic: Hann-windowed, mean-removed FFT magnitude of the
 *  newest pow2-sized stretch of the ring. Returns null until ~5s of data. */
function micSpectrum(m) {
  let n = 256;
  while (n * 2 <= Math.min(m.count, MIC_N)) n *= 2;
  if (n > Math.min(m.count, MIC_N)) return null;
  const re = new Float32Array(n), im = new Float32Array(n);
  const start = (m.head - n + MIC_N) % MIC_N;
  let mean = 0;
  for (let i = 0; i < n; i++) mean += m.buf[(start + i) % MIC_N];
  mean /= n;
  for (let i = 0; i < n; i++) {
    const w = 0.5 - 0.5 * Math.cos((2 * Math.PI * i) / (n - 1));   // Hann
    re[i] = (m.buf[(start + i) % MIC_N] - mean) * w;
  }
  fft(re, im);
  const mag = new Float32Array(n / 2);
  for (let i = 1; i < n / 2; i++) mag[i] = Math.hypot(re[i], im[i]);
  return { mag, n };
}

function drawSpectrum() {
  const c = $('spectrum');
  if (!mics.length) { c.style.display = 'none'; return; }
  c.style.display = 'block';
  const ctx = c.getContext('2d');
  const w = c.width, h = c.height;
  ctx.clearRect(0, 0, w, h);
  ctx.font = '11px ui-monospace, Consolas, monospace';

  const L = Module._lbm_get_char_length();
  const u0 = parseFloat($('vel').value);
  const dt = stepsPerFrame;          // sim steps per recorded sample
  const plotL = 42, plotR = w - 8, plotT = 22, plotB = h - 26;
  ctx.strokeStyle = 'rgba(255,255,255,0.25)';
  ctx.strokeRect(plotL, plotT, plotR - plotL, plotB - plotT);
  ctx.fillStyle = '#8d97a5';
  ctx.fillText('|p′(f)|  (log)', 8, 14);

  let anyData = false;
  mics.forEach((m, mi) => {
    const spec = micSpectrum(m);
    const color = MIC_COLORS[mi];
    if (!spec) {
      ctx.fillStyle = color;
      ctx.fillText(`M${mi + 1} collecting…`, plotL + 6 + mi * 90, h - 8);
      return;
    }
    anyData = true;
    const { mag, n } = spec;
    const nShow = Math.min(mag.length, 256);     // low-frequency band
    let hi = 1e-12;
    for (let i = 1; i < nShow; i++) if (mag[i] > hi) hi = mag[i];
    const lo = hi * 1e-4;
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let i = 1; i < nShow; i++) {
      const x = plotL + ((i - 1) / (nShow - 2)) * (plotR - plotL);
      const t = Math.max(0, Math.log(Math.max(mag[i], lo) / lo)
                            / Math.log(hi / lo));
      const y = plotB - t * (plotB - plotT);
      if (i === 1) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    // Dominant peak -> Strouhal number St = f L / u0 (f in cycles/step).
    let pk = 2;
    for (let i = 3; i < nShow; i++) if (mag[i] > mag[pk]) pk = i;
    const f = pk / (n * dt);
    const st = u0 > 0 ? (f * L) / u0 : 0;
    ctx.fillStyle = color;
    ctx.fillText(`M${mi + 1} St=${st.toFixed(2)}`, plotL + 6 + mi * 90, h - 8);
  });
  if (anyData) {
    ctx.fillStyle = '#8d97a5';
    ctx.textAlign = 'right';
    ctx.fillText(`f → 0 … ${(256 / (MIC_N * dt)).toExponential(1)} cyc/step`,
                 plotR, 14);
    ctx.textAlign = 'left';
  }
}

function drawMicMarkers(ctx, w, h) {
  mics.forEach((m, mi) => {
    const px = m.u * w, py = h - m.v * h;
    ctx.strokeStyle = MIC_COLORS[mi];
    ctx.fillStyle = MIC_COLORS[mi];
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(px, py, 5, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(px, py, 1.6, 0, 2 * Math.PI);
    ctx.fill();
    ctx.font = `${Math.max(10, w / 90)}px ui-monospace, Consolas, monospace`;
    ctx.fillText(`M${mi + 1}`, px + 7, py - 6);
  });
}

/** Dashed translucent contour of the armed preset, following the cursor. */
function drawGhost(ctx, w, h) {
  const angle = parseFloat($('angle').value);
  const pts = outlinePreset(placing, Ny, angle);
  if (!pts.length) return;
  const cx = ghostPos[0] * Nx, cy = ghostPos[1] * Ny;
  ctx.beginPath();
  for (let i = 0; i < pts.length; i++) {
    const px = ((cx + pts[i][0]) / Nx) * w;
    const py = h - ((cy + pts[i][1]) / Ny) * h;
    if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
  }
  ctx.closePath();
  ctx.fillStyle = 'rgba(255, 255, 255, 0.16)';
  ctx.fill();
  ctx.strokeStyle = 'rgba(255, 255, 255, 0.75)';
  ctx.lineWidth = Math.max(1, w / 1500);
  ctx.setLineDash([6, 4]);
  ctx.stroke();
  ctx.setLineDash([]);
}

function drawOverlay() {
  const c = $('overlay');
  const dpr = window.devicePixelRatio || 1;
  const w = Math.max(1, Math.round(c.clientWidth * dpr));
  const h = Math.max(1, Math.round(c.clientHeight * dpr));
  if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
  if (!octx) octx = c.getContext('2d');
  octx.clearRect(0, 0, w, h);
  if (overlayMode === 'streamlines') drawStreamlines(octx, w, h);
  else if (overlayMode === 'arrows') drawArrows(octx, w, h);
  else if (overlayMode === 'particles') drawParticles(octx, w, h);
  if (mics.length) drawMicMarkers(octx, w, h);
  if (placing && ghostPos) drawGhost(octx, w, h);
}

/* ------------------------------------------------------------------ */
/* Render loop                                                          */
/* ------------------------------------------------------------------ */

function fieldBounds(view) {
  if (scaleMode === 'fixed') {
    return [parseFloat($('scale-min').value), parseFloat($('scale-max').value)];
  }
  // Dynamic: single pass over the float view (cheap: <= 240k elements).
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < view.length; i++) {
    const v = view[i];
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  if (!isFinite(lo) || !isFinite(hi) || lo === hi) { lo = 0; hi = 1; }
  if (fieldMode === 'vorticity' || fieldMode === 'dilatation') {
    // Symmetric bounds keep the diverging colormap centered on zero.
    const a = Math.max(Math.abs(lo), Math.abs(hi));
    return [-a, a];
  }
  return [lo, hi];
}

function updateHUD() {
  Module._lbm_get_forces(forcesPtr, forcesPtr + 4);
  const Cl = Module.getValue(forcesPtr, 'float');
  const Cd = Module.getValue(forcesPtr + 4, 'float');
  const tau = parseFloat($('tau').value);
  const re = derivedRe();
  $('hud').textContent =
    `Re = ${re.toFixed(0)}    τ = ${tau.toFixed(3)}    step ${stepCount}\n` +
    `Cl = ${Cl >= 0 ? '+' : ''}${Cl.toFixed(3)}   Cd = ${Cd.toFixed(3)}\n` +
    `${msFrame > 0 ? (1000 / msFrame).toFixed(0) : '--'} fps · ` +
    `${threadCount} thread${threadCount > 1 ? 's' : ''}`;
  $('re-val').textContent = `Re ≈ ${re.toFixed(0)}`;
}

function tick(now) {
  for (let s = 0; s < stepsPerFrame; s++) Module._lbm_step();
  stepCount += stepsPerFrame;
  Module._lbm_compute_fields();

  ensureViews();
  if (overlayMode === 'particles') advectParticles();
  if (mics.length) {
    recordMics();
    if (frame % 30 === 0) drawSpectrum();
  }
  const view = currentFieldView();
  const [lo, hi] = fieldBounds(view);
  renderer.uploadField(view);
  if (obstacleDirty) {
    renderer.uploadObstacle(views.obstacle);
    obstacleDirty = false;
  }
  renderer.draw(lo, hi);
  drawOverlay();

  if (lastFrameT) {
    const dt = now - lastFrameT;
    msFrame = msFrame ? msFrame + 0.1 * (dt - msFrame) : dt;
  }
  lastFrameT = now;
  if (frame % 10 === 0) updateHUD();
  frame++;
  requestAnimationFrame(tick);
}

/* ------------------------------------------------------------------ */
/* UI wiring                                                            */
/* ------------------------------------------------------------------ */

function setupUI() {
  // Live parameter sliders.
  for (const id of ['vel', 'angle', 'tau']) {
    $(id).addEventListener('input', pushParams);
  }

  // Grid resolution: resize on release only (full realloc).
  $('res').addEventListener('input', () => {
    $('res-val').textContent = `${$('res').value} × …`;
  });
  $('res').addEventListener('change', () => resizeGrid(parseInt($('res').value, 10)));

  // Domain length: re-letterbox live, regrid on release (full realloc).
  $('domlen').addEventListener('input', () => {
    domainAspect = parseFloat($('domlen').value);
    $('domlen-val').textContent = `${domainAspect.toFixed(1)} : 1`;
    layoutCanvases();
  });
  $('domlen').addEventListener('change', () => resizeGrid(Nx));
  window.addEventListener('resize', layoutCanvases);

  $('spin').addEventListener('input', () => {
    $('spin-val').textContent = parseFloat($('spin').value).toFixed(3);
  });

  $('steps').addEventListener('input', () => {
    stepsPerFrame = parseInt($('steps').value, 10);
    $('steps-val').textContent = String(stepsPerFrame);
  });

  // Visualization mode (segmented buttons).
  for (const btn of document.querySelectorAll('#viz-mode button')) {
    btn.addEventListener('click', () => {
      document.querySelectorAll('#viz-mode button')
        .forEach((b) => b.classList.toggle('active', b === btn));
      fieldMode = btn.dataset.mode;
      Module._lbm_set_active_field(FIELD_MASKS[fieldMode]);
      const cmap = DEFAULT_CMAP[fieldMode];
      $('cmap').value = cmap;
      renderer.setColormap(COLORMAPS[cmap]);
      if (scaleMode === 'fixed') {
        $('scale-min').value = FIXED_DEFAULTS[fieldMode].min;
        $('scale-max').value = FIXED_DEFAULTS[fieldMode].max;
      }
    });
  }

  $('cmap').addEventListener('change', () => {
    renderer.setColormap(COLORMAPS[$('cmap').value]);
  });

  $('scale-mode').addEventListener('change', () => {
    scaleMode = $('scale-mode').checked ? 'fixed' : 'dynamic';
    $('fixed-scale-row').style.display = scaleMode === 'fixed' ? '' : 'none';
    if (scaleMode === 'fixed') {
      $('scale-min').value = FIXED_DEFAULTS[fieldMode].min;
      $('scale-max').value = FIXED_DEFAULTS[fieldMode].max;
    }
  });

  $('walls').addEventListener('change', () => {
    Module._lbm_set_wall_mode($('walls').value === 'noslip' ? 1 : 0);
  });

  $('les').addEventListener('change', () => {
    Module._lbm_set_les($('les').checked ? 1 : 0, 0.18);
  });

  $('reg').addEventListener('change', () => {
    Module._lbm_set_regularized($('reg').checked ? 1 : 0);
  });

  $('profile').addEventListener('change', () => {
    Module._lbm_set_flow_profile(parseInt($('profile').value, 10));
    Module._lbm_reset_flow();
    stepCount = 0;
  });
  $('reset-flow').addEventListener('click', () => {
    Module._lbm_reset_flow();
    stepCount = 0;
  });

  $('flow-overlay').addEventListener('change', () => {
    overlayMode = $('flow-overlay').value;
    $('particle-opts').style.display = overlayMode === 'particles' ? '' : 'none';
  });

  $('psize').addEventListener('input', () => {
    particleSize = parseInt($('psize').value, 10);
    $('psize-val').textContent = String(particleSize);
  });
  $('pcolor').addEventListener('input', () => {
    particleColor = $('pcolor').value;
  });

  // Draw / erase toggle.
  for (const id of ['mode-draw', 'mode-erase']) {
    $(id).addEventListener('click', () => {
      paintMode = id === 'mode-draw' ? 'draw' : 'erase';
      $('mode-draw').classList.toggle('active', paintMode === 'draw');
      $('mode-erase').classList.toggle('active', paintMode === 'erase');
    });
  }

  // Preset obstacle buttons arm click-to-place (click again to cancel).
  for (const btn of document.querySelectorAll('[data-preset]')) {
    btn.addEventListener('click', () =>
      setPlacing(placing === btn.dataset.preset ? null : btn.dataset.preset));
  }
  $('clear-obstacles').addEventListener('click', clearObstacles);

  // Microphones: arm click-to-place; clear empties the list.
  $('mic-add').addEventListener('click', () => {
    placingMic = !placingMic && mics.length < MIC_MAX;
    setPlacing(null);
    $('mic-add').classList.toggle('placing', placingMic);
  });
  $('mic-clear').addEventListener('click', () => {
    mics.length = 0;
    $('spectrum').style.display = 'none';
  });

  // Collapsible panel.
  $('panel-toggle').addEventListener('click', () => {
    $('panel').classList.toggle('collapsed');
  });
}

/* ------------------------------------------------------------------ */
/* Boot                                                                 */
/* ------------------------------------------------------------------ */

async function boot() {
  const createLBM = await loadEngine();
  Module = await createLBM();
  forcesPtr = Module._malloc(8);

  const canvas = $('glcanvas');
  layoutCanvases();
  renderer = new Renderer(canvas);
  renderer.resizeToDisplay();

  // Initial grid: 300 on the long axis, height from the domain aspect.
  Ny = Math.min(512, Math.max(48, Math.round(Nx / domainAspect)));
  if (Module._lbm_init(Nx, Ny) !== 0) {
    throw new Error('lbm_init failed');
  }
  $('res-val').textContent = `${Nx} × ${Ny}`;
  initParticles();

  pushParams();
  Module._lbm_set_wall_mode(0);          // free-slip tunnel walls
  Module._lbm_set_les(1, 0.18);          // Smagorinsky on by default
  Module._lbm_set_regularized(1);        // Hermite-regularized collision
  Module._lbm_set_active_field(FIELD_MASKS[fieldMode]);
  if (engineMT) {
    Module._lbm_set_threads(Math.min(8, navigator.hardwareConcurrency || 4));
  }
  threadCount = Module._lbm_get_threads();
  Module._lbm_add_preset(0, 0.18, 0.5, 0.25, 0);  // starter cylinder, upstream
  cShapes.push({ name: 'cylinder', omega: 0 });
  rebuildCylinderList();

  renderer.setGrid(Nx, Ny);
  renderer.setColormap(COLORMAPS[DEFAULT_CMAP[fieldMode]]);
  refreshViews();
  renderer.uploadObstacle(views.obstacle);

  setupUI();
  setupPainting(canvas);
  requestAnimationFrame(tick);
}

boot().catch((err) => {
  console.error(err);
  const hud = $('hud');
  if (hud) hud.textContent = `Failed to start: ${err.message}\n` +
    'Run build.sh first and serve over HTTP (wasm cannot load from file://).';
});
