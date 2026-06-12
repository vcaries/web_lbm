# web_lbm

2D interactive wind tunnel: a Lattice Boltzmann (D2Q9) engine in C compiled
to WebAssembly, rendered with WebGL2, vanilla ES6 UI. Static-host friendly
(GitHub Pages) — no JS frameworks, no runtime build step.

## Stack
- Engine: C (C99), compiled with Emscripten (`emcc`) to wasm + ES-module glue
- Frontend: vanilla ES6 modules, WebGL2 (no frameworks, no bundler)
- Build: `build.sh` / `build.bat` (fall back to the `emscripten/emsdk`
  Docker image when `emcc` is not on PATH)

## Commands
- Build: `bash build.sh` (or `build.bat`); `--mrt` flag selects the MRT
  collision operator. Produces **two** engine variants: `dist/lbm_engine.js`
  (single-thread) and `dist/lbm_engine_mt.js` (`-pthread`); the app
  feature-detects `crossOriginIsolated` and picks at runtime
  (`app/coi-serviceworker.min.js` provides the COOP/COEP headers on
  GitHub Pages)
- Run: `python -m http.server 8000` from repo root, open `http://localhost:8000/app/`
  (wasm cannot load from `file://`)
- Test: `emcc src/lbm.c -O2 -DLBM_TEST_MAIN -s ENVIRONMENT=node -s ALLOW_MEMORY_GROWTH=1 -o tmp_test/test.js && node tmp_test/test.js`
  (physics smoke tests embedded in `src/lbm.c` behind `LBM_TEST_MAIN`;
  `ALLOW_MEMORY_GROWTH` is required by the resize test; add `-DUSE_MRT` to
  test the MRT operator)

## Layout
- `src/lbm.c` — entire C engine: D2Q9 BGK (+MRT via `-DUSE_MRT`), Zou-He inlet,
  zero-gradient outlet, bounce-back/specular walls, momentum-exchange forces,
  Smagorinsky LES (`lbm_set_les`, on by default, Cs = 0.18),
  Hermite-regularized collision (`lbm_set_regularized`, on by default,
  BGK only) with bulk-viscosity damping of the acoustic stress
  (`BULK_DAMP`), absorbing sponge layers at inlet/outlet,
  400-step inlet soft-start after (re)init, inlet flow
  profiles (`lbm_set_flow_profile`: 0 uniform, 1 shear layer, 2 jet) with
  `lbm_reset_flow`, rotating cylinders via moving-wall bounce-back
  (`lbm_set_shape_omega(index, surface_speed)`, resolution-independent,
  replayed on resize), derived fields incl. schlieren `|∇ρ|` and
  dilatation `∇·u` (acoustics), analytic shape rasterizers,
  resize-with-obstacle-persistence
- `app/` — `index.html`, `main.js` (init/RAF/UI), `renderer.js` (WebGL2),
  `colormaps.js` (LUTs), `obstacles.js` (analytic NACA + presets), `style.css`
- `dist/` — emcc output (`lbm_engine.js` + `.wasm`), gitignored during dev

## Invariants (do not break)
- **Zero-copy data path**: JS reads engine fields via TypedArray views over
  `Module.wasmMemory.buffer` and uploads straight to WebGL textures. Never
  introduce intermediate copies; rebuild views when the buffer detaches
  (memory growth) or after `lbm_resize()`. One sanctioned exception: with
  the `-pthread` engine the heap is a SharedArrayBuffer, and texture
  uploads stage through a small ordinary buffer (`renderer.js _staged`)
  because not every browser accepts SAB-backed views in `texSubImage2D`.
- **Threading model**: the interior sweep and fixup pass are row-partitioned
  across a pthread pool (`lbm_set_threads`, one-shot) with barriers per
  step; every f_tmp slot has exactly one producer, so slices are race-free
  with no locks. Forces reduce in tid order — results are deterministic for
  a fixed thread count. Border/boundary/sponge passes stay on the main
  thread.
- **SoA layout** in C: `float *f[9]`, `f_tmp[9]`, flat `Nx*Ny` arrays, `idx = y*Nx + x`.
- Obstacle mask semantics: 0 fluid, 1 hand-painted (resampled on resize),
  2 preset (re-rasterized analytically from registries on resize).
- New C exports must be added to `EXPORTED_FUNCTIONS` in **both** build scripts.
- All asset paths stay relative (GitHub Pages, no rewrite rules).
- **The interior collide+stream sweep stays branchless** (`restrict` pointer
  copies, no obstacle checks, LES folded into the arithmetic): that is what
  lets LLVM auto-vectorize it (~3x). Solid-cell handling belongs in the
  border path and the fixup pass, never in the hot loop.

## Tooling & conventions (token-efficient workflow — keep this)
- **Navigate with serena, not full-file reads.** Use serena's symbol/reference tools to find and edit code; only read whole files when you genuinely need the full context. This is the main token saver.
- **Look up library/API usage with context7** rather than guessing or relying on memory — it returns version-correct docs.
- **Use GitHub via the github MCP** for repo/PR/issue operations.
- **Ask project-wide questions with `/graphify`** — build/refresh the knowledge graph for "how does X work across the codebase" questions.
- **Lean on superpowers skills** as appropriate: test-driven-development for new features, systematic-debugging for bugs, brainstorming/writing-plans before large changes, verification-before-completion before declaring done.
- Match existing code style, naming, and structure. Run tests/lint before claiming a change works.
