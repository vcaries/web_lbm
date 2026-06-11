#!/usr/bin/env bash
# build.sh — compile the LBM engine to WebAssembly with Emscripten.
#
# Usage:
#   bash build.sh           # BGK collision (default)
#   bash build.sh --mrt     # MRT collision (-DUSE_MRT), more stable at high Re
#
# If emcc is not on PATH but Docker is available, the official
# emscripten/emsdk image is used instead (same flags, same output).
set -euo pipefail
cd "$(dirname "$0")"

EXTRA_CFLAGS=""
if [[ "${1:-}" == "--mrt" ]]; then
  EXTRA_CFLAGS="-DUSE_MRT"
  echo "Building with MRT collision operator."
fi

mkdir -p dist

# Every emcc flag, explained:
#
#   -O3                  Maximum optimization. The collision/streaming loops
#                        dominate runtime; -O3 enables vectorization and
#                        aggressive inlining of feq().
#   -msimd128            Emit WebAssembly SIMD (128-bit). Lets LLVM auto-
#                        vectorize the per-cell loops; supported by all
#                        current browsers.
#   -ffast-math          Relaxed FP semantics (no NaN/Inf special-casing in
#                        codegen). Safe here: the engine clamps velocity and
#                        tau so fields stay finite by construction.
#   -s WASM=1            Emit a .wasm binary (not asm.js).
#   -s EXPORTED_FUNCTIONS
#                        The exact C symbols kept callable from JS. Anything
#                        not listed is eligible for dead-code elimination.
#                        _malloc/_free are included so JS can allocate the
#                        2-float scratch buffer lbm_get_forces() writes into.
#   -s EXPORTED_RUNTIME_METHODS
#                        getValue/setValue for typed reads of out-params, and
#                        wasmMemory so JS can build zero-copy TypedArray
#                        views over the engine's heap (and rebuild them after
#                        memory growth, which detaches old ArrayBuffers).
#   -s ALLOW_MEMORY_GROWTH=1
#                        REQUIRED: lbm_resize() frees and reallocates every
#                        field array at runtime. Without growth, a resize to
#                        a larger grid than the initial heap accommodates
#                        would make malloc fail. With it, the wasm memory
#                        grows transparently (JS views are rebuilt on demand).
#   -s INITIAL_MEMORY=33554432
#                        32 MB up front — enough for the largest grid
#                        (600x300 needs ~18 MB), so growth events are rare.
#   -s MAXIMUM_MEMORY=268435456
#                        Cap growth at 256 MB so a runaway allocation can
#                        never exhaust the host machine.
#   -s MODULARIZE=1      Wrap the runtime in a factory function instead of
#                        polluting globals: allows multiple independent
#                        engine instances on one page and clean namespacing.
#   -s EXPORT_ES6=1      Emit the factory as an ES module default export, so
#                        the vanilla-JS app can `import createLBM from ...`
#                        with no bundler.
#   -s ENVIRONMENT=web   Drop Node/worker/shell support code from the glue:
#                        smaller output, and GitHub Pages only serves web.
#   -s FILESYSTEM=0      The engine does no I/O; omitting Emscripten's
#                        virtual FS shrinks the JS glue considerably.
#   -o dist/lbm_engine.js
#                        Emits dist/lbm_engine.js (ES module glue) and
#                        dist/lbm_engine.wasm side by side. The glue locates
#                        the .wasm relative to its own URL, so the pair can
#                        be served from any static path.

EMCC_ARGS=(
  src/lbm.c
  -O3 -msimd128 -ffast-math
  -Wall -Wextra
  $EXTRA_CFLAGS
  -s WASM=1
  -s "EXPORTED_FUNCTIONS=[\"_lbm_init\",\"_lbm_step\",\"_lbm_set_params\",\"_lbm_compute_fields\",\"_lbm_set_obstacle\",\"_lbm_add_preset\",\"_lbm_get_forces\",\"_lbm_resize\",\"_lbm_get_ux_ptr\",\"_lbm_get_uy_ptr\",\"_lbm_get_rho_ptr\",\"_lbm_get_speed_ptr\",\"_lbm_get_vorticity_ptr\",\"_lbm_get_schlieren_ptr\",\"_lbm_get_dilatation_ptr\",\"_lbm_get_obstacle_ptr\",\"_lbm_set_wall_mode\",\"_lbm_set_les\",\"_lbm_set_regularized\",\"_lbm_set_flow_profile\",\"_lbm_reset_flow\",\"_lbm_set_shape_omega\",\"_lbm_clear_obstacles\",\"_lbm_get_char_length\",\"_malloc\",\"_free\"]"
  -s "EXPORTED_RUNTIME_METHODS=[\"getValue\",\"setValue\",\"wasmMemory\"]"
  -s ALLOW_MEMORY_GROWTH=1
  -s INITIAL_MEMORY=33554432
  -s MAXIMUM_MEMORY=268435456
  -s MODULARIZE=1
  -s EXPORT_ES6=1
  -s ENVIRONMENT=web
  -s FILESYSTEM=0
  -o dist/lbm_engine.js
)

if command -v emcc >/dev/null 2>&1; then
  emcc "${EMCC_ARGS[@]}"
elif command -v docker >/dev/null 2>&1; then
  echo "emcc not found — falling back to the emscripten/emsdk Docker image."
  # pwd -W yields a Windows-style path under Git Bash (no-op elsewhere);
  # MSYS_NO_PATHCONV stops Git Bash from rewriting the /src mount point.
  HOST_DIR="$(pwd -W 2>/dev/null || pwd)"
  MSYS_NO_PATHCONV=1 docker run --rm -v "${HOST_DIR}:/src" -w /src \
    emscripten/emsdk emcc "${EMCC_ARGS[@]}"
else
  echo "ERROR: neither emcc nor docker found. Install the Emscripten SDK" >&2
  echo "  https://emscripten.org/docs/getting_started/downloads.html" >&2
  exit 1
fi

echo "Build OK:"
ls -l dist/lbm_engine.js dist/lbm_engine.wasm
