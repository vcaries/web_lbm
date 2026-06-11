@echo off
rem build.bat — compile the LBM engine to WebAssembly on Windows.
rem
rem Usage:
rem   build.bat          (BGK collision, default)
rem   build.bat --mrt    (MRT collision, -DUSE_MRT, more stable at high Re)
rem
rem Flag-by-flag rationale lives in build.sh (the canonical build script);
rem the flags here are identical. Highlights:
rem   ALLOW_MEMORY_GROWTH=1  required because lbm_resize() reallocates all
rem                          field arrays at runtime; the wasm heap must be
rem                          able to grow past its initial size.
rem   MODULARIZE=1           wraps the runtime in a factory function for
rem                          clean namespacing / multiple instances.
rem   EXPORT_ES6=1           default-exports that factory as an ES module.
setlocal

set EXTRA=
if "%~1"=="--mrt" set EXTRA=-DUSE_MRT
if "%~1"=="--mrt" echo Building with MRT collision operator.

if not exist dist mkdir dist

rem Same settings as build.sh, written in emcc's bracket-less list syntax
rem (name1,name2,...) so no quoting is needed — cmd, docker.exe and emcc.bat
rem all parse quote-free arguments identically.
set EXPORTS=EXPORTED_FUNCTIONS=_lbm_init,_lbm_step,_lbm_set_params,_lbm_compute_fields,_lbm_set_obstacle,_lbm_add_preset,_lbm_get_forces,_lbm_resize,_lbm_get_ux_ptr,_lbm_get_uy_ptr,_lbm_get_rho_ptr,_lbm_get_speed_ptr,_lbm_get_vorticity_ptr,_lbm_get_schlieren_ptr,_lbm_get_dilatation_ptr,_lbm_get_obstacle_ptr,_lbm_set_wall_mode,_lbm_set_les,_lbm_set_regularized,_lbm_set_flow_profile,_lbm_reset_flow,_lbm_set_shape_omega,_lbm_clear_obstacles,_lbm_get_char_length,_malloc,_free
set RUNTIME=EXPORTED_RUNTIME_METHODS=getValue,setValue,wasmMemory

where emcc >nul 2>nul
if errorlevel 1 goto try_docker

call emcc src/lbm.c -O3 -msimd128 -ffast-math -Wall -Wextra %EXTRA% ^
  -s WASM=1 ^
  -s %EXPORTS% ^
  -s %RUNTIME% ^
  -s ALLOW_MEMORY_GROWTH=1 ^
  -s INITIAL_MEMORY=33554432 ^
  -s MAXIMUM_MEMORY=268435456 ^
  -s MODULARIZE=1 ^
  -s EXPORT_ES6=1 ^
  -s ENVIRONMENT=web ^
  -s FILESYSTEM=0 ^
  -o dist/lbm_engine.js
if errorlevel 1 exit /b 1
goto done

:try_docker
where docker >nul 2>nul
if errorlevel 1 goto no_tools
echo emcc not found - falling back to the emscripten/emsdk Docker image.
docker run --rm -v "%cd%:/src" -w /src emscripten/emsdk emcc src/lbm.c ^
  -O3 -msimd128 -ffast-math -Wall -Wextra %EXTRA% ^
  -s WASM=1 ^
  -s %EXPORTS% ^
  -s %RUNTIME% ^
  -s ALLOW_MEMORY_GROWTH=1 ^
  -s INITIAL_MEMORY=33554432 ^
  -s MAXIMUM_MEMORY=268435456 ^
  -s MODULARIZE=1 ^
  -s EXPORT_ES6=1 ^
  -s ENVIRONMENT=web ^
  -s FILESYSTEM=0 ^
  -o dist/lbm_engine.js
if errorlevel 1 exit /b 1
goto done

:no_tools
echo ERROR: neither emcc nor docker found. Install the Emscripten SDK:
echo   https://emscripten.org/docs/getting_started/downloads.html
exit /b 1

:done
echo Build OK:
dir dist\lbm_engine.js dist\lbm_engine.wasm
endlocal
