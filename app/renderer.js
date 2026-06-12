/*
 * renderer.js — WebGL2 renderer for the LBM fields.
 *
 * Data path (zero-copy): the C engine's field arrays live inside
 * WebAssembly.Memory; main.js wraps them in Float32Array/Uint8Array views
 * and this module hands those views straight to texSubImage2D. No staging
 * buffers, no per-frame allocation.
 *
 * Textures:
 *   uField    R32F,  NEAREST — the scalar field (speed / rho / vorticity),
 *                              re-uploaded every frame.
 *   uObstacle R8UI,  NEAREST — the obstacle mask, re-uploaded only when it
 *                              changes (usampler2D in the shader).
 *   uLUT      RGBA8, LINEAR  — 256x1 colormap lookup table.
 *
 * The fragment shader normalizes the scalar with (v - min) / (max - min),
 * clamps, and samples the LUT; obstacle cells are drawn flat gray.
 */

const VS = `#version 300 es
out vec2 vUV;
void main() {
  // Fullscreen triangle: 3 vertices, no buffers needed.
  vec2 pos = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
  vUV = pos;
  gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);
}`;

const FS = `#version 300 es
precision highp float;
precision highp int;
uniform sampler2D uField;
uniform highp usampler2D uObstacle;
uniform sampler2D uLUT;
uniform float uMin;
uniform float uMax;
in vec2 vUV;
out vec4 outColor;
void main() {
  uint solid = texture(uObstacle, vUV).r;
  if (solid > 0u) {
    outColor = vec4(0.32, 0.33, 0.36, 1.0);   // obstacle overlay
    return;
  }
  float v = texture(uField, vUV).r;
  float t = clamp((v - uMin) / max(uMax - uMin, 1e-20), 0.0, 1.0);
  outColor = texture(uLUT, vec2(t, 0.5));
}`;

function compile(gl, type, src) {
  const sh = gl.createShader(type);
  gl.shaderSource(sh, src);
  gl.compileShader(sh);
  if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
    throw new Error('Shader compile error: ' + gl.getShaderInfoLog(sh));
  }
  return sh;
}

export class Renderer {
  constructor(canvas) {
    const gl = canvas.getContext('webgl2', {
      antialias: false, depth: false, stencil: false, alpha: false,
    });
    if (!gl) throw new Error('WebGL2 is required but not available.');
    this.gl = gl;
    this.canvas = canvas;
    this.Nx = 0;
    this.Ny = 0;

    const prog = gl.createProgram();
    gl.attachShader(prog, compile(gl, gl.VERTEX_SHADER, VS));
    gl.attachShader(prog, compile(gl, gl.FRAGMENT_SHADER, FS));
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
      throw new Error('Program link error: ' + gl.getProgramInfoLog(prog));
    }
    this.prog = prog;
    gl.useProgram(prog);
    this.uMin = gl.getUniformLocation(prog, 'uMin');
    this.uMax = gl.getUniformLocation(prog, 'uMax');
    gl.uniform1i(gl.getUniformLocation(prog, 'uField'), 0);
    gl.uniform1i(gl.getUniformLocation(prog, 'uObstacle'), 1);
    gl.uniform1i(gl.getUniformLocation(prog, 'uLUT'), 2);

    // An empty VAO must still be bound to issue draw calls.
    gl.bindVertexArray(gl.createVertexArray());

    // Field rows are tightly packed floats; obstacle rows are single
    // bytes, so default 4-byte row alignment would corrupt odd widths.
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);

    this.fieldTex = null;
    this.obstacleTex = null;

    this.lutTex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, this.lutTex);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA8, 256, 1);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  }

  /** (Re)create grid-sized textures. Called on init and after lbm_resize. */
  setGrid(Nx, Ny) {
    const gl = this.gl;
    this.Nx = Nx;
    this.Ny = Ny;

    if (this.fieldTex) gl.deleteTexture(this.fieldTex);
    this.fieldTex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.fieldTex);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.R32F, Nx, Ny);
    // Float textures are not filterable without an extension: use NEAREST.
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

    if (this.obstacleTex) gl.deleteTexture(this.obstacleTex);
    this.obstacleTex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, this.obstacleTex);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.R8UI, Nx, Ny);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  }

  /** Upload a 256x1 RGBA Uint8Array LUT. */
  setColormap(lut) {
    const gl = this.gl;
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, this.lutTex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, 256, 1,
                     gl.RGBA, gl.UNSIGNED_BYTE, lut);
  }

  /** With the -pthread engine the wasm heap is a SharedArrayBuffer, and
   *  not every browser accepts SAB-backed views in texSubImage2D — stage
   *  through a small ordinary buffer in that case (the one exception to
   *  the zero-copy data path; ~0.05 ms for a 300x150 field). */
  _staged(view, Ctor, slot) {
    if (typeof SharedArrayBuffer === 'undefined'
        || !(view.buffer instanceof SharedArrayBuffer)) return view;
    if (!this[slot] || this[slot].length !== view.length) {
      this[slot] = new Ctor(view.length);
    }
    this[slot].set(view);
    return this[slot];
  }

  /** Per-frame scalar upload, straight from the wasm-memory view. */
  uploadField(view) {
    const gl = this.gl;
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.fieldTex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, this.Nx, this.Ny,
                     gl.RED, gl.FLOAT,
                     this._staged(view, Float32Array, '_stageF'));
  }

  /** Obstacle mask upload — only called when the mask changed. */
  uploadObstacle(view) {
    const gl = this.gl;
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, this.obstacleTex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, this.Nx, this.Ny,
                     gl.RED_INTEGER, gl.UNSIGNED_BYTE,
                     this._staged(view, Uint8Array, '_stageO'));
  }

  /** Match the drawing buffer to CSS size * devicePixelRatio. */
  resizeToDisplay() {
    const c = this.canvas;
    const dpr = window.devicePixelRatio || 1;
    const w = Math.max(1, Math.round(c.clientWidth * dpr));
    const h = Math.max(1, Math.round(c.clientHeight * dpr));
    if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
  }

  draw(min, max) {
    const gl = this.gl;
    this.resizeToDisplay();
    gl.viewport(0, 0, this.canvas.width, this.canvas.height);
    gl.useProgram(this.prog);
    gl.uniform1f(this.uMin, min);
    gl.uniform1f(this.uMax, max);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }
}
