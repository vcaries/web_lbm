/*
 * colormaps.js — 256x1 RGBA lookup tables for the fragment shader.
 *
 * Each LUT is a Uint8Array of 256*4 bytes, built once at module load by
 * piecewise-linear interpolation between anchor colors. The anchors for
 * Viridis and Inferno are sampled from the reference matplotlib maps;
 * Jet is the classic MATLAB rainbow; Coolwarm is Moreland's diverging
 * blue-white-red, ideal for signed fields such as vorticity.
 */

const N = 256;

/** Interpolate [t, r, g, b] anchor rows (t in [0,1]) into a 256x1 RGBA LUT. */
function buildLUT(anchors) {
  const lut = new Uint8Array(N * 4);
  let seg = 0;
  for (let i = 0; i < N; i++) {
    const t = i / (N - 1);
    while (seg < anchors.length - 2 && t > anchors[seg + 1][0]) seg++;
    const [t0, r0, g0, b0] = anchors[seg];
    const [t1, r1, g1, b1] = anchors[seg + 1];
    const u = t1 > t0 ? Math.min(1, Math.max(0, (t - t0) / (t1 - t0))) : 0;
    lut[i * 4 + 0] = Math.round(r0 + (r1 - r0) * u);
    lut[i * 4 + 1] = Math.round(g0 + (g1 - g0) * u);
    lut[i * 4 + 2] = Math.round(b0 + (b1 - b0) * u);
    lut[i * 4 + 3] = 255;
  }
  return lut;
}

export const COLORMAPS = {
  JET: buildLUT([
    [0.000, 0, 0, 131], [0.125, 0, 60, 170], [0.375, 5, 255, 255],
    [0.625, 255, 255, 0], [0.875, 250, 50, 0], [1.000, 128, 0, 0],
  ]),
  INFERNO: buildLUT([
    [0.000, 0, 0, 4], [0.130, 31, 12, 72], [0.250, 85, 15, 109],
    [0.380, 136, 34, 106], [0.500, 186, 54, 85], [0.630, 227, 89, 51],
    [0.750, 249, 140, 10], [0.880, 249, 201, 50], [1.000, 252, 255, 164],
  ]),
  COOLWARM: buildLUT([
    [0.000, 59, 76, 192], [0.250, 124, 159, 249], [0.500, 221, 221, 221],
    [0.750, 245, 156, 125], [1.000, 180, 4, 38],
  ]),
  VIRIDIS: buildLUT([
    [0.000, 68, 1, 84], [0.130, 71, 44, 122], [0.250, 59, 81, 139],
    [0.380, 44, 113, 142], [0.500, 33, 144, 141], [0.630, 39, 173, 129],
    [0.750, 92, 200, 99], [0.880, 170, 220, 50], [1.000, 253, 231, 37],
  ]),
  // Plain luminance ramp: the schlieren look (bright wave fronts and
  // shear layers on black, like a negative of the knife-edge photograph).
  GRAY: buildLUT([
    [0.000, 0, 0, 0], [1.000, 255, 255, 255],
  ]),
};
