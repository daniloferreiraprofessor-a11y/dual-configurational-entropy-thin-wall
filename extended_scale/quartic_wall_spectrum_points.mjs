// Collect the two independent large-scale quartic-bubble points in one CSV.
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const directory = dirname(fileURLToPath(import.meta.url));
const cases = [
  ['0p0015', '0p00003125', 8],
  ['0p001', '0p000015625', 4],
];
const rows = cases.map(([eps, dk, kMax]) => {
  const base = `quartic_wall_bvp_profile_eps${eps}_L15_h0p000625.json`;
  const spectrum = `quartic_wall_spectrum_eps${eps}_L15_h0p000625_stride8_dk${dk}_k${kMax}.json`;
  const bvp = JSON.parse(readFileSync(join(directory, base), 'utf8'));
  const measure = JSON.parse(readFileSync(join(directory, spectrum), 'utf8'));
  return {
    epsilon: bvp.epsilon,
    R: bvp.R_phase,
    ell: bvp.ell,
    Lambda_3: measure.Lambda_3,
    Sx_3: measure.Sx_3,
    Sk_3: measure.Sk_3,
    U_3: measure.U_3,
    H: measure.H,
    derrick_residual: bvp.derrick_residual,
    wall_window_L: bvp.L,
    bvp_h: bvp.h,
    quadrature_wall_h: measure.wall_h,
    spectral_dk: measure.spectral_dk,
    spectral_k_max: measure.kMax,
    profile_method: 'wall_centered_matching',
    spectrum_method: 'analytic_ball_plus_direct_wall_quadrature',
  };
});
const keys = Object.keys(rows[0]);
writeFileSync(join(directory, 'quartic_wall_spectrum_points.csv'),
  [keys.join(','), ...rows.map((row) => keys.map((key) => row[key]).join(','))]
    .join('\n') + '\n');
