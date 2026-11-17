// Independent wall-centered BVP for the quartic critical bubble.
// The radius is a shooting parameter for two half-wall boundary-value
// problems with a fixed phase, rather than a central-field shooting datum.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const directory = path.dirname(fileURLToPath(import.meta.url));
const EPSILONS = process.argv[2]
  ? process.argv[2].split(',').map(Number)
  : [0.003, 0.0015];
const DEFAULT_CONFIGS = [
  { L: 10, h: 0.02 },
  { L: 12, h: 0.02 },
  { L: 15, h: 0.02 },
  { L: 12, h: 0.01 },
  { L: 15, h: 0.01 },
  { L: 20, h: 0.01 },
  { L: 15, h: 0.005 },
  { L: 15, h: 0.0025 },
  { L: 12, h: 0.0025 },
  { L: 20, h: 0.0025 },
  { L: 15, h: 0.00125 },
  { L: 15, h: 0.000625 },
];
const CONFIGS = process.argv[3]
  ? process.argv[3].split(',').map((entry) => {
    const [L, h] = entry.split(':').map(Number);
    if (!(L > 0 && h > 0)) throw new Error(`invalid L:h configuration: ${entry}`);
    return { L, h };
  })
  : DEFAULT_CONFIGS;

function W(phi, eps) {
  return 0.25 * phi ** 2 * (phi ** 2 - 4 * (1 + eps) * phi + 4);
}
function Wp(phi, eps) {
  return phi * (phi ** 2 - 3 * (1 + eps) * phi + 2);
}
function Wpp(phi, eps) {
  return 3 * phi ** 2 - 6 * (1 + eps) * phi + 2;
}
function vacua(eps) {
  const discriminant = Math.sqrt(9 * (1 + eps) ** 2 - 8);
  const barrier = (3 * (1 + eps) - discriminant) / 2;
  const trueVacuum = (3 * (1 + eps) + discriminant) / 2;
  return { barrier, trueVacuum };
}

function tridiagonal(lower, diagonal, upper, right) {
  const n = right.length;
  const c = new Float64Array(n);
  const d = new Float64Array(n);
  let pivot = diagonal[0];
  if (!Number.isFinite(pivot) || Math.abs(pivot) < 1e-14)
    throw new Error(`bad tridiagonal pivot 0: ${pivot}`);
  c[0] = upper[0] / pivot;
  d[0] = right[0] / pivot;
  for (let j = 1; j < n; j++) {
    pivot = diagonal[j] - lower[j] * c[j - 1];
    if (!Number.isFinite(pivot) || Math.abs(pivot) < 1e-14)
      throw new Error(`bad tridiagonal pivot ${j}: ${pivot}`);
    c[j] = j < n - 1 ? upper[j] / pivot : 0;
    d[j] = (right[j] - lower[j] * d[j - 1]) / pivot;
  }
  const solution = new Float64Array(n);
  solution[n - 1] = d[n - 1];
  for (let j = n - 2; j >= 0; j--)
    solution[j] = d[j] - c[j] * solution[j + 1];
  return solution;
}

function halfResidual(y, eps, R, L, h, side, trueVacuum, phase, withJacobian) {
  const n = y.length - 1;
  const f = new Float64Array(n + 1);
  const lower = withJacobian ? new Float64Array(n + 1) : null;
  const diagonal = withJacobian ? new Float64Array(n + 1) : null;
  const upper = withJacobian ? new Float64Array(n + 1) : null;
  const mTrue = Math.sqrt(Wpp(trueVacuum, eps));
  const sLeft = -L;
  const rLeft = R - L;
  const coth = 1 + 2 / Math.expm1(2 * mTrue * rLeft);
  const a = mTrue * coth - 1 / rLeft;
  const rRight = R + L;
  const b = Math.SQRT2 + 1 / rRight;

  if (side === 'left') {
    const p = 2 / rLeft;
    const difference = y[0] - trueVacuum;
    f[0] = (y[1] - y[0]) / h - a * difference
      - 0.5 * h * (Wp(y[0], eps) - p * a * difference);
    if (withJacobian) {
      diagonal[0] = -1 / h - a - 0.5 * h * (Wpp(y[0], eps) - p * a);
      upper[0] = 1 / h;
    }
  } else {
    f[0] = y[0] - phase;
    if (withJacobian) diagonal[0] = 1;
  }

  for (let j = 1; j < n; j++) {
    const s = side === 'left' ? sLeft + j * h : j * h;
    const r = R + s;
    const low = 1 / h ** 2 - 1 / (r * h);
    const high = 1 / h ** 2 + 1 / (r * h);
    const forward = y[j + 1] - y[j];
    const backward = y[j] - y[j - 1];
    f[j] = (forward - backward) / h ** 2
      + (forward + backward) / (r * h) - Wp(y[j], eps);
    if (withJacobian) {
      lower[j] = low;
      diagonal[j] = -2 / h ** 2 - Wpp(y[j], eps);
      upper[j] = high;
    }
  }

  if (side === 'left') {
    f[n] = y[n] - phase;
    if (withJacobian) diagonal[n] = 1;
  } else {
    const p = 2 / rRight;
    f[n] = (y[n] - y[n - 1]) / h + b * y[n]
      + 0.5 * h * (Wp(y[n], eps) + p * b * y[n]);
    if (withJacobian) {
      lower[n] = -1 / h;
      diagonal[n] = 1 / h + b
        + 0.5 * h * (Wpp(y[n], eps) + p * b);
    }
  }

  let norm = 0;
  for (const value of f) norm = Math.max(norm, Math.abs(value));
  return { f, lower, diagonal, upper, norm };
}

function solveHalf(eps, R, L, h, side, trueVacuum, phase) {
  const n = Math.round(L / h);
  if (Math.abs(n * h - L) > 1e-10) throw new Error('L/h must be integer');
  const shift = -Math.log(Math.SQRT2 - 1) / Math.SQRT2;
  let y = Float64Array.from({ length: n + 1 }, (_, j) => {
    const s = side === 'left' ? -L + j * h : j * h;
    return trueVacuum / (1 + Math.exp(Math.SQRT2 * (s - shift)));
  });
  if (side === 'left') y[n] = phase;
  else y[0] = phase;

  let current = halfResidual(y, eps, R, L, h, side, trueVacuum, phase, false);
  let iterations = 0;
  for (; iterations < 45; iterations++) {
    if (current.norm < 1e-9) break;
    const jac = halfResidual(y, eps, R, L, h, side, trueVacuum, phase, true);
    const correction = tridiagonal(
      jac.lower, jac.diagonal, jac.upper,
      Float64Array.from(jac.f, (value) => -value),
    );
    let accepted = false;
    for (let factor = 1; factor >= 1 / 4096; factor /= 2) {
      const trial = Float64Array.from(y, (value, j) => value + factor * correction[j]);
      if (trial.some((value) => value < -1e-10 || value > trueVacuum * (1 + 1e-8)))
        continue;
      const next = halfResidual(trial, eps, R, L, h, side, trueVacuum, phase, false);
      if (next.norm < current.norm) {
        y = trial;
        current = next;
        accepted = true;
        break;
      }
    }
    if (!accepted) break;
  }
  if (current.norm > 3e-8)
    throw new Error(`${side} half-wall did not converge at R=${R}: residual=${current.norm}`);
  return { y, norm: current.norm, iterations };
}

function matching(eps, R, L, h, trueVacuum, phase) {
  const left = solveHalf(eps, R, L, h, 'left', trueVacuum, phase);
  const right = solveHalf(eps, R, L, h, 'right', trueVacuum, phase);
  const m = left.y.length - 1;
  const leftDerivative = (3 * left.y[m] - 4 * left.y[m - 1] + left.y[m - 2]) / (2 * h);
  const rightDerivative = (-3 * right.y[0] + 4 * right.y[1] - right.y[2]) / (2 * h);
  return { mismatch: leftDerivative - rightDerivative, leftDerivative,
    rightDerivative, left, right };
}

function crossing(s, y, target) {
  for (let j = 1; j < y.length; j++) {
    if (y[j] <= target) {
      const t = (y[j - 1] - target) / (y[j - 1] - y[j]);
      return s[j - 1] + t * (s[j] - s[j - 1]);
    }
  }
  throw new Error(`level ${target} not crossed`);
}

function virialFromWall(eps, R, L, h, phi, trueVacuum) {
  const n = phi.length - 1;
  const mTrue = Math.sqrt(Wpp(trueVacuum, eps));
  const rLeft = R - L;
  const rRight = R + L;
  const coth = 1 + 2 / Math.expm1(2 * mTrue * rLeft);
  const a = mTrue * coth - 1 / rLeft;
  const b = Math.SQRT2 + 1 / rRight;
  let kineticWall = 0;
  let potentialWall = 0;
  for (let j = 0; j <= n; j++) {
    const r = rLeft + j * h;
    const derivative = j === 0 ? a * (phi[0] - trueVacuum)
      : j === n ? -b * phi[n]
        : (phi[j + 1] - phi[j - 1]) / (2 * h);
    const weight = j === 0 || j === n ? 1 : j % 2 ? 4 : 2;
    kineticWall += weight * r ** 2 * derivative ** 2;
    potentialWall += weight * r ** 2 * W(phi[j], eps);
  }
  const kinetic = 2 * Math.PI * h * kineticWall / 3;
  const potential = 4 * Math.PI * (
    W(trueVacuum, eps) * rLeft ** 3 / 3 + h * potentialWall / 3
  );
  return { kinetic, potential, derrick_residual: Math.abs(kinetic + 3 * potential)
    / (Math.abs(kinetic) + 3 * Math.abs(potential)) };
}

function solve(eps, L, h) {
  const { barrier, trueVacuum } = vacua(eps);
  const phase = trueVacuum / Math.SQRT2;
  const gap = -W(trueVacuum, eps);
  const sigma = 4 / (3 * Math.SQRT2);
  const Rtw = 2 * sigma / gap;
  let low = 0.8 * Rtw;
  let high = 1.2 * Rtw;
  let leftResult = matching(eps, low, L, h, trueVacuum, phase);
  let rightResult = matching(eps, high, L, h, trueVacuum, phase);
  if (leftResult.mismatch * rightResult.mismatch > 0)
    throw new Error(`radius not bracketed: Rtw=${Rtw}, g(low)=${leftResult.mismatch}, g(high)=${rightResult.mismatch}`);
  let result;
  let rootIterations = 0;
  for (; rootIterations < 60; rootIterations++) {
    const mid = (low + high) / 2;
    result = matching(eps, mid, L, h, trueVacuum, phase);
    if (Math.abs(result.mismatch) < 2e-10 && high - low < 1e-4) {
      low = high = mid;
      break;
    }
    if (leftResult.mismatch * result.mismatch <= 0) {
      high = mid;
      rightResult = result;
    } else {
      low = mid;
      leftResult = result;
    }
    if (high - low < 2e-7) break;
  }
  const R = (low + high) / 2;
  result = matching(eps, R, L, h, trueVacuum, phase);
  const m = Math.round(L / h);
  const s = Float64Array.from({ length: 2 * m + 1 }, (_, j) => -L + j * h);
  const phi = Float64Array.from({ length: 2 * m + 1 }, (_, j) =>
    j <= m ? result.left.y[j] : result.right.y[j - m]);
  let increasingCount = 0;
  for (let j = 1; j < phi.length; j++)
    if (phi[j] > phi[j - 1] + 1e-11) increasingCount++;
  const halfCenter = trueVacuum / Math.SQRT2;
  const width = crossing(s, phi, trueVacuum * 0.5)
    - crossing(s, phi, trueVacuum * Math.sqrt(0.75));
  const rPhase = R;
  const rHalfCenter = crossing(s, phi, halfCenter) + R;
  const ell = width;
  const virial = virialFromWall(eps, R, L, h, phi, trueVacuum);
  const mTrue = Math.sqrt(Wpp(trueVacuum, eps));
  const centerGapLog10 = Math.log10(Math.abs(trueVacuum - phi[0]))
    + Math.log10(2 * (R - L) * mTrue)
    - mTrue * (R - L) / Math.LN10;
  const summary = { epsilon: eps, phi_tv: trueVacuum, phi_barrier: barrier,
    R_tw: Rtw, R_phase: rPhase, R_half_center: rHalfCenter,
    ell, Lambda_3: R / ell, L, h, phase,
    derivative_mismatch: result.mismatch,
    left_derivative: result.leftDerivative,
    right_derivative: result.rightDerivative,
    left_residual: result.left.norm, right_residual: result.right.norm,
    left_iterations: result.left.iterations,
    right_iterations: result.right.iterations,
    root_iterations: rootIterations, monotonicity_violations: increasingCount,
    phi_inner_edge: phi[0], phi_outer_edge: phi.at(-1),
    inner_edge_deviation: trueVacuum - phi[0],
    estimated_center_gap_log10: centerGapLog10,
    ...virial,
  };
  return { summary, s, phi };
}

function label(value) { return String(value).replace('.', 'p'); }

const summaries = [];
for (const eps of EPSILONS) {
  for (const { L, h } of CONFIGS) {
    const result = solve(eps, L, h);
    summaries.push(result.summary);
    console.log(JSON.stringify(result.summary));
    const saveMain = L === 15
      && (h === 0.01 || h === 0.005 || h === 0.0025
        || h === 0.00125 || h === 0.000625);
    const saveWindowCheck = eps <= 0.0015 && h === 0.0025
      && (L === 12 || L === 20);
    if (saveMain || saveWindowCheck || process.argv[3]) {
      const stem = `quartic_wall_bvp_profile_eps${label(eps)}_L${L}_h${label(h)}`;
      const filename = `${stem}.csv`;
      fs.writeFileSync(path.join(directory, filename),
        ['r,phi', ...Array.from(result.phi, (value, j) =>
          `${(result.summary.R_phase + result.s[j]).toPrecision(17)},${value.toPrecision(17)}`)]
          .join('\n') + '\n');
      fs.writeFileSync(path.join(directory, `${stem}.json`),
        JSON.stringify(result.summary, null, 2) + '\n');
    }
  }
}
const summaryName = process.argv[3]
  ? `quartic_wall_bvp_summary_eps${EPSILONS.map(label).join('_')}_${CONFIGS.map(({ L, h }) => `L${label(L)}_h${label(h)}`).join('_')}.json`
  : process.argv[2]
    ? `quartic_wall_bvp_summary_eps${EPSILONS.map(label).join('_')}.json`
    : 'quartic_wall_bvp_summary.json';
fs.writeFileSync(path.join(directory, summaryName),
  JSON.stringify(summaries, null, 2) + '\n');
