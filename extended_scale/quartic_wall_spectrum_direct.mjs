// Independent radial Fourier quadrature for a wall-centered quartic bubble.
// The sharp-ball contribution is analytic; only the localized wall is sampled.
// Usage: node extended_scale/quartic_wall_spectrum_direct.mjs PROFILE.csv SUMMARY.json [dk] [kmax] [stride]
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';

const profilePath = resolve(process.argv[2] ??
  'extended_scale/quartic_wall_bvp_profile_eps0p0015_L15_h0p000625.csv');
const summaryPath = resolve(process.argv[3] ??
  'extended_scale/quartic_wall_bvp_profile_eps0p0015_L15_h0p000625.json');
const targetDk = Number(process.argv[4] ?? 0.001);
const kMax = Number(process.argv[5] ?? 12);
const stride = Number(process.argv[6] ?? 1);
if (!(targetDk > 0 && kMax > 0)) throw new Error('dk and kmax must be positive');
if (!Number.isInteger(stride) || stride < 1) throw new Error('stride must be a positive integer');

const constants = JSON.parse(readFileSync('numerics/thin_wall_summary.json', 'utf8')).constants;
const summaryFile = JSON.parse(readFileSync(summaryPath, 'utf8'));
const lines = readFileSync(profilePath, 'utf8').trim().split(/\r?\n/);
if (lines.shift() !== 'r,phi') throw new Error('Expected r,phi CSV header');
if ((lines.length - 1) % stride !== 0) {
  throw new Error('The profile interval count must be divisible by stride');
}
const profile = lines.filter((_, index) => index % stride === 0)
  .map((line) => line.split(',').map(Number));
const n = profile.length - 1;
const h = profile[1][0] - profile[0][0];
const summary = Array.isArray(summaryFile) ? summaryFile.find((item) =>
  Math.abs(item.R_phase - item.L - profile[0][0]) < 1e-6 &&
  Math.abs(item.h - h) < 1e-8) : summaryFile;
if (!summary) throw new Error('No matching summary row for profile grid');
const R = summary.R_phase;
const radius = summary.R_half_center ?? R;
const L = summary.L;
const A = summary.phi_tv;
const eps = summary.epsilon;
const ell = summary.ell;
if (![R, radius, L, A, eps, ell].every(Number.isFinite)) {
  throw new Error('Summary must contain epsilon, phi_tv, R_phase, ell, L');
}
if (n % 2 !== 0) throw new Error('Wall CSV needs an even number of intervals');
const middle = n / 2;
if (Math.abs(profile[middle][0] - R) > 1e-6 ||
    Math.abs(n * h - 2 * L) > 1e-6) {
  throw new Error('CSV must be a uniform R-L to R+L grid including R');
}
for (let j = 1; j <= n; j++) {
  if (Math.abs(profile[j][0] - (profile[0][0] + j * h)) > 1e-6) {
    throw new Error(`Nonuniform grid at row ${j}`);
  }
}
if (middle % 2 !== 0) throw new Error('Each half-wall needs an even number of intervals');

const pi = Math.PI;
function entropy(q) { return q > 0 ? -q * Math.log(q) : 0; }
function simpsonWeight(j, begin, end) {
  return (j === begin || j === end ? 1 : ((j - begin) % 2 ? 4 : 2)) / 3;
}
function ballKernel(x) {
  if (Math.abs(x) < 0.3) {
    const y = x * x;
    return 1 / 3 - y / 30 + y * y / 840 - y ** 3 / 45360 + y ** 4 / 3991680;
  }
  return (Math.sin(x) - x * Math.cos(x)) / x ** 3;
}
function potentialPrime(phi) {
  return phi * (phi * phi - 3 * (1 + eps) * phi + 2);
}

const r = new Float64Array(n + 1);
const phi = new Float64Array(n + 1);
const left = new Float64Array(middle + 1);
const right = new Float64Array(middle + 1);
const ode = new Float64Array(n + 1);
let sx = 0;
let maxPhiOverA = 0;
for (let j = 0; j <= n; j++) {
  r[j] = profile[j][0];
  phi[j] = profile[j][1];
  maxPhiOverA = Math.max(maxPhiOverA, phi[j] / A);
  const qx = Math.min(1, (phi[j] / A) ** 2);
  sx += simpsonWeight(j, 0, n) * r[j] ** 2 * entropy(qx);
  ode[j] = h * simpsonWeight(j, 0, n) * r[j] * potentialPrime(phi[j]);
}
sx *= 4 * pi * h;
for (let j = 0; j <= middle; j++) {
  left[j] = h * simpsonWeight(j, 0, middle) * r[j] ** 2 * (phi[j] - A);
  const p = middle + j;
  right[j] = h * simpsonWeight(p, middle, n) * r[p] ** 2 * phi[p];
}
const wallZero = left.reduce((sum, v) => sum + v, 0) +
  right.reduce((sum, v) => sum + v, 0);
const transformZero = A * R ** 3 / 3 + wallZero;

function sinc(x) {
  return Math.abs(x) < 1e-4 ? 1 - x * x / 6 + x ** 4 / 120 : Math.sin(x) / x;
}
function transform(k) {
  if (k === 0) return transformZero;
  let correction = 0;
  for (let j = 0; j <= middle; j++) {
    correction += left[j] * sinc(k * r[j]);
    correction += right[j] * sinc(k * r[middle + j]);
  }
  return A * R ** 3 * ballKernel(k * R) + correction;
}
function odeTransform(k) {
  let sum = 0;
  for (let j = 0; j <= n; j++) sum += ode[j] * Math.sin(k * r[j]);
  return -sum / k ** 3;
}

const nK = 2 * Math.ceil(kMax / targetDk / 2);
const dk = kMax / nK;
let skSum = 0;
let maxQ = 0;
let overshootCount = 0;
for (let j = 1; j <= nK; j++) {
  const k = j * dk;
  const ratio = transform(k) / transformZero;
  const rawQ = ratio * ratio;
  maxQ = Math.max(maxQ, rawQ);
  if (rawQ > 1 + 1e-10) overshootCount++;
  const q = Math.min(1, rawQ);
  skSum += simpsonWeight(j, 0, nK) * k * k * entropy(q);
}
const sk = 4 * pi * dk * skSum;
const U = sx * sk;
const lambda = radius / ell;
const C = constants.C_Lambda_3;
const D = constants.D_Lambda_3;
const H = lambda ** 2 * U - C * lambda + D * Math.log(lambda);

const odeChecks = [0.5, 1, 2, 3, 5].map((k) => {
  const f = transform(k);
  const fOde = odeTransform(k);
  return { k, f, fOde, relative_difference: (fOde - f) / f,
    normalized_absolute_difference: Math.abs(fOde - f) / transformZero };
});

console.log(JSON.stringify({
  profilePath,
  summaryPath,
  epsilon: eps,
  R_phase: R,
  R_half_center: radius,
  ell,
  Lambda_3: lambda,
  phi_tv: A,
  left_endpoint_phi: phi[0],
  right_endpoint_phi: phi[n],
  max_phi_over_phi_tv: maxPhiOverA,
  wall_h: h,
  profile_stride: stride,
  spectral_dk: dk,
  kMax,
  transformZero,
  Sx_3: sx,
  Sk_3: sk,
  U_3: U,
  H,
  max_raw_q: maxQ,
  q_overshoot_count: overshootCount,
  odeChecks,
}, null, 2));
