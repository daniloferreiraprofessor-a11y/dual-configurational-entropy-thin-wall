"""Recheck numerical entropies from fresh continuous BVP interpolants."""

from __future__ import annotations

import csv
from pathlib import Path

import reconstruct_thin_wall as model

import numpy as np
from scipy.fft import dst, next_fast_len
from scipy.integrate import simpson, solve_bvp
from scipy.interpolate import PchipInterpolator


ROOT = Path(__file__).resolve().parents[1]
R0 = model.R0


def solve_again(eps: float, extra_radius: float):
    archive = np.load(ROOT / "numerics" / "thin_wall_profiles.npz")
    prefix = f"eps_{eps:.4f}"
    source_r = archive[f"{prefix}_r"]
    source_phi = archive[f"{prefix}_phi"]
    seed = PchipInterpolator(source_r, source_phi)
    r_tw = model.thin_wall_radius(eps)
    # A common multiple of all tested spacings avoids changing the endpoint.
    r_max = np.ceil((max(25.0, r_tw + 18.0) + extra_radius) / 0.003) * 0.003
    left = np.linspace(R0, max(R0 + 1.0e-4, r_tw - 8.0), 250)
    wall = np.linspace(max(R0 + 2.0e-4, r_tw - 8.0), r_tw + 8.0, 900)
    right = np.linspace(r_tw + 8.0, r_max, 350)
    mesh = np.unique(np.concatenate((left, wall, right)))

    old_edge = float(source_r[-1])
    old_value = max(float(seed(old_edge)), 0.0)
    guess_phi = np.empty_like(mesh)
    guess_derivative = np.empty_like(mesh)
    within = mesh <= old_edge
    guess_phi[within] = seed(mesh[within])
    guess_derivative[within] = seed.derivative()(mesh[within])
    beyond = ~within
    if np.any(beyond):
        distance = mesh[beyond] - old_edge
        tail = old_value * old_edge / mesh[beyond] * np.exp(-np.sqrt(2.0) * distance)
        guess_phi[beyond] = tail
        guess_derivative[beyond] = -(
            np.sqrt(2.0) + 1.0 / mesh[beyond]
        ) * tail
    guess = np.vstack((guess_phi, guess_derivative))

    def ode(r, y):
        return np.vstack((y[1], model.potential_prime(y[0], eps) - 2.0 * y[1] / r))

    def boundary(ya, yb):
        return np.array((
            ya[1] - model.potential_prime(ya[0], eps) * R0 / 3.0,
            yb[1] + (np.sqrt(2.0) + 1.0 / r_max) * yb[0],
        ))

    solution = solve_bvp(ode, boundary, mesh, guess, tol=1.0e-8, max_nodes=300_000)
    if not solution.success:
        raise RuntimeError(f"BVP failed at epsilon={eps}, extra_radius={extra_radius}: {solution.message}")
    phi_bar, phi_tv = model.stationary_points(eps)
    phi_c = float(solution.sol(R0)[0])
    if not (phi_bar < phi_c <= phi_tv * (1.0 + 1.0e-10)):
        raise RuntimeError(f"Unexpected branch at epsilon={eps}: phi_c={phi_c}")
    return solution, float(r_max)


def entropies(solution, eps: float, r_max: float, dr: float, padding: int, k_max: float):
    count = round(r_max / dr)
    r = np.arange(count + 1, dtype=float) * dr
    phi = solution.sol(np.maximum(r, R0))[0]
    at_origin = float(solution.sol(R0)[0])
    phi[0] = at_origin - float(model.potential_prime(np.array(at_origin), eps)) * R0**2 / 6.0
    phi_c = float(phi[0])
    q_x = np.clip((phi / phi_c) ** 2, 0.0, 1.0)
    spatial = np.zeros_like(q_x)
    positive = q_x > 0.0
    spatial[positive] = -4.0 * np.pi * r[positive] ** 2 * q_x[positive] * np.log(q_x[positive])
    s_x = float(simpson(spatial, x=r))

    transform_n = next_fast_len(padding * count)
    samples = np.zeros(transform_n, dtype=float)
    samples[:count] = r[1:] * phi[1:]
    sine_sum = 0.5 * dr * dst(samples, type=1)
    k = np.pi * np.arange(1, transform_n + 1) / ((transform_n + 1) * dr)
    keep = k <= k_max
    k = k[keep]
    transform = np.sqrt(2.0 / np.pi) * sine_sum[keep] / k
    zero_mode = np.sqrt(2.0 / np.pi) * simpson(r**2 * phi, x=r)
    q_k = np.clip((transform / zero_mode) ** 2, 0.0, 1.0)
    full_k = np.concatenate(([0.0], k))
    full_q = np.concatenate(([1.0], q_k))
    spectral = np.zeros_like(full_q)
    positive = full_q > 0.0
    spectral[positive] = -4.0 * np.pi * full_k[positive] ** 2 * full_q[positive] * np.log(full_q[positive])
    s_k = float(simpson(spectral, x=full_k))
    return s_x, s_k, s_x * s_k


def main():
    tests = (
        (0.003, 32, 20.0),
        (0.0015, 16, 20.0),
        (0.0015, 32, 10.0),
        (0.0015, 32, 20.0),
        (0.0015, 32, 30.0),
        (0.0015, 64, 20.0),
        (0.00075, 32, 20.0),
    )
    rows = []
    for eps in (0.006, 0.003):
        for extra_radius in (0.0, 6.0):
            solution, r_max = solve_again(eps, extra_radius)
            selected = tests if extra_radius == 0.0 else ((0.0015, 64, 20.0),)
            for dr, padding, k_max in selected:
                s_x, s_k, product = entropies(solution, eps, r_max, dr, padding, k_max)
                row = dict(
                    epsilon=eps, extra_radius=extra_radius, r_max=r_max,
                    nodes=solution.x.size, dr=dr, padding=padding, k_max=k_max,
                    Sx_3=s_x, Sk_3=s_k, U_3=product,
                )
                rows.append(row)
                print(row, flush=True)
    path = ROOT / "numerics" / "direct_bvp_convergence.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
