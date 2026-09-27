from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

import numpy as np
from scipy.fft import dst, next_fast_len
from scipy.integrate import simpson
from scipy.interpolate import PchipInterpolator


def entropies(
    source_r: np.ndarray,
    source_phi: np.ndarray,
    dr: float,
    padding: int,
    k_max: float,
) -> tuple[float, float, float]:
    r_max = float(source_r[-1])
    count = int(np.floor(r_max / dr))
    r = np.arange(count + 1, dtype=float) * dr
    phi = PchipInterpolator(source_r, source_phi)(r)
    phi_c = float(phi[0])

    q_x = np.clip((phi / phi_c) ** 2, 0.0, 1.0)
    h_x = np.zeros_like(q_x)
    positive_x = q_x > 0.0
    h_x[positive_x] = (
        -4.0
        * np.pi
        * r[positive_x] ** 2
        * q_x[positive_x]
        * np.log(q_x[positive_x])
    )
    s_x = float(simpson(h_x, x=r))

    physical_n = len(r) - 1
    transform_n = next_fast_len(padding * physical_n)
    samples = np.zeros(transform_n, dtype=float)
    j = np.arange(1, physical_n + 1)
    samples[j - 1] = r[j] * phi[j]
    sine_sum = 0.5 * dr * dst(samples, type=1)
    k = np.pi * np.arange(1, transform_n + 1) / (
        (transform_n + 1) * dr
    )
    keep = k <= k_max
    k = k[keep]
    transform = np.sqrt(2.0 / np.pi) * sine_sum[keep] / k
    transform_zero = np.sqrt(2.0 / np.pi) * simpson(r**2 * phi, x=r)
    q_k = np.clip((transform / transform_zero) ** 2, 0.0, 1.0)
    k_full = np.concatenate(([0.0], k))
    q_full = np.concatenate(([1.0], q_k))
    h_k = np.zeros_like(q_full)
    positive_k = q_full > 0.0
    h_k[positive_k] = (
        -4.0
        * np.pi
        * k_full[positive_k] ** 2
        * q_full[positive_k]
        * np.log(q_full[positive_k])
    )
    s_k = float(simpson(h_k, x=k_full))
    return s_x, s_k, s_x * s_k


def main() -> None:
    archive = np.load(ROOT / "numerics" / "thin_wall_profiles.npz")
    tests = [
        (3.0e-3, 16, 20.0),
        (1.5e-3, 16, 20.0),
        (1.5e-3, 32, 10.0),
        (1.5e-3, 32, 20.0),
        (1.5e-3, 32, 30.0),
        (1.5e-3, 64, 20.0),
        (7.5e-4, 32, 20.0),
    ]
    rows: list[dict[str, float | int]] = []
    for eps in (0.006, 0.003):
        prefix = f"eps_{eps:.4f}"
        source_r = archive[f"{prefix}_r"]
        source_phi = archive[f"{prefix}_phi"]
        for dr, padding, k_max in tests:
            s_x, s_k, product = entropies(
                source_r, source_phi, dr, padding, k_max
            )
            rows.append(
                {
                    "epsilon": eps,
                    "dr": dr,
                    "padding": padding,
                    "k_max": k_max,
                    "Sx_3": s_x,
                    "Sk_3": s_k,
                    "U_3": product,
                }
            )

    path = ROOT / "numerics" / "quadrature_convergence.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"eps={row['epsilon']:.3f}, dr={row['dr']:.5g}, "
            f"pad={row['padding']:d}, kmax={row['k_max']:.0f}: "
            f"U={row['U_3']:.10f}"
        )


if __name__ == "__main__":
    main()
