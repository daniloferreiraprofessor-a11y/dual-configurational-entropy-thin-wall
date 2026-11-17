"""Check a proposed asymptotic activation-barrier bound against saved data.

This script is a diagnostic, not an independent bubble or spectrum solver.
It reads the wall-centered solutions and spectral results already saved in
this directory. In particular, finite-radius values may lie below the
asymptotic bound without contradicting it.

Run from any directory with ``python check_activation_bound.py``.
Only the Python standard library is required.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


DIRECTORY = Path(__file__).resolve().parent
CONSTANTS = json.loads(
    (DIRECTORY.parent / "numerics" / "thin_wall_summary.json").read_text(encoding="utf-8")
)["constants"]
J_BALL = CONSTANTS["J_ball"]
ASYMPTOTIC_BOUND = math.pi**4 * J_BALL**1.5 / (3 * math.sqrt(2))
SIGMA_I_QUARTIC = (math.pi**2 - 6) / 18
SHARP_SIGMA_I = math.pi / 16


def potential(phi: float, epsilon: float) -> float:
    return phi**2 * (phi**2 - 4 * (1 + epsilon) * phi + 4) / 4


def check_row(row: dict[str, str]) -> tuple[float, float, float, float, float, float]:
    epsilon = float(row["epsilon"])
    tag = "0p0015" if epsilon == 0.0015 else "0p001" if epsilon == 0.001 else None
    if tag is None:
        raise ValueError(f"No saved BVP summary for epsilon={epsilon}")
    path = DIRECTORY / f"quartic_wall_bvp_profile_eps{tag}_L15_h0p000625.json"
    bvp = json.loads(path.read_text(encoding="utf-8"))
    if not math.isclose(epsilon, bvp["epsilon"], rel_tol=0, abs_tol=1e-12):
        raise ValueError(f"Epsilon mismatch in {path.name}")
    if not math.isclose(float(row["R"]), bvp["R_phase"], rel_tol=1e-10):
        raise ValueError(f"Radius mismatch in {path.name}")

    sx, sk, dual = (float(row[key]) for key in ("Sx_3", "Sk_3", "U_3"))
    if not math.isclose(sx * sk, dual, rel_tol=1e-10):
        raise ValueError(f"Entropy product mismatch for epsilon={epsilon}")

    amplitude = bvp["phi_tv"]
    delta_w = -potential(amplitude, epsilon)
    barrier = bvp["kinetic"] + bvp["potential"]
    if delta_w <= 0 or barrier <= 0:
        raise ValueError(f"Nonpositive physical quantity for epsilon={epsilon}")
    if not math.isclose(barrier, 2 * bvp["kinetic"] / 3, rel_tol=1e-6):
        raise ValueError(f"Derrick identity failed for epsilon={epsilon}")

    # Dimensionless combination for a canonical O(3) critical bubble.
    diagnostic = barrier * math.sqrt(delta_w) * dual**1.5 / amplitude**3
    return (
        epsilon,
        float(row["Lambda_3"]),
        barrier,
        delta_w,
        diagnostic,
        diagnostic / ASYMPTOTIC_BOUND,
    )


def main() -> None:
    path = DIRECTORY / "quartic_wall_spectrum_points.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    original = json.loads(
        (DIRECTORY.parent / "numerics" / "thin_wall_summary.json").read_text(
            encoding="utf-8"
        )
    )["smallest_epsilon"]
    original_epsilon = original["epsilon"]
    original_amplitude = original["phi_tv"]
    original_barrier = original["T"] + original["V"]
    if not math.isclose(original_barrier, 2 * original["T"] / 3, rel_tol=1e-6):
        raise ValueError("Derrick identity failed for original sequence")
    original_delta_w = -potential(original_amplitude, original_epsilon)
    original_q = (
        original_barrier
        * math.sqrt(original_delta_w)
        * original["U_3"] ** 1.5
        / original_amplitude**3
    )
    checks = [
        (original_epsilon, original["Lambda_3"], original_barrier,
         original_q / ASYMPTOTIC_BOUND, "original")
    ]
    for row in rows:
        epsilon, lam, barrier, _delta_w, _diagnostic, ratio = check_row(row)
        checks.append((epsilon, lam, barrier, ratio, "extended"))

    print(f"Asymptotic universal bound: {ASYMPTOTIC_BOUND:.6f}")
    print(
        "Quartic asymptotic ratio to bound: "
        f"{(SIGMA_I_QUARTIC / SHARP_SIGMA_I)**1.5:.6f}"
    )
    print("epsilon    Lambda    barrier B3       Q/bound   log prediction   residual    method")
    limit_ratio = (SIGMA_I_QUARTIC / SHARP_SIGMA_I) ** 1.5
    log_coefficient = 1.5 * CONSTANTS["D_Lambda_3"] / CONSTANTS["C_Lambda_3"]
    for epsilon, lam, barrier, ratio, method in checks:
        # Candidate first-log prediction: the action and bulk factors must
        # independently be shown to have no larger logarithmic correction.
        prediction = limit_ratio * (1 - log_coefficient * math.log(lam) / lam)
        print(
            f"{epsilon:0.4f}    {lam:7.3f}   {barrier:12.6f}    "
            f"{ratio:0.6f}      {prediction:0.6f}      "
            f"{ratio - prediction:+0.6f}    {method}"
        )


if __name__ == "__main__":
    main()
