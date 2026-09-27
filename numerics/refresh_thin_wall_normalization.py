"""Refresh derived ratios, fit, and figures without resolving the 49 BVPs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import reconstruct_thin_wall as model


ROOT = Path(__file__).resolve().parents[1]
OBSERVABLES = ROOT / "numerics" / "thin_wall_observables.csv"


def main() -> None:
    with OBSERVABLES.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if fieldnames is None or len(rows) != len(model.EPSILONS):
        raise RuntimeError("Unexpected source table; do not replace it")

    for row in rows:
        radius = float(row["R"])
        ell = float(row["ell"])
        spatial = float(row["Sx_3"])
        spectral = float(row["Sk_3"])
        product = float(row["U_3"])
        lambda_3 = float(row["Lambda_3"])
        epsilon = float(row["epsilon"])
        ratios = {
            "ratio_radius": radius / float(row["R_tw"]),
            "ratio_width": ell / model.ELL_WALL,
            "ratio_sx": spatial / (4.0 * model.np.pi * model.I_WALL * radius**2),
            "ratio_sk": radius**3 * spectral / model.J_BALL,
            "ratio_u_radius": radius * product / model.A_R,
            "ratio_u_lambda": lambda_3 * product / model.C_LAMBDA_3,
            "ratio_u_epsilon": product / (model.A_EPS * epsilon),
        }
        for key, value in ratios.items():
            row[key] = repr(float(value))

    temporary = OBSERVABLES.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(OBSERVABLES)

    numeric_rows = [
        {key: int(value) if key == "nodes" else float(value)
         for key, value in row.items()}
        for row in rows
    ]
    model.plot_ratios(numeric_rows)
    fit = model.plot_scale_law(numeric_rows)
    summary = {
        "constants": {
            "I_wall": model.I_WALL,
            "ell_wall": model.ELL_WALL,
            "J_ball": model.J_BALL,
            "C_Lambda_3": model.C_LAMBDA_3,
            "A_R": model.A_R,
            "A_epsilon": model.A_EPS,
        },
        "fit": fit,
        "smallest_epsilon": numeric_rows[-1],
    }
    summary_path = ROOT / "numerics" / "thin_wall_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
