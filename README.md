# Dual configurational entropy and thin-wall scaling: reproducibility package

This repository contains the numerical code and data accompanying
the manuscript *Dual configurational entropy and universal thin-wall scaling*.
The author reports having reviewed the manuscript's scientific claims and
independently cross-checked the numerical results represented in the figures
using Wolfram Mathematica 15.0 and Fortran. Those separate implementations
are not included here; the reproducibility route provided here uses Python.
The repository is available at
https://github.com/daniloferreiraprofessor-a11y/dual-configurational-entropy-thin-wall.
No archival DOI has been assigned yet.

## Licensing

Copyright (c) 2026 Danilo Cardoso Ferreira. The Python source code and its
documentation are licensed under the MIT license in `LICENSE`. The numerical
data (`numerics/*.csv`, `numerics/*.json`, `numerics/*.npz`) and manuscript
figure files (`figures/*.pdf`, `figures/*.png`) are licensed under Creative
Commons Attribution 4.0 International (CC BY 4.0) as specified in
`LICENSE-DATA.md`. The manuscript itself is not included or licensed by this
package.

## Contents

- `numerics/reconstruct_thin_wall.py`: solves the radial boundary-value
  problem for 49 values of `epsilon`, evaluates the observables, fits the
  finite-scale law, and generates all three manuscript figures.
- `numerics/evaluate_ball_spectral_constant.py`: independently evaluates the
  sharp-ball spectral constant by quadrature and tail extrapolation.
- `numerics/check_quadrature_convergence.py`: checks entropy quadrature using
  the saved continuous profiles.
- `numerics/verify_direct_bvp_convergence.py`: resolves two profiles and
  checks sensitivity to numerical resolution and radial-domain size.
- `numerics/refresh_thin_wall_normalization.py`: recalculates derived
  normalization columns and related figures from the saved observable table.
- `numerics/thin_wall_profiles.npz`: radial meshes and fields for the 49
  computed profiles.
- `numerics/thin_wall_observables.csv`: observables for the 49 profiles.
- `numerics/thin_wall_summary.json`: constants, finite-scale fit, and the
  smallest-`epsilon` profile's observables.
- `numerics/quadrature_convergence.csv` and
  `numerics/direct_bvp_convergence.csv`: saved convergence checks.
- `figures/*.pdf` and `figures/*.png`: the three manuscript figures.

## Requirements

Python 3.12 and the packages in `requirements.txt`. In a fresh environment,
install them with `python -m pip install -r requirements.txt` from this
directory. The original calculations used NumPy 2.5.3, SciPy 1.18.1, and
Matplotlib 3.11.2; these versions are pinned for reproducibility. The
plotting scripts use Matplotlib's built-in math-text renderer and do not
require a LaTeX installation. They select the non-interactive `Agg` backend,
so no desktop graphics toolkit is required.

## Reproduction

Run the following commands from this directory. They write output only inside
this package, replacing the corresponding saved CSV, JSON, NPZ, PDF, and PNG
files; retain a separate copy of the released data before rerunning them.

```text
python numerics/evaluate_ball_spectral_constant.py
python numerics/reconstruct_thin_wall.py
python numerics/check_quadrature_convergence.py
python numerics/verify_direct_bvp_convergence.py
```

The full reconstruction solves all 49 boundary-value problems and can take
substantially longer than the two convergence checks. The diagnostic scripts
read `numerics/thin_wall_profiles.npz` and save their output tables under
`numerics/`.

This package records the Python calculations. The author's
cross-software checks are reported separately in `VALIDATION.md`; they cannot
be rerun from this package. No public DOI has been assigned yet.

The local reconstruction checks, including a fresh-environment test from
`requirements.txt`, are documented in `VALIDATION.md`.

The plotting code was developed with assistance from OpenAI ChatGPT and
Codex. The author directed the analysis and reviewed the results.
