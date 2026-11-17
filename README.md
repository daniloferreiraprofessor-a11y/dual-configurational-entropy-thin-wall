# Dual configurational entropy: radial thin-wall scaling and critical bubbles — reproducibility package

Author: [Danilo Cardoso Ferreira](https://orcid.org/0000-0002-9329-3978)
(ORCID: 0000-0002-9329-3978).

This repository contains the numerical code and data accompanying
the manuscript *Dual configurational entropy: radial thin-wall scaling and critical bubbles*.
The author reports having reviewed the manuscript's scientific claims and
independently cross-checked the numerical results represented in the figures
using Wolfram Mathematica 15.0 and Fortran. The original 49-profile sequence
is reproduced with Python; two independent larger-scale checks use Node.js,
and the new action-correction checks use an independent Fortran radial solver.
The Mathematica notebooks and earlier cross-check implementations are not
included here.
The repository is available at
https://github.com/daniloferreiraprofessor-a11y/dual-configurational-entropy-thin-wall.
The current v1.1.0 package, including the additional wall-centered and
Fortran checks, has the version-specific DOI
[10.5281/zenodo.24318129](https://doi.org/10.5281/zenodo.24318129).
The archived v1.0.0 code and data for the original sequence have the
version-specific DOI
[10.5281/zenodo.23001558](https://doi.org/10.5281/zenodo.23001558).

## Licensing

Copyright (c) 2026 Danilo Cardoso Ferreira. The Python, Node.js, and Fortran source code and their
documentation are licensed under the MIT license in `LICENSE`. The numerical
data (`numerics/*.csv`, `numerics/*.json`, `numerics/*.npz`,
`extended_scale/*.csv`, `extended_scale/*.json`, `fortran/*.csv`) and manuscript
figure files (`figures/*.pdf`, `figures/*.png`) are licensed under Creative
Commons Attribution 4.0 International (CC BY 4.0) as specified in
`LICENSE-DATA.md`. The manuscript itself is not included or licensed by this
package.

## Contents

- `numerics/reconstruct_thin_wall.py`: solves the radial boundary-value
  problem for 49 values of `epsilon`, evaluates the observables, records
  exploratory finite-scale fits, and generates Matplotlib versions of
  the three manuscript figures.
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
- `numerics/thin_wall_summary.json`: constants, exploratory finite-scale fits, and the
  smallest-`epsilon` profile's observables.
- `numerics/quadrature_convergence.csv` and
  `numerics/direct_bvp_convergence.csv`: saved convergence checks.
- `figures/*.pdf` and `figures/*.png`: the three manuscript figures.
- `figures/dual_entropy_scale_separation_plot.tex`: exact PGFPlots source
  of the scale-separation figure, reading the saved observable table.
- `extended_scale/`: separate wall-centered quartic BVP checks at
  `epsilon=0.0015` and `0.0010`, with Node.js scripts, profiles, a
  two-point observable table, and validation notes.
- `fortran/bubble_action_check.f90`: independent quartic and sextic radial
  boundary-value solver and action integrator. It does not read saved profiles.
- `fortran/mesh_refinement_check.f90`: runs the two smallest-gap cases on
  three meshes and compares consecutive Richardson extrapolations.
- `fortran/action_grid_summary.csv`: the eight action-correction checks,
  including exact gaps, analytic first-order predictions, and numerical
  extrapolations. The corresponding source and commands are documented in
  `fortran/README.md`.

## Requirements

Python 3.12 and the packages in `requirements.txt`. In a fresh environment,
install them with `python -m pip install -r requirements.txt` from this
directory. The original calculations used NumPy 2.5.3, SciPy 1.18.1, and
Matplotlib 3.11.2; these versions are pinned for reproducibility. The
plotting scripts use Matplotlib's built-in math-text renderer and do not
require a LaTeX installation. They select the non-interactive `Agg` backend,
so no desktop graphics toolkit is required. The archived scale-separation
figure was typeset from its PGFPlots source; rebuilding that exact PDF
requires a LaTeX installation with PGFPlots.
The independent checks in `extended_scale/` require Node.js 18 or newer;
they use only its standard library. See that directory's README for the
reproduction commands.
The action checks in `fortran/` require a Fortran 2008 compiler; they have
been tested with GNU Fortran 16.1.0. See that directory's README.

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

The original sequence uses Python, the independent larger-scale checks use
Node.js, and the action-correction check uses Fortran. The Mathematica
cross-checks are reported separately in `VALIDATION.md` and cannot be rerun
from this package. The version-specific DOIs above identify the respective
archived packages; the v1.0.0 DOI does not identify the later additions.

The local reconstruction checks, including a fresh-environment test from
`requirements.txt`, are documented in `VALIDATION.md`.

The plotting code and parts of the numerical validation code were developed
with assistance from OpenAI ChatGPT and Codex. The author directed the
analysis and reviewed the results.
