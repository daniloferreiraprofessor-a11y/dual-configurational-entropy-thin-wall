# Local reconstruction check

Status on 2026-09-27: the copied package was exercised in a fresh Python
virtual environment before any public release. This is a software
reproducibility check, not an independent scientific validation of the
manuscript.

- A new Python 3.12.14 environment initially contained only `pip`. The
  pinned requirements installed successfully: NumPy 2.5.3, SciPy 1.18.1,
  and Matplotlib 3.11.2. Its module path did not include `.numerical-deps`.
- The first clean-environment run exposed a missing graphical toolkit in
  Matplotlib's default backend. Both source and package scripts now select
  the non-interactive `Agg` backend before importing `pyplot`.
- All five Python files parsed successfully.
- The independent sharp-ball quadrature returned `J_ball = 159.7247` at four
  decimal places across its truncation and step-size checks.
- The saved-profile quadrature check reproduced the original 14-row CSV
  byte-for-byte.
- The direct boundary-value re-solve check reproduced the original 16-row
  CSV byte-for-byte.
- A full run of `reconstruct_thin_wall.py` regenerated all 49 profiles.
  The observables CSV, profile NPZ, and summary JSON matched the manuscript
  copies byte-for-byte (SHA-256 hashes matched).
- The three regenerated PNG figures matched the manuscript copies byte-for-
  byte. The regenerated PDFs had different file hashes because their embedded
  creation dates differ; rendering each original/reconstructed PDF pair at
  120 dpi produced identical PNG hashes. The three PDF figures were also
  inspected visually.

All five scripts subsequently ran directly in the clean environment, without
access to the original workspace's private package directory.

Separately, the author reports having inspected the manuscript's scientific
claims and found no mathematical errors or unsupported general statements.
The author also reports that the numerical results represented in the three
figures were independently cross-checked using Wolfram Mathematica 15.0 and
Fortran. The Mathematica notebooks and the earlier cross-check implementations
are not part of this package, so this document does not claim to reproduce
those cross-software checks.

## Additional action-correction validation, 2026-11-17

The independent Fortran source in `fortran/bubble_action_check.f90` solves the
quartic and sextic radial boundary-value problems without reading saved
profiles. It was compiled with GNU Fortran 16.1.0 using `-O2 -std=f2008
-Wall -Wextra -fcheck=all`, with no compiler warnings. The action is calculated
as `T+V` by composite Simpson quadrature and compared with `2T/3`.

At four exact vacuum gaps per model, `fortran/action_grid_summary.csv`
records the analytic first-order term, two-grid Richardson estimate of
`Q=B_3/B_TW-1`, and `rho=(Q-Q_first)/(Delta W)^2`. The resulting `rho`
sequence is `1.34970, 1.35200, 1.35314, 1.35371` for the quartic and
`16.08533, 16.16303, 16.20062, 16.21907` for the sextic model as the
gap decreases.

For the smallest gap of each model, a third spacing `h=0.003125` gives
successive-change ratios `3.999987` and `3.999989`. The Richardson
estimates from the two consecutive grid pairs differ by `2.55e-11` and
`1.84e-11` in `Q`. These checks support second-order grid convergence and
the observed quadratic action remainder. They do not prove a uniform
asymptotic expansion or bound all systematic errors. Build and run commands
are in `fortran/README.md`.
