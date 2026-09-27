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
Fortran. Those independent implementations and their raw outputs are not part
of this Python package, so this document does not claim to reproduce those
cross-software checks.
