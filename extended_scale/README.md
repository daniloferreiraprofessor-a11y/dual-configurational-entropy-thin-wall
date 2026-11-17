# Independent large-scale quartic checks

This directory contains two quartic critical-bubble solutions calculated
with a wall-centered boundary-value method, separately from the original
49-profile Python sequence. The overlap run at `epsilon=0.003` calibrates
the method against that sequence. The two additional points are at
`epsilon=0.0015` and `0.0010`.

The `quartic_wall_bvp_solver.mjs` script solves the two half-wall problems
with an unknown radius. The `quartic_wall_spectrum_direct.mjs` script
calculates the position entropy and evaluates the Fourier entropy from
the analytic sharp-ball transform plus direct quadrature of the localized
wall correction. Neither script reads the archived Python profiles.

Run these commands from the package root with Node.js 18 or newer:

```text
node extended_scale/quartic_wall_spectrum_direct.mjs extended_scale/quartic_wall_bvp_profile_eps0p0015_L15_h0p000625.csv extended_scale/quartic_wall_bvp_profile_eps0p0015_L15_h0p000625.json 0.00003125 8 8
node extended_scale/quartic_wall_spectrum_direct.mjs extended_scale/quartic_wall_bvp_profile_eps0p001_L15_h0p000625.csv extended_scale/quartic_wall_bvp_profile_eps0p001_L15_h0p000625.json 0.000015625 4 8
node extended_scale/quartic_wall_spectrum_points.mjs
```

The final command regenerates `quartic_wall_spectrum_points.csv` from the
saved result JSON files. To recompute the profiles, run
`node extended_scale/quartic_wall_bvp_solver.mjs 0.003,0.0015,0.001`;
this may take longer and writes additional convergence profiles in this
directory. The saved BVP and spectral validation notes document window,
mesh, spectral-step, and cutoff checks. The most extreme attempted BVP
refinement in double precision was rejected after its residuals worsened.

The optional standard-library Python script
`extended_scale/check_activation_bound.py` combines the saved BVP actions
and entropy products, including the original-sequence overlap point, to
inspect a proposed asymptotic activation-barrier bound and its first
logarithmic correction. Run `python extended_scale/check_activation_bound.py`
from the package root. It is a post-processing diagnostic, not an
independent solution of the boundary-value or spectral problems. Its
finite-radius numbers may fall below the asymptotic limit.

These checks were added after the v1.0.0 archive and are therefore not
identified by its version-specific DOI.
