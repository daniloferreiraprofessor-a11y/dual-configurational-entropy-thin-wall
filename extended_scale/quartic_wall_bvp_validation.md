# Independent wall-centered quartic bubble BVP

Numerical validation, 2026-11-12. Script:
`extended_scale/quartic_wall_bvp_solver.mjs`. It does not read the archived
profiles or the original Python BVP solution.

The tilted quartic potential is
`W_epsilon(phi)=phi^2[phi^2-4(1+epsilon)phi+4]/4`.
At each trial radius `R`, two nonlinear finite-difference BVPs solve
`phi_ss+2 phi_s/(R+s)=W_epsilon'(phi)` on `[-L,0]` and `[0,L]`.
The phase is `phi(0)=phi_tv/sqrt(2)`. The inner Robin condition is the
regular linear mode around the true vacuum at `r=R-L`; the outer Robin
condition is the decaying mode around zero at `r=R+L`. The scalar root
condition matches the left and right derivatives at the phase. This
avoids shooting from a central value whose true-vacuum gap is
exponentially below double precision. The reported `R_phase` differs
from the radius defined using `phi_center/sqrt(2)` by an exponentially
small amount; the estimated center gap is about `10^-95` at
`epsilon=0.0015`.

## Calibration against archived quartic data

| `epsilon` | `L` | `h` | `R_phase` | `ell` | `Lambda_3` | virial residual |
|---:|---:|---:|---:|---:|---:|---:|
| 0.003 | 15 | 0.01 | 78.989217186 | 1.308237615 | 60.378341284 | 3.40e-5 |
| 0.003 | 15 | 0.005 | 78.988420133 | 1.308241044 | 60.377573769 | 5.31e-6 |
| 0.003 | 15 | 0.0025 | 78.988332341 | 1.308241648 | 60.377478783 | 9.18e-7 |
| 0.003 | 15 | 0.00125 | 78.988324486 | 1.308241667 | 60.377471914 | 1.73e-7 |
| 0.003 | 15 | 0.000625 | 78.988323100 | 1.308241604 | 60.377473741 | 6.01e-10 |

The archived Python result is `R=78.988324907`, `ell=1.308241583`,
`Lambda_3=60.377476098`. The finest independent solve differs by
`1.81e-6`, `2.1e-8`, and `2.36e-6`, respectively. The virial residual
is `|T+3V|/(|T|+3|V|)`, with the constant true-vacuum core included
in `V`. These comparisons calibrate the independent geometry and
action integrals before the larger-scale run.

## New larger-scale point

| `epsilon` | `L` | `h` | `R_phase` | `ell` | `Lambda_3` | virial residual |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0015 | 15 | 0.01 | 157.567656098 | 1.313834868 | 119.929574030 | 5.88e-5 |
| 0.0015 | 15 | 0.005 | 157.564351612 | 1.313830645 | 119.927444392 | 8.43e-6 |
| 0.0015 | 15 | 0.0025 | 157.563964601 | 1.313829653 | 119.927240389 | 1.36e-6 |
| 0.0015 | 15 | 0.00125 | 157.563919945 | 1.313829386 | 119.927230716 | 2.51e-7 |
| 0.0015 | 15 | 0.000625 | 157.563908782 | 1.313829330 | 119.927227356 | 2.72e-8 |

Changing `L` among 10, 12, 15, and 20 at fixed `h=0.01` changes
`Lambda_3` by less than `5e-11` in the saved run. The interior Robin
condition uses the regular solution proportional to `sinh(m_tv r)/r`;
the outer one uses `exp(-sqrt(2)r)/r`. The full window and mesh sweep is
in `extended_scale/quartic_wall_bvp_summary.json`.

Profiles at `L=15` are saved as
`extended_scale/quartic_wall_bvp_profile_eps0p003_L15_h0p000625.csv` and
`extended_scale/quartic_wall_bvp_profile_eps0p0015_L15_h0p000625.csv`, with
`r,phi` columns and corresponding per-profile JSON files. They cover
`R-15 <= r <= R+15`. For the omitted core, the regular linearized
extension is
`phi(r)-phi_tv=[phi(R-L)-phi_tv](R-L)/r * sinh(m_tv r)/sinh(m_tv(R-L))`;
the central limit is understood at `r=0`. Beyond `R+L`, use
`phi(r)=phi(R+L)(R+L)/r * exp[-sqrt(2)(r-R-L)]`.

## Second larger-scale point

After the first point passed the geometry and virial checks, the same
independent solver was run at `epsilon=0.001`. With `L=15` and
`h=0.000625`, it gives `R_phase=236.134088016`, `ell=1.315742551`,
`Lambda_3=179.468306939`, and virial residual `5.44e-8`. Changing
`h=0.00125` to `0.000625` shifts `Lambda_3` by `2.40e-7`.
At fixed `h=0.0025`, changing `L` from 12 to 15 or 20 shifts
`Lambda_3` by less than `1e-9`. The true-vacuum gap at the center is
estimated near `10^-143`. The `L=20` inner boundary is already at the
double-precision floor, so `L=12` and `L=15` are more informative
for sensitivity tests. Results are in
`extended_scale/quartic_wall_bvp_summary_eps0p001.json`; the finest profile
is `extended_scale/quartic_wall_bvp_profile_eps0p001_L15_h0p000625.csv`.

The geometry and virial checks do not by themselves validate the two
configurational entropies. Their independent quadratures and an error
budget for `H=Lambda_3^2(U_3-U_asym)` are the next calculation.

## Limit of further double-precision refinement

A targeted attempt at `epsilon=0.001`, `h=0.0003125` did not improve the
validation. With `L=15`, the inner-edge field rose about `1.0e-11` above
the true-vacuum value and the virial residual worsened to `4.66e-7`.
With `L=12`, the inner-edge deviation remained positive, but the BVP
residual reached `1.76e-8` and the virial residual was `4.65e-7`, also
worse than at `h=0.000625`. These runs are retained as diagnostics and
are excluded from the reported observables. Further digits require a
better-conditioned or higher-precision solve, not finer double-precision
spacing alone.
