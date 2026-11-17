# Independent Fortran check of the O(3) bubble action

`bubble_action_check.f90` solves the radial boundary-value equation

`phi'' + 2 phi'/r = W'(phi),  phi'(0)=0,  phi(r_max)=0`

on a uniform grid. It uses a tridiagonal Newton iteration, rather than
reading a saved profile. The quartic and sextic potentials are the ones
in the paper. The initial wall shape and its `O(1)` radius offset only
start the iteration; the reported profile must satisfy the discretized
equation and the monotonicity check.

The action is calculated independently as `B_3=T+V`, with composite
Simpson quadrature for both radial integrals. A second evaluation,
`2T/3`, checks the Derrick identity. The program computes the **exact**
gap `Delta W=-W(phi_tv)` before forming
`B_TW=16*pi*sigma_0^3/(3*(Delta W)^2)`. Its input parameter is the
potential's `epsilon` or sextic `delta`, not `Delta W`.

## Build and run

In PowerShell, from this folder:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
gfortran -O2 -std=f2008 -Wall -Wextra -fcheck=all bubble_action_check.f90 -o bubble_action_check.exe
.\bubble_action_check.exe quartic 0.003 0.0125
.\bubble_action_check.exe quartic 0.003 0.00625
.\bubble_action_check.exe quartic 0.001 0.0125
.\bubble_action_check.exe quartic 0.001 0.00625
.\bubble_action_check.exe sextic 0.005 0.0125
.\bubble_action_check.exe sextic 0.005 0.00625
```

The same source can be compiled with another Fortran 2008 compiler.
An optional fourth argument names a CSV file for the computed
`r,phi` profile.

## Numerical audit

Read `observed B/B_TW - 1`, `first-order prediction`, `difference`,
`maximum equation residual`, and `Derrick residual` from each run.
Halve `h` at fixed potential parameter. The second-order grid error
can be estimated with Richardson extrapolation:

`Q_0 = [4 Q(h/2) - Q(h)]/3`, where `Q=B_3/B_TW-1`.

Results from GNU Fortran 16.1.0 on 2026-11-16:

| Potential parameter | `Q(0.0125)` | `Q(0.00625)` | extrapolated `Q_0` | first-order term |
|:--|--:|--:|--:|--:|
| Quartic `epsilon=0.003` | 0.082856319 | 0.082882139 | 0.082890746 | 0.082093558 |
| Quartic `epsilon=0.001` | 0.027176519 | 0.027200736 | 0.027208809 | 0.027121501 |
| Sextic `delta=0.005` | 0.070796111 | 0.070818209 | 0.070825575 | 0.070414948 |

The extrapolated quartic values agree with the independently saved
action values 0.082890746 and 0.027208728 at the displayed accuracy
of the physical comparison. The sextic value quoted in the current
manuscript, 0.070796, is the **finite-grid** value at `h=0.0125`;
its extrapolated value is about 0.070826. The difference is a grid
effect and should not be attributed to the asymptotic remainder.

The code checks the action calculation. It does not independently
prove the thin-wall coefficient or the uniform remainder in Appendix D.

## Three-grid check at the two smallest gaps

`mesh_refinement_check.f90` runs `bubble_action_check.exe` at
`h=0.0125`, `0.00625`, and `0.003125` for quartic `epsilon=0.00025`
and sextic `delta=0.00125`. Run it from the directory containing
`bubble_action_check.exe`. It writes six `mesh_*.log` files and prints
two Richardson extrapolations from consecutive grid pairs and their
difference.

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
gfortran -O2 -std=f2008 -Wall -Wextra -fcheck=all mesh_refinement_check.f90 -o mesh_refinement_check.exe
.\mesh_refinement_check.exe
```

Observed on 2026-11-17:

| Model | `Q_0` from `h=0.0125,0.00625` | `Q_0` from `h=0.00625,0.003125` | `[(Q_0-Q_1)/Delta W^2]` from finer pair |
|:--|--:|--:|--:|
| Quartic `epsilon=0.00025` | 0.006763020744786 | 0.006763020770238 | 1.353706794 |
| Sextic `delta=0.00125` | 0.017523588764930 | 0.017523588783288 | 16.219065298 |

The two extrapolations differ by `2.55e-11` and `1.84e-11`,
respectively, in `Q`. This tests grid convergence at fixed potential
parameters; it is not an estimate of all systematic errors in the
continuum action or a proof of the `O((Delta W)^2)` law. The compact
results for all eight tested parameters are in `action_grid_summary.csv`.
Full radial profiles are not bundled: to regenerate one, pass an output
CSV name as the optional fourth argument, for example
`bubble_action_check.exe quartic 0.00025 0.003125 profile_quartic_h003125.csv`.
The CSV then has one header row and one row for each radial grid node.
