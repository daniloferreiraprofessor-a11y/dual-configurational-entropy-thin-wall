# Independent large-scale quartic-bubble observables

Date: 2026-11-12. This is a numerical validation, not an
additional asymptotic proof. The profile comes from the independent
wall-centered solver `extended_scale/quartic_wall_bvp_solver.mjs`, which fixes the
phase at $\phi(R)=\phi_{\mathrm{tv}}/\sqrt{2}$ and determines $R$ by matching
the two wall halves. It does not use the archived SciPy profiles.

The independent observable routine is
`extended_scale/quartic_wall_spectrum_direct.mjs`. For $A=\phi_{\mathrm{tv}}$ it
uses the exact sharp-ball transform plus a direct Simpson integral of the
localized wall difference $w(r)=\phi(r)-A\mathbf 1_{r<R}$:

```text
F(k) = A [sin(kR)-kR cos(kR)]/k^3
     + integral_{R-L}^{R+L} r^2 w(r) sinc(kr) dr,
F(0) = A R^3/3 + integral_{R-L}^{R+L} r^2 w(r) dr.
```

The integral is split at $R$ so that the jump in $w$ is integrated with
separate left and right endpoint values. The code uses a series for the
sharp-ball term near $kR=0$, then integrates
$S_k=4\pi\int k^2[-q(k)\ln q(k)]\,dk$, with $q=[F(k)/F(0)]^2$.
$S_x$ is integrated over the wall with $M=\phi_{\mathrm{tv}}$; the omitted
core/tail and the difference between $M$ and the actual central value are
exponentially small. Varying $L$ tests the truncation directly.

## New independent points

These values use the $L=15$, $h_{\mathrm{BVP}}=0.000625$ profiles. For
$\epsilon=0.0015$, the spectral quadrature uses effective wall step
$0.005$, $\Delta k=0.00003125$, and $k_{\max}=8$; for $\epsilon=0.001$ it
uses the same wall step, $\Delta k=0.000015625$, and $k_{\max}=4$.

| $\epsilon$ | $R$ | $\ell$ | $\Lambda_3$ | $S_x^{(3)}$ | $S_k^{(3)}$ | $\mathcal U_3$ | $H$ | Derrick residual |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0015 | 157.563909 | 1.313829 | 119.927227 | 283800.314 | $3.6939133\times10^{-5}$ | 10.4833375 | 238.675 | $2.72\times10^{-8}$ |
| 0.0010 | 236.134088 | 1.315743 | 179.468307 | 637948.255 | $1.1291216\times10^{-5}$ | 7.2032113 | 203.839 | $5.44\times10^{-8}$ |

Here the tabulated numerical diagnostic
$H_{\mathrm{num}}=\Lambda_3^2\mathcal U_3-C_{\mathrm{num}}\Lambda_3
+D_4\ln\Lambda_3$ uses the analytic $D_4=3307.053668612$ and the
same quadrature-based $C_{\mathrm{num}}=1387.248201576098$ used for
the earlier data. No fit is involved. The corrected asymptote gives $\mathcal U_{\mathrm{asym}}
=10.4667427$ and $7.1968827$, respectively. The residuals in $\mathcal U_3$
are about 0.158% and 0.088% of the computed values.

## Overlap with the archived sequence

At $\epsilon=0.003$, the independent wall solve gives
$\Lambda_3=60.377473741$, compared with the archived $60.377476098$.
The direct wall-correction quadrature gives
$\mathcal U_3=19.3437108355$ and $H=318.73997$.
The earlier padded-sine calculation with padding 64 gives
$\mathcal U_3=19.3437113701$, a difference of $5.35\times10^{-7}$;
the manuscript dataset used padding 32 and gave $19.3435597266$.
Thus the independent spectrum reproduces the finer padding result at the
overlap point. This comparison also makes clear why a few decimal places
of $H$ should not be inferred from the original 49-point table alone.

## Numerical sensitivity

For the $\epsilon=0.0015$ fine profile, decreasing $\Delta k$ from
$0.000125$ to $0.0000625$ to $0.00003125$ at $k_{\max}=8$ gives
$H=238.66610$, $238.67565$, and $238.67534$. At fixed
$\Delta k=0.00003125$, using $k_{\max}=4$ rather than 8 changes $H$
by 0.00124; $k_{\max}=2$ misses about 6.39.

For the $\epsilon=0.001$ fine profile, $\Delta k=0.0000625$,
$0.00003125$, and $0.000015625$ at $k_{\max}=4$ give
$H=203.84400$, $203.83867$, and $203.83881$. At the first step,
raising $k_{\max}$ from 4 to 8 changes $H$ by 0.00127; $k_{\max}=2$
misses about 6.70. Sampling the same fine BVP profile at effective wall
steps $0.005$, $0.0025$, and $0.00125$ changes $H$ by less than
$2\times10^{-8}$ at fixed spectral step.

The BVP discretization matters more: at $\epsilon=0.0015$, decreasing
$h_{\mathrm{BVP}}$ from $0.00125$ to $0.000625$ changes $H$ by about
0.011 at fixed spectral settings; at $\epsilon=0.001$ the corresponding
change is about 0.029. For the latter at $h_{\mathrm{BVP}}=0.0025$,
changing the wall window $L=12,15,20$ gives $H=203.75853$,
$203.76296$, $203.76301$, respectively. Thus $L=15$ is adequate at
the precision discussed here.

The computed normalized spectral power never exceeded one. As a separate
cross-check, the field equation gives, for $k>0$,
$F(k)=-k^{-3}\int r W_\epsilon'(\phi(r))\sin(kr)\,dr$.
At $k=0.5,1,2,3,5$, this formula agrees with the ball-plus-wall
transform to relative error below $6\times10^{-6}$ at both new points;
the largest absolute difference normalized by $F(0)$ is below
$6\times10^{-12}$.

An additional BVP refinement to $h_{\mathrm{BVP}}=0.0003125$ was
discarded because its equation and virial residuals worsened; it does
not justify additional displayed digits. The displayed raw $H$ values
are not error bars. BVP mesh effects,
spectral resolution, and the four-decimal input for $J_{\mathrm{ball}}$
support reporting roughly $H\simeq239$ and $H\simeq204$, rather than
interpreting their final digits physically. The proven statement is
$H_{\mathrm{exact}}=O(1)$ with the exact $J_{\mathrm{ball}}$; it does not
imply that $H$ approaches a constant. At finite $\Lambda_3$ the
four-decimal rounding of $J_{\mathrm{ball}}$ affects $H_{\mathrm{num}}$
by roughly $0.08$ or less if that numerical precision is accurate;
a fixed rounded coefficient would eventually introduce a linear drift.
These new
finite-scale points strengthen the numerical check over a threefold
larger scale range without establishing an additional asymptotic law.
