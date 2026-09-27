from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter
from scipy.fft import dst, next_fast_len
from scipy.integrate import quad, simpson, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_bvp
from scipy.optimize import brentq


EPSILONS = np.array(
    [0.050, 0.040, 0.030, 0.025, 0.020, 0.018, 0.017, 0.016,
     0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.009, 0.0088,
     0.0086, 0.0084, 0.0082, 0.0080, 0.0078, 0.0076, 0.0074,
     0.0072, 0.0070, 0.0068, 0.0066, 0.0064, 0.0062, 0.0060,
     0.0058, 0.0056, 0.0054, 0.0052, 0.0050, 0.0048, 0.0046,
     0.0044, 0.0042, 0.0040, 0.0038, 0.0037, 0.0036, 0.0035,
     0.0034, 0.0033, 0.0032, 0.0031, 0.0030],
    dtype=float,
)
R0 = 1.0e-5
DR = 1.5e-3
K_MAX = 20.0
PAD_FACTOR = 32
# Four-decimal quadrature documented in evaluate_ball_spectral_constant.py.
J_BALL = 159.7247
I_WALL = (np.pi**2 / 3.0 - 2.0) / np.sqrt(2.0)
ELL_WALL = np.log(3.0 + 2.0 * np.sqrt(3.0)) / np.sqrt(2.0)
C_LAMBDA_3 = (
    4.0 * np.pi * (np.pi**2 / 3.0 - 2.0) * J_BALL
    / np.log(3.0 + 2.0 * np.sqrt(3.0))
)
A_R = 4.0 * np.pi * I_WALL * J_BALL
A_EPS = 4.0 * np.pi * (np.pi**2 - 6.0) * J_BALL

DATA_DIR = ROOT / "numerics"
FIG_DIR = ROOT / "figures"


def potential(phi: np.ndarray, eps: float) -> np.ndarray:
    return 0.25 * phi**2 * (phi**2 - 4.0 * (1.0 + eps) * phi + 4.0)


def potential_prime(phi: np.ndarray, eps: float) -> np.ndarray:
    return phi * (phi**2 - 3.0 * (1.0 + eps) * phi + 2.0)


def stationary_points(eps: float) -> tuple[float, float]:
    disc = np.sqrt(9.0 * (1.0 + eps) ** 2 - 8.0)
    return (3.0 * (1.0 + eps) - disc) / 2.0, (
        3.0 * (1.0 + eps) + disc
    ) / 2.0


def delta_w(eps: float) -> float:
    _, phi_tv = stationary_points(eps)
    return -float(potential(np.array(phi_tv), eps))


def thin_wall_radius(eps: float) -> float:
    sigma = 4.0 / (3.0 * np.sqrt(2.0))
    return 2.0 * sigma / delta_w(eps)


def shooting_seed(eps: float, r_max: float, mesh: np.ndarray) -> np.ndarray:
    phi_bar, phi_tv = stationary_points(eps)

    def integrate(phi_c: float, dense_output: bool = False):
        y0 = np.array(
            [
                phi_c + potential_prime(np.array(phi_c), eps) * R0**2 / 6.0,
                potential_prime(np.array(phi_c), eps) * R0 / 3.0,
            ]
        )

        def ode(r, y):
            return [y[1], potential_prime(y[0], eps) - 2.0 * y[1] / r]

        def crossed_zero(r, y):
            return y[0]

        crossed_zero.terminal = True
        crossed_zero.direction = -1

        def turned_around(r, y):
            return y[1]

        turned_around.terminal = True
        turned_around.direction = 1

        return solve_ivp(
            ode,
            (R0, r_max),
            y0,
            method="DOP853",
            rtol=2.0e-11,
            atol=2.0e-13,
            dense_output=dense_output,
            events=None if dense_output else (crossed_zero, turned_around),
            max_step=0.05,
        )

    low = phi_bar * (1.0 + 1.0e-10)
    high = np.nextafter(phi_tv, phi_bar)
    low_result = integrate(low)
    high_result = integrate(high)
    if len(low_result.t_events[1]) == 0:
        raise RuntimeError("The lower shooting bracket is not an undershoot")
    if len(high_result.t_events[0]) == 0:
        raise RuntimeError("The upper shooting bracket is not an overshoot")

    for _ in range(72):
        midpoint = 0.5 * (low + high)
        result = integrate(midpoint)
        if len(result.t_events[0]) > 0:
            high = midpoint
        else:
            low = midpoint
    phi_c = 0.5 * (low + high)
    result = integrate(phi_c, dense_output=True)
    guess = result.sol(mesh)
    guess[0] = np.maximum(guess[0], 0.0)
    return guess


def continued_seed(
    previous: dict[str, object],
    eps: float,
    phi_tv: float,
    r_tw: float,
    mesh: np.ndarray,
) -> np.ndarray:
    shift = r_tw - float(previous["r_tw"])
    amplitude = phi_tv / float(previous["phi_tv"])
    old_coordinate = mesh - shift
    guess = np.zeros((2, mesh.size), dtype=float)
    inside = old_coordinate < R0
    overlap = (old_coordinate >= R0) & (old_coordinate <= float(previous["r_max"]))
    guess[0, inside] = phi_tv
    guess[1, inside] = 0.0
    if np.any(overlap):
        values = previous["solution"].sol(old_coordinate[overlap])
        guess[:, overlap] = amplitude * values
    return guess


def solve_profile(eps: float, previous: dict[str, object] | None = None):
    phi_bar, phi_tv = stationary_points(eps)
    r_tw = thin_wall_radius(eps)
    r_max = max(25.0, r_tw + 18.0)

    left = np.linspace(R0, max(R0 + 1.0e-4, r_tw - 8.0), 250)
    wall = np.linspace(max(R0 + 2.0e-4, r_tw - 8.0), r_tw + 8.0, 900)
    right = np.linspace(r_tw + 8.0, r_max, 350)
    mesh = np.unique(np.concatenate([left, wall, right]))

    if previous is None:
        guess = shooting_seed(eps, r_max, mesh)
    else:
        guess = continued_seed(previous, eps, phi_tv, r_tw, mesh)

    def ode(r, y):
        return np.vstack(
            [y[1], potential_prime(y[0], eps) - 2.0 * y[1] / r]
        )

    def bc(ya, yb):
        origin = ya[1] - potential_prime(ya[0], eps) * R0 / 3.0
        outer = yb[1] + (np.sqrt(2.0) + 1.0 / r_max) * yb[0]
        return np.array([origin, outer])

    solution = solve_bvp(
        ode,
        bc,
        mesh,
        guess,
        tol=1.0e-8,
        max_nodes=300_000,
        verbose=0,
    )
    if not solution.success:
        raise RuntimeError(f"BVP failed for epsilon={eps}: {solution.message}")
    phi_c = float(solution.sol(R0)[0])
    phi_outer = float(solution.sol(r_max)[0])
    physical_center = phi_bar < phi_c <= phi_tv * (1.0 + 1.0e-10)
    physical_wall = abs(phi_outer) < 1.0e-6 * phi_tv
    if not (physical_center and physical_wall):
        raise RuntimeError(
            f"Nonphysical branch for epsilon={eps}: phi_c={phi_c}, "
            f"phi_outer={phi_outer}, expected phi_c near ({phi_bar}, {phi_tv}] "
            "and a decayed exterior"
        )
    continuation = {
        "solution": solution,
        "r_max": r_max,
        "r_tw": r_tw,
        "phi_tv": phi_tv,
    }
    return solution, r_max, r_tw, phi_tv, continuation


def level_radius(solution, phi_c: float, theta: float, r_max: float) -> float:
    target = phi_c * np.sqrt(theta)

    def residual(r):
        return float(solution.sol(r)[0] - target)

    return brentq(residual, R0, r_max, xtol=2.0e-12, rtol=2.0e-12)


def radial_grid(r_max: float) -> np.ndarray:
    count = int(np.ceil(r_max / DR))
    return np.linspace(0.0, count * DR, count + 1)


def spectral_entropy(r: np.ndarray, phi: np.ndarray) -> tuple[float, float, float]:
    physical_n = len(r) - 1
    transform_n = next_fast_len(PAD_FACTOR * physical_n)
    samples = np.zeros(transform_n, dtype=float)
    j = np.arange(1, physical_n + 1)
    samples[j - 1] = r[j] * phi[j]

    sine_sum = 0.5 * DR * dst(samples, type=1)
    frequencies = np.pi * np.arange(1, transform_n + 1) / (
        (transform_n + 1) * DR
    )
    keep = frequencies <= K_MAX
    k = frequencies[keep]
    transform = np.sqrt(2.0 / np.pi) * sine_sum[keep] / k
    transform_zero = np.sqrt(2.0 / np.pi) * simpson(r**2 * phi, x=r)
    q = np.clip((transform / transform_zero) ** 2, 0.0, 1.0)

    k_full = np.concatenate([[0.0], k])
    q_full = np.concatenate([[1.0], q])
    integrand = np.zeros_like(q_full)
    positive = q_full > 0.0
    integrand[positive] = -4.0 * np.pi * k_full[positive] ** 2 * q_full[
        positive
    ] * np.log(q_full[positive])
    s_k = float(simpson(integrand, x=k_full))
    return s_k, transform_zero, float(k[1] - k[0])


def analyze_profile(
    eps: float, previous: dict[str, object] | None = None
) -> tuple[dict[str, float], dict[str, np.ndarray], dict[str, object]]:
    solution, r_max, r_tw, phi_tv, continuation = solve_profile(eps, previous)
    r = radial_grid(r_max)
    phi = solution.sol(np.maximum(r, R0))[0]
    phi[0] = solution.sol(R0)[0] - potential_prime(
        solution.sol(R0)[0], eps
    ) * R0**2 / 6.0
    phi_c = float(phi[0])

    radii = {
        theta: level_radius(solution, phi_c, theta, r_max)
        for theta in (0.25, 0.50, 0.75)
    }
    radius = radii[0.50]
    ell = radii[0.25] - radii[0.75]
    lambda_3 = radius / ell

    q_x = np.clip((phi / phi_c) ** 2, 0.0, 1.0)
    h_x = np.zeros_like(q_x)
    positive_x = q_x > 0.0
    h_x[positive_x] = -4.0 * np.pi * r[positive_x] ** 2 * q_x[
        positive_x
    ] * np.log(q_x[positive_x])
    s_x = float(simpson(h_x, x=r))
    s_k, transform_zero, delta_k = spectral_entropy(r, phi)
    u_3 = s_x * s_k

    integration_points = [
        point for point in (radius - 8.0, radius, radius + 8.0)
        if R0 < point < r_max
    ]

    def kinetic_density(radius_value: float) -> float:
        derivative = float(solution.sol(radius_value)[1])
        return 2.0 * np.pi * radius_value**2 * derivative**2

    def potential_density(radius_value: float) -> float:
        field = float(solution.sol(radius_value)[0])
        return 4.0 * np.pi * radius_value**2 * float(potential(field, eps))

    kinetic = float(quad(
        kinetic_density, R0, r_max, points=integration_points,
        epsabs=1.0e-8, epsrel=2.0e-12, limit=500,
    )[0])
    potential_term = float(quad(
        potential_density, R0, r_max, points=integration_points,
        epsabs=1.0e-8, epsrel=2.0e-12, limit=500,
    )[0])
    derrick = abs(kinetic + 3.0 * potential_term) / (
        abs(kinetic) + 3.0 * abs(potential_term)
    )

    ratios = {
        "ratio_radius": radius / r_tw,
        "ratio_width": ell / ELL_WALL,
        "ratio_sx": s_x / (4.0 * np.pi * I_WALL * radius**2),
        "ratio_sk": radius**3 * s_k / J_BALL,
        "ratio_u_radius": radius * u_3 / A_R,
        "ratio_u_lambda": lambda_3 * u_3 / C_LAMBDA_3,
        "ratio_u_epsilon": u_3 / (A_EPS * eps),
    }
    row = {
        "epsilon": eps,
        "phi_tv": phi_tv,
        "phi_c": phi_c,
        "r_max": r_max,
        "nodes": int(solution.x.size),
        "delta_k": delta_k,
        "R_tw": r_tw,
        "R": radius,
        "ell": ell,
        "Lambda_3": lambda_3,
        "Sx_3": s_x,
        "Sk_3": s_k,
        "U_3": u_3,
        "transform_zero": transform_zero,
        "T": kinetic,
        "V": potential_term,
        "Delta_D": derrick,
        **ratios,
    }
    profile = {"r": r, "phi": phi, "phi_c": np.array(phi_c), "R": np.array(radius)}
    return row, profile, continuation


def save_csv(rows: list[dict[str, float]]) -> None:
    path = DATA_DIR / "thin_wall_observables.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def wall_profile(s: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + (np.sqrt(2.0) - 1.0) * np.exp(np.sqrt(2.0) * s))


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "font.size": 9,
            "axes.labelsize": 9,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def selected_indices(rows: list[dict[str, float]]) -> list[int]:
    targets = [0.050, 0.030, 0.015, 0.003]
    eps = np.array([row["epsilon"] for row in rows])
    return [int(np.argmin(abs(eps - target))) for target in targets]


def plot_profiles(rows, profiles) -> None:
    configure_plotting()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, 4))
    indices = selected_indices(rows)
    for color, index in zip(colors, indices):
        row = rows[index]
        profile = profiles[index]
        r = profile["r"]
        normalized = profile["phi"] / float(profile["phi_c"])
        label = rf"$\epsilon={row['epsilon']:.3f}$"
        axes[0].plot(r / row["R"], normalized, color=color, lw=1.5, label=label)

        s = r - row["R"]
        mask = (s >= -5.0) & (s <= 5.0)
        axes[1].plot(s[mask], normalized[mask], color=color, lw=1.5, label=label)

    y = np.linspace(0.0, 1.4, 400)
    axes[0].plot(y, (y <= 1.0).astype(float), "k--", lw=1.0, label="sharp ball")
    s_ref = np.linspace(-5.0, 5.0, 500)
    axes[1].plot(s_ref, wall_profile(s_ref), "k--", lw=1.0, label=r"$F_{\rm wall}$")

    axes[0].set(xlabel=r"$r/R$", ylabel=r"$\phi(r)/\phi_c$", xlim=(0.0, 1.35), ylim=(-0.03, 1.05))
    axes[1].set(xlabel=r"$s=r-R$", ylabel=r"$\phi(R+s)/\phi_c$", xlim=(-5.0, 5.0), ylim=(-0.03, 1.05))
    axes[0].set_title("(a)", loc="left", pad=7)
    axes[1].set_title("(b)", loc="left", pad=7)
    axes[0].legend(frameon=False, loc="lower left")
    axes[1].legend(frameon=False, loc="upper right")
    for suffix in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"critical_bubble_profile_collapse.{suffix}", bbox_inches="tight")
    plt.close(fig)


def plot_ratios(rows) -> None:
    configure_plotting()
    lambda_3 = np.array([row["Lambda_3"] for row in rows])
    order = np.argsort(lambda_3)
    lambda_3 = lambda_3[order]
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.45), constrained_layout=True)

    series = [
        (("ratio_radius", r"$R/R_{\rm TW}$"), ("ratio_width", r"$\ell/\ell_{\rm wall}$")),
        (("ratio_sx", r"$\mathcal{R}_{x}$"), ("ratio_sk", r"$\mathcal{R}_{k}$")),
        (("ratio_u_radius", r"$\mathcal{R}_{U}$"), ("ratio_u_lambda", r"$\mathcal{R}_{\Lambda}$")),
    ]
    markers = ("o", "s")
    colors = ("#31688e", "#d1495b")
    for panel, pairs in enumerate(series):
        ax = axes[panel]
        ax.axhline(1.0, color="0.25", lw=0.9, ls="--")
        for marker, color, (key, label) in zip(markers, colors, pairs):
            values = np.array([rows[i][key] for i in order])
            ax.plot(lambda_3, values, marker=marker, ms=3.2, lw=1.1, color=color, label=label)
        ax.set_xlabel(r"$\Lambda_{3}$")
        ax.legend(frameon=False)
        ax.set_title(f"({chr(97 + panel)})", loc="left", pad=7)
    axes[0].set_ylabel("normalized ratio")
    for suffix in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"thin_wall_ratio_convergence.{suffix}", bbox_inches="tight")
    plt.close(fig)


def plot_scale_law(rows) -> dict[str, float]:
    configure_plotting()
    lambda_3 = np.array([row["Lambda_3"] for row in rows])
    u_3 = np.array([row["U_3"] for row in rows])
    order = np.argsort(lambda_3)
    lambda_3 = lambda_3[order]
    u_3 = u_3[order]
    fit_count = min(6, len(lambda_3))
    fit_x = np.log(lambda_3[-fit_count:])
    fit_y = np.log(u_3[-fit_count:])
    slope, intercept = np.polyfit(fit_x, fit_y, 1)
    residual = fit_y - (slope * fit_x + intercept)
    slope_error = np.sqrt(
        np.sum(residual**2) / (fit_count - 2) / np.sum((fit_x - fit_x.mean()) ** 2)
    )

    tail = lambda_3 >= 20.0
    tail_design = np.column_stack(
        [
            np.ones(np.count_nonzero(tail)),
            np.log(lambda_3[tail]) / lambda_3[tail],
            1.0 / lambda_3[tail],
        ]
    )
    tail_ratio = (
        lambda_3[tail] * u_3[tail] / C_LAMBDA_3
    )
    tail_coefficients = np.linalg.lstsq(
        tail_design, tail_ratio, rcond=None
    )[0]
    tail_residual = tail_ratio - tail_design @ tail_coefficients
    tail_rmse = float(np.sqrt(np.mean(tail_residual**2)))
    limit, log_coefficient, inverse_coefficient = tail_coefficients

    fig, ax = plt.subplots(figsize=(3.45, 2.8), constrained_layout=True)
    ax.loglog(lambda_3, u_3, "o", ms=4.0, color="#31688e", label="numerical bubbles")
    line_x = np.geomspace(lambda_3.min() * 0.9, lambda_3.max() * 1.12, 300)
    ax.loglog(line_x, C_LAMBDA_3 / line_x, "k--", lw=1.2, label=r"$C_{\Lambda,3}/\Lambda_3$")
    correction_x = np.geomspace(20.0, lambda_3.max(), 200)
    correction_ratio = (
        limit
        + log_coefficient * np.log(correction_x) / correction_x
        + inverse_coefficient / correction_x
    )
    ax.loglog(
        correction_x,
        C_LAMBDA_3 * correction_ratio / correction_x,
        color="#d1495b",
        ls=":",
        lw=1.3,
        label="tail-motivated fit",
    )
    ax.set(xlabel=r"$\Lambda_{3}$", ylabel=r"$\mathcal{U}_{3}$")
    ax.set_xticks([5.0, 10.0, 20.0, 40.0, 60.0], labels=["5", "10", "20", "40", "60"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(frameon=False)
    for suffix in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"dual_entropy_scale_separation.{suffix}", bbox_inches="tight")
    plt.close(fig)
    return {
        "fit_count": fit_count,
        "fit_slope": float(slope),
        "fit_slope_error": float(slope_error),
        "fit_intercept": float(intercept),
        "tail_fit_minimum_Lambda": 20.0,
        "tail_fit_limit": float(limit),
        "tail_fit_log_coefficient": float(log_coefficient),
        "tail_fit_inverse_coefficient": float(inverse_coefficient),
        "tail_fit_rmse": tail_rmse,
    }


def main() -> None:
    DATA_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(exist_ok=True)
    rows = []
    profiles = []
    continuation = None
    for eps in EPSILONS:
        print(f"Solving epsilon={eps:.4f}", flush=True)
        row, profile, continuation = analyze_profile(float(eps), continuation)
        rows.append(row)
        profiles.append(profile)
        print(
            f"  R={row['R']:.8f}, ell={row['ell']:.8f}, "
            f"U={row['U_3']:.8f}, Delta_D={row['Delta_D']:.3e}",
            flush=True,
        )

    save_csv(rows)
    np.savez_compressed(
        DATA_DIR / "thin_wall_profiles.npz",
        **{
            f"eps_{row['epsilon']:.4f}_{key}": value
            for row, profile in zip(rows, profiles)
            for key, value in profile.items()
        },
    )
    plot_profiles(rows, profiles)
    plot_ratios(rows)
    fit = plot_scale_law(rows)
    summary = {
        "constants": {
            "I_wall": I_WALL,
            "ell_wall": ELL_WALL,
            "J_ball": J_BALL,
            "C_Lambda_3": C_LAMBDA_3,
            "A_R": A_R,
            "A_epsilon": A_EPS,
        },
        "fit": fit,
        "smallest_epsilon": rows[-1],
    }
    (DATA_DIR / "thin_wall_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
