"""Independent quadrature of the sharp three-dimensional ball entropy.

The tail estimate uses the period average of the leading large-z integrand.
Taking Z at integer multiples of pi suppresses its leading oscillatory
endpoint correction. The convergence of both Z and the Simpson spacing is
reported so the four-decimal value used in the manuscript can be checked.
"""

from __future__ import annotations

from math import cos, log, pi, sin


def entropy_integrand(z: float) -> float:
    if z == 0.0:
        return 0.0
    if z < 0.03:
        q = 1.0 - z * z / 5.0 + 3.0 * z**4 / 175.0
    else:
        amplitude = 3.0 * (sin(z) - z * cos(z)) / z**3
        q = amplitude * amplitude
    return -4.0 * pi * z * z * q * log(q) if q > 0.0 else 0.0


def truncated_integral(z_max: float, target_step: float) -> float:
    intervals = round(z_max / target_step)
    intervals += intervals % 2
    step = z_max / intervals
    total = entropy_integrand(0.0) + entropy_integrand(z_max)
    for index in range(1, intervals):
        total += (4 if index % 2 else 2) * entropy_integrand(index * step)
    return total * step / 3.0


def leading_tail(z_max: float) -> float:
    return 36.0 * pi * (2.0 * log(z_max) + 1.5 - log(1.5)) / z_max


def main() -> None:
    estimates = []
    for periods in (500, 1000, 2000):
        z_max = periods * pi
        estimate = truncated_integral(z_max, 0.005) + leading_tail(z_max)
        estimates.append(estimate)
        print(f"Z={periods} pi: J_ball={estimate:.9f}")
    finer = truncated_integral(1000 * pi, 0.0025) + leading_tail(1000 * pi)
    print(f"Z=1000 pi, half spacing: J_ball={finer:.9f}")
    if max(estimates + [finer]) - min(estimates + [finer]) > 5.0e-5:
        raise RuntimeError("The displayed four-decimal value is not stable")
    print(f"Value at four decimal places: {estimates[-1]:.4f}")


if __name__ == "__main__":
    main()
