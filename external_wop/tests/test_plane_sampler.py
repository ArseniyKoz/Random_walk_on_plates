from __future__ import annotations

import numpy as np

from wop.sampling import sample_hit_on_plane_from_point


def _theoretical_quantile_r(p: float, h: float) -> float:
    return float(h * np.sqrt(1.0 / (1.0 - p) ** 2 - 1.0))


def test_plane_sampler_quantiles_match_poisson_kernel() -> None:
    rng = np.random.default_rng(20260216)

    h = 1.7
    x = np.array([0.0, 0.0, h], dtype=float)
    nu = np.array([0.0, 0.0, 1.0], dtype=float)
    b = 0.0

    n = 60000
    radial = np.empty(n, dtype=float)

    for k in range(n):
        y, _, _ = sample_hit_on_plane_from_point(x=x, nu=nu, b=b, rng=rng)
        radial[k] = np.hypot(y[0], y[1])

    probs = np.array([0.25, 0.5, 0.75, 0.9], dtype=float)
    empirical = np.quantile(radial, probs)
    theoretical = np.array([_theoretical_quantile_r(float(p), h) for p in probs], dtype=float)

    # Quantiles are robust and avoid external dependencies for full KS tests.
    tol = 0.06 * theoretical + 0.02
    assert np.all(np.abs(empirical - theoretical) <= tol)
