from __future__ import annotations

import numpy as np

from wop.geometry import make_axis_aligned_box, orient_normals
from wop.wop import estimate_wop


def test_wop_matches_exact_harmonic_solution_outside_box() -> None:
    box_min = np.array([-1.0, -1.0, -1.0], dtype=float)
    box_max = np.array([1.0, 1.0, 1.0], dtype=float)
    interior = np.array([0.0, 0.0, 0.0], dtype=float)

    poly = orient_normals(make_axis_aligned_box(box_min, box_max), interior)
    a = np.array([0.2, -0.1, 0.3], dtype=float)

    def exact_u(x: np.ndarray) -> float:
        return float(1.0 / np.linalg.norm(x - a))

    def boundary_f(y: np.ndarray, _face: int | None) -> float:
        return exact_u(y)

    rng = np.random.default_rng(314159)

    x0_list = [
        np.array([3.0, 0.0, 0.0], dtype=float),
        np.array([0.0, 0.0, 4.0], dtype=float),
        np.array([2.5, 1.7, 0.2], dtype=float),
    ]

    for x0 in x0_list:
        result = estimate_wop(
            poly=poly,
            x0=x0,
            boundary_f=boundary_f,
            n_paths=4000,
            rng=rng,
            max_steps=200000,
            u_inf=0.0,
            r_max=1e6,
        )
        assert result.n_total == 4000
        assert 0 <= result.n_truncated <= result.n_total
        assert result.S2 >= 0.0
        assert np.isclose(result.eps, 3.0 * np.sqrt(result.S2 / result.n_total))

        exact = exact_u(x0)
        err = abs(result.J - exact)

        assert err <= 1.2 * result.eps + 2e-3
