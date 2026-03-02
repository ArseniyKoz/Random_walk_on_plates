from __future__ import annotations

import numpy as np

from wop.geometry import make_axis_aligned_box, orient_normals
from wop.wop import estimate_wop
from wop.wos_box import estimate_wos_box


def test_wop_agrees_with_wos_for_box_exterior() -> None:
    box_min = np.array([-1.0, -1.0, -1.0], dtype=float)
    box_max = np.array([1.0, 1.0, 1.0], dtype=float)
    interior = np.array([0.0, 0.0, 0.0], dtype=float)

    poly = orient_normals(make_axis_aligned_box(box_min, box_max), interior)
    L = poly.characteristic_length()

    a = np.array([0.2, -0.1, 0.3], dtype=float)

    def exact_u(x: np.ndarray) -> float:
        return float(1.0 / np.linalg.norm(x - a))

    def boundary_f(y: np.ndarray, _face: int | None) -> float:
        return exact_u(y)

    eps = 1e-3 * L
    delta_min = float(np.min(np.concatenate([a - box_min, box_max - a])))
    if delta_min <= eps:
        raise AssertionError("Configuration must satisfy delta_min > eps.")
    lipschitz_bound = 1.0 / (delta_min - eps) ** 2
    bias_eps = lipschitz_bound * eps

    x0_list = [
        np.array([3.0, 0.0, 0.0], dtype=float),
        np.array([2.0, 2.0, 2.0], dtype=float),
    ]

    for idx, x0 in enumerate(x0_list):
        rng_wop = np.random.default_rng(10_000 + idx)
        rng_wos = np.random.default_rng(20_000 + idx)

        res_wop = estimate_wop(
            poly=poly,
            x0=x0,
            boundary_f=boundary_f,
            n_paths=4000,
            rng=rng_wop,
            max_steps=200000,
            u_inf=0.0,
            r_max=1e6,
        )
        res_wos = estimate_wos_box(
            x0=x0,
            box_min=box_min,
            box_max=box_max,
            boundary_f=boundary_f,
            n_paths=4000,
            rng=rng_wos,
            eps=eps,
            max_steps=200000,
            u_inf=0.0,
            r_max=1e6,
        )
        assert res_wop.n_total == 4000
        assert res_wos.n_total == 4000
        assert 0 <= res_wop.n_truncated <= res_wop.n_total
        assert 0 <= res_wos.n_truncated <= res_wos.n_total
        assert np.isclose(res_wop.eps, 3.0 * np.sqrt(res_wop.S2 / res_wop.n_total))
        assert np.isclose(res_wos.eps, 3.0 * np.sqrt(res_wos.S2 / res_wos.n_total))

        diff = abs(res_wop.J - res_wos.J)
        se_wop = np.sqrt(res_wop.S2 / res_wop.n_total)
        se_wos = np.sqrt(res_wos.S2 / res_wos.n_total)
        combined = np.sqrt(se_wop**2 + se_wos**2 + bias_eps**2)

        assert diff <= 3.0 * combined + 2e-3
