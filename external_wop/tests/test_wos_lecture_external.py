from __future__ import annotations

import numpy as np

from wop.wos_box import estimate_wos_box, trace_wos_box_trajectory


def test_trajectory_stops_on_boundary_delta_neighborhood() -> None:
    box_min = np.array([-1.0, -1.0, -1.0], dtype=float)
    box_max = np.array([1.0, 1.0, 1.0], dtype=float)
    x0 = np.array([3.0, 0.0, 0.0], dtype=float)
    rng = np.random.default_rng(20260304)

    tr = trace_wos_box_trajectory(
        x0=x0,
        box_min=box_min,
        box_max=box_max,
        boundary_f=lambda _y, _face: 1.0,
        rng=rng,
        delta=1e-3,
        rho_scale=2.0,
        rho1_scale=4.0,
        max_steps=100000,
        u_inf=0.0,
    )

    assert tr.status == "hit_face"
    assert tr.steps > 0


def test_weighted_estimator_is_not_constant_one_for_far_start() -> None:
    box_min = np.array([-1.0, -1.0, -1.0], dtype=float)
    box_max = np.array([1.0, 1.0, 1.0], dtype=float)
    x0 = np.array([10.0, 0.0, 0.0], dtype=float)
    rng = np.random.default_rng(20260305)

    res = estimate_wos_box(
        x0=x0,
        box_min=box_min,
        box_max=box_max,
        boundary_f=lambda _y, _face: 1.0,
        n_paths=2000,
        rng=rng,
        delta=1e-3,
        rho_scale=2.0,
        rho1_scale=4.0,
        max_steps=200000,
        u_inf=0.0,
    )

    assert 0 <= res.n_truncated <= res.n_total
    assert res.J < 0.95


def test_harmonic_reference_solution_is_matched() -> None:
    box_min = np.array([-1.0, -1.0, -1.0], dtype=float)
    box_max = np.array([1.0, 1.0, 1.0], dtype=float)
    x0 = np.array([3.0, 0.0, 0.0], dtype=float)
    a = np.array([0.2, -0.1, 0.3], dtype=float)

    def exact_u(x: np.ndarray) -> float:
        return float(1.0 / np.linalg.norm(x - a))

    rng = np.random.default_rng(20260306)
    res = estimate_wos_box(
        x0=x0,
        box_min=box_min,
        box_max=box_max,
        boundary_f=lambda y, _face: exact_u(y),
        n_paths=3000,
        rng=rng,
        delta=1e-3,
        rho_scale=2.0,
        rho1_scale=4.0,
        max_steps=200000,
        u_inf=0.0,
    )

    exact = exact_u(x0)
    err = abs(res.J - exact)
    assert err <= 2.5 * res.eps + 1e-2
