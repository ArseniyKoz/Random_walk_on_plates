from __future__ import annotations

import numpy as np
import pytest

from wop.geometry import make_axis_aligned_box, orient_normals
from wop.wop import estimate_wop
from wop.wos_box import estimate_wos_box


BOX_MIN = np.array([-1.0, -1.0, -1.0], dtype=float)
BOX_MAX = np.array([1.0, 1.0, 1.0], dtype=float)
INTERIOR = np.array([0.0, 0.0, 0.0], dtype=float)
A_REF = np.array([0.2, -0.1, 0.3], dtype=float)

# Unified with notebooks: canonical start point + deterministic seed.
X0_CASES = [np.array([3.0, 0.0, 0.0], dtype=float)]
SEED = 314159
N_PATHS = 1000
MAX_STEPS = 100000


def _exact_u(x: np.ndarray) -> float:
    return float(1.0 / np.linalg.norm(x - A_REF))


def _boundary_f(y: np.ndarray, _face: int | None) -> float:
    return _exact_u(y)


@pytest.mark.parametrize("x0", X0_CASES, ids=["x0_3_0_0"])
def test_wop_agrees_with_wos_for_box_exterior(x0: np.ndarray) -> None:
    poly = orient_normals(make_axis_aligned_box(BOX_MIN, BOX_MAX), INTERIOR)
    L = poly.characteristic_length()
    delta = 1e-3 * L

    delta_min = float(np.min(np.concatenate([A_REF - BOX_MIN, BOX_MAX - A_REF])))
    if delta_min <= delta:
        raise AssertionError("Configuration must satisfy delta_min > delta.")
    lipschitz_bound = 1.0 / (delta_min - delta) ** 2
    bias_delta = lipschitz_bound * delta

    rng_wop = np.random.default_rng(SEED)
    rng_wos = np.random.default_rng(SEED + 1)

    res_wop = estimate_wop(
        poly=poly,
        x0=x0,
        boundary_f=_boundary_f,
        n_paths=N_PATHS,
        rng=rng_wop,
        max_steps=MAX_STEPS,
        u_inf=0.0,
        r_max=1e6,
    )
    res_wos = estimate_wos_box(
        x0=x0,
        box_min=BOX_MIN,
        box_max=BOX_MAX,
        boundary_f=_boundary_f,
        n_paths=N_PATHS,
        rng=rng_wos,
        delta=delta,
        rho_scale=1.0,
        rho1_scale=2.0,
        max_steps=MAX_STEPS,
        u_inf=0.0,
    )

    assert res_wop.n_total == N_PATHS
    assert res_wos.n_total == N_PATHS
    assert 0 <= res_wop.n_truncated <= res_wop.n_total
    assert 0 <= res_wos.n_truncated <= res_wos.n_total

    assert np.isclose(
        res_wop.eps,
        3.0 * np.sqrt(res_wop.S2 / res_wop.n_total),
        rtol=0.0,
        atol=1e-12,
    )
    assert np.isclose(
        res_wos.eps,
        3.0 * np.sqrt(res_wos.S2 / res_wos.n_total),
        rtol=0.0,
        atol=1e-12,
    )

    diff = abs(res_wop.J - res_wos.J)
    se_wop = np.sqrt(res_wop.S2 / res_wop.n_total)
    se_wos = np.sqrt(res_wos.S2 / res_wos.n_total)
    combined = np.sqrt(se_wop**2 + se_wos**2 + bias_delta**2)

    assert diff <= 4.0 * combined + 2e-3

    exact = _exact_u(x0)
    err_wop = abs(res_wop.J - exact)
    err_wos = abs(res_wos.J - exact)
    assert err_wop <= 2.0 * res_wop.eps + 3e-3
    assert err_wos <= 2.0 * res_wos.eps + bias_delta + 3e-3
