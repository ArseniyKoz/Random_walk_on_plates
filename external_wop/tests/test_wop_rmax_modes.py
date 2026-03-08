from __future__ import annotations

import numpy as np
import pytest

from wop.geometry import make_axis_aligned_box, orient_normals
from wop.wop import _resolve_r_max_projection, trace_wop_trajectory


def _make_unit_box_poly():
    poly = make_axis_aligned_box(
        np.array([-1.0, -1.0, -1.0], dtype=float),
        np.array([1.0, 1.0, 1.0], dtype=float),
    )
    return orient_normals(poly, np.array([0.0, 0.0, 0.0], dtype=float))


def _boundary_const(y: np.ndarray, _face: int | None) -> float:
    return float(np.sum(y) * 0.0 + 1.0)


def test_escape_mode_keeps_legacy_escape_behavior() -> None:
    poly = _make_unit_box_poly()
    rng = np.random.default_rng(20260303)

    tr = trace_wop_trajectory(
        poly=poly,
        x0=np.array([3.0, 0.0, 0.0], dtype=float),
        boundary_f=_boundary_const,
        rng=rng,
        eps_in=1e-12,
        eps_plane=1e-12,
        min_abs_denom=1e-3,
        max_steps=1000,
        u_inf=0.0,
        r_max=2.0,
        r_max_mode="escape",
        r_max_factor=3.0,
    )

    assert tr.status == "escaped"
    assert tr.steps == 0


def test_project_mode_does_not_escape_on_rmax_crossing() -> None:
    poly = _make_unit_box_poly()
    rng = np.random.default_rng(20260303)

    tr = trace_wop_trajectory(
        poly=poly,
        x0=np.array([3.0, 0.0, 0.0], dtype=float),
        boundary_f=_boundary_const,
        rng=rng,
        eps_in=1e-12,
        eps_plane=1e-12,
        min_abs_denom=1e-3,
        max_steps=1000,
        u_inf=0.0,
        r_max=2.0,
        r_max_mode="project",
        r_max_factor=3.0,
    )

    assert tr.status != "escaped"


def test_project_mode_autormax_is_enabled_when_rmax_is_none() -> None:
    poly = _make_unit_box_poly()
    rng = np.random.default_rng(20260303)

    tr = trace_wop_trajectory(
        poly=poly,
        x0=np.array([4.0, 0.5, 0.0], dtype=float),
        boundary_f=_boundary_const,
        rng=rng,
        eps_in=1e-12,
        eps_plane=1e-12,
        min_abs_denom=1e-3,
        max_steps=1000,
        u_inf=0.0,
        r_max=None,
        r_max_mode="project",
        r_max_factor=3.0,
    )

    assert tr.status != "escaped"


def test_manual_rmax_overrides_auto_rmax() -> None:
    poly = _make_unit_box_poly()
    x0 = np.array([4.0, 0.0, 0.0], dtype=float)

    auto_cfg = _resolve_r_max_projection(
        poly=poly,
        x0=x0,
        r_max=None,
        r_max_mode="project",
        r_max_factor=3.0,
    )
    manual_cfg = _resolve_r_max_projection(
        poly=poly,
        x0=x0,
        r_max=2.5,
        r_max_mode="project",
        r_max_factor=3.0,
    )

    assert auto_cfg.enabled is True
    assert manual_cfg.enabled is True
    assert manual_cfg.radius == pytest.approx(2.5)
    assert auto_cfg.radius > manual_cfg.radius


def test_invalid_rmax_factor_raises() -> None:
    poly = _make_unit_box_poly()
    rng = np.random.default_rng(20260303)

    with pytest.raises(ValueError, match="r_max_factor"):
        trace_wop_trajectory(
            poly=poly,
            x0=np.array([3.0, 0.0, 0.0], dtype=float),
            boundary_f=_boundary_const,
            rng=rng,
            eps_in=1e-12,
            eps_plane=1e-12,
            min_abs_denom=1e-3,
            max_steps=1000,
            u_inf=0.0,
            r_max=None,
            r_max_mode="project",
            r_max_factor=1.0,
        )
