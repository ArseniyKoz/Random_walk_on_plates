from __future__ import annotations

import numpy as np

from wop.wop import _scan_signed_distances


def test_scan_signed_distances_outside_and_minima() -> None:
    d = np.array([-0.2, 0.3, -0.1, 0.05], dtype=float)

    info = _scan_signed_distances(d=d, eps_in=0.0, target_idx=3, eps_plane=1e-3)

    assert info.any_outside is True
    assert info.all_inside is False
    assert info.argmin_outside == 3
    assert info.argmin_abs == 3
    assert info.hit_target_face is False


def test_scan_signed_distances_hit_on_target_face() -> None:
    d = np.array([-0.2, -1e-5, -0.1, -1e-6], dtype=float)

    info = _scan_signed_distances(d=d, eps_in=0.0, target_idx=1, eps_plane=1e-4)

    assert info.any_outside is False
    assert info.all_inside is True
    assert info.argmin_outside is None
    assert info.argmin_abs == 3
    assert info.hit_target_face is True
