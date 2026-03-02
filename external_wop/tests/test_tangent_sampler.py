from __future__ import annotations

import numpy as np

from wop.sampling import sample_tangent_direction


def _normalize(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def test_tangent_sampler_produces_unit_tangent_vectors() -> None:
    rng = np.random.default_rng(20260225)

    normals = [
        np.array([0.0, 0.0, 1.0], dtype=float),
        _normalize(np.array([1.0, -2.0, 0.5], dtype=float)),
        _normalize(np.array([-0.7, 0.1, 2.0], dtype=float)),
    ]

    for nu in normals:
        for _ in range(2000):
            w = sample_tangent_direction(nu=nu, rng=rng)
            assert np.isclose(np.linalg.norm(w), 1.0, atol=1e-12, rtol=0.0)
            assert abs(float(np.dot(w, nu))) <= 1e-12


def test_tangent_sampler_isotropic_second_moment() -> None:
    rng = np.random.default_rng(20260225)
    nu = _normalize(np.array([1.0, -2.0, 0.5], dtype=float))

    n = 50000
    samples = np.empty((n, 3), dtype=float)
    for k in range(n):
        samples[k] = sample_tangent_direction(nu=nu, rng=rng)

    mean = np.mean(samples, axis=0)
    second_moment = (samples.T @ samples) / n
    target = 0.5 * (np.eye(3, dtype=float) - np.outer(nu, nu))

    assert np.all(np.abs(mean) <= 1.5e-2)
    assert np.all(np.abs(second_moment - target) <= 2.5e-2)
