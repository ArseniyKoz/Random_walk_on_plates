from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from wop.geometry import make_axis_aligned_box, orient_normals
from wop.wop import estimate_wop


def _find_cpp_cli() -> Path | None:
    env_path = os.environ.get("WOP_CPP_CLI")
    if env_path:
        p = Path(env_path).expanduser().resolve()
        if p.exists():
            return p
        return None

    py_root = Path(__file__).resolve().parents[1]
    workspace_root = py_root.parent
    cpp_root = workspace_root / "external_wop_cpp"
    candidates = [
        cpp_root / "build" / "wop_cli",
        cpp_root / "build" / "wop_cli.exe",
        cpp_root / "build" / "Release" / "wop_cli.exe",
        cpp_root / "build" / "Debug" / "wop_cli.exe",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


CPP_CLI = _find_cpp_cli()
pytestmark = pytest.mark.skipif(
    CPP_CLI is None,
    reason="C++ CLI not found. Build external_wop_cpp and set WOP_CPP_CLI.",
)


def _run_cpp_cli(x0: np.ndarray, seed: int, n_paths: int = 4000, max_steps: int = 200000, r_max: float = 1e6) -> dict[str, float]:
    assert CPP_CLI is not None
    cmd = [
        str(CPP_CLI),
        "--example",
        "box",
        "--x0",
        f"{float(x0[0])} {float(x0[1])} {float(x0[2])}",
        "--n",
        str(int(n_paths)),
        "--seed",
        str(int(seed)),
        "--max-steps",
        str(int(max_steps)),
        "--r-max",
        str(float(r_max)),
        "--json",
    ]
    out = subprocess.check_output(cmd, text=True)
    return json.loads(out)


def test_cpp_cli_reproducible_for_fixed_seed() -> None:
    x0 = np.array([3.0, 0.0, 0.0], dtype=float)
    run1 = _run_cpp_cli(x0=x0, seed=314159, n_paths=2000)
    run2 = _run_cpp_cli(x0=x0, seed=314159, n_paths=2000)
    assert run1 == run2


def test_cpp_cli_matches_python_estimator_statistics() -> None:
    box_min = np.array([-1.0, -1.0, -1.0], dtype=float)
    box_max = np.array([1.0, 1.0, 1.0], dtype=float)
    interior = np.array([0.0, 0.0, 0.0], dtype=float)
    poly = orient_normals(make_axis_aligned_box(box_min, box_max), interior)

    a = np.array([0.2, -0.1, 0.3], dtype=float)

    def exact_u(x: np.ndarray) -> float:
        return float(1.0 / np.linalg.norm(x - a))

    def boundary_f(y: np.ndarray, _face: int | None) -> float:
        return exact_u(y)

    x0_list = [
        np.array([3.0, 0.0, 0.0], dtype=float),
        np.array([2.0, 2.0, 2.0], dtype=float),
    ]

    for idx, x0 in enumerate(x0_list):
        seed = 10_000 + idx
        rng = np.random.default_rng(seed)
        py_res = estimate_wop(
            poly=poly,
            x0=x0,
            boundary_f=boundary_f,
            n_paths=4000,
            rng=rng,
            max_steps=200000,
            u_inf=0.0,
            r_max=1e6,
        )
        cpp_res = _run_cpp_cli(x0=x0, seed=seed, n_paths=4000, max_steps=200000, r_max=1e6)

        assert int(cpp_res["n_total"]) == 4000
        assert 0 <= int(cpp_res["n_truncated"]) <= int(cpp_res["n_total"])
        assert float(cpp_res["S2"]) >= 0.0
        assert np.isclose(
            float(cpp_res["eps"]),
            3.0 * np.sqrt(float(cpp_res["S2"]) / float(cpp_res["n_total"])),
            rtol=0.0,
            atol=1e-12,
        )

        diff = abs(float(cpp_res["J"]) - py_res.J)
        se_cpp = np.sqrt(float(cpp_res["S2"]) / float(cpp_res["n_total"]))
        se_py = np.sqrt(py_res.S2 / py_res.n_total)
        combined = np.sqrt(se_cpp**2 + se_py**2)

        assert diff <= 4.0 * combined + 3e-3

        exact = exact_u(x0)
        assert abs(float(cpp_res["J"]) - exact) <= 1.5 * float(cpp_res["eps"]) + 3e-3
