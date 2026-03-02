from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class TrajectoryResult:
    value: float
    steps: int
    status: str  # hit_face | timeout | escaped


@dataclass(frozen=True)
class EstimateResult:
    J: float
    S2: float
    eps: float
    n_total: int
    n_truncated: int
    mean_steps: float


class _EstimateAccumulator:
    def __init__(self) -> None:
        self.j_sum = 0.0
        self.s2_raw_sum = 0.0
        self.steps_sum = 0
        self.n_timeout = 0
        self.n_escaped = 0

    def add(self, tr: TrajectoryResult) -> None:
        ksi = float(tr.value)
        self.j_sum += ksi
        self.s2_raw_sum += ksi * ksi
        self.steps_sum += int(tr.steps)

        if tr.status == "timeout":
            self.n_timeout += 1
        elif tr.status == "escaped":
            self.n_escaped += 1

    def finalize(self, n_total: int) -> EstimateResult:
        if n_total <= 0:
            raise ValueError("n_total must be positive.")

        j = self.j_sum / n_total
        s2 = self.s2_raw_sum / n_total - j * j
        s2 = max(float(s2), 0.0)
        eps = 3.0 * float(np.sqrt(s2 / n_total))
        n_truncated = self.n_timeout + self.n_escaped

        return EstimateResult(
            J=j,
            S2=s2,
            eps=eps,
            n_total=n_total,
            n_truncated=n_truncated,
            mean_steps=float(self.steps_sum / n_total),
        )


def estimate_from_trajectories(
    n_paths: int,
    trace_once: Callable[[], TrajectoryResult],
) -> EstimateResult:
    """Generic Monte Carlo aggregation for trajectory-based estimators."""
    if n_paths <= 0:
        raise ValueError("n_paths must be positive.")

    acc = _EstimateAccumulator()
    add = acc.add
    for _ in range(n_paths):
        add(trace_once())
    return acc.finalize(n_paths)
