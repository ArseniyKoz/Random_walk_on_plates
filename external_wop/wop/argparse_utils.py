from __future__ import annotations

import argparse

import numpy as np


def parse_vec3_arg(text: str, example: str = "3 0 0") -> np.ndarray:
    parts = text.replace(",", " ").split()
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            f"Expected exactly three numbers, e.g. '{example}'."
        )
    try:
        return np.asarray([float(p) for p in parts], dtype=float)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
