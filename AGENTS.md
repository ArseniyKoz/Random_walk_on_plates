Always first write tests then implement features.

# Repository Guidelines

## Project Structure & Module Organization
Primary code lives in `external_wop/`. The Python package is in `external_wop/wop/`:
- `geometry.py`, `sampling.py`, `wop.py`, `wos_box.py`: core math and Monte Carlo algorithms.
- `estimation.py`: trajectory aggregation/statistics.
- `cli.py`: command-line entry point.

Tests are in `external_wop/tests/` (`test_*.py`). Supporting materials:
- `external_wop/docs/`: derivations and usage notes.
- `external_wop/try/`: exploratory scripts (not part of the package API).
- `external_wop/pyproject.toml`: dependencies, packaging, and pytest config.

## Build, Test, and Development Commands
Run commands from `external_wop/`.

```bash
python -m pip install -e .[test]
```
Install package in editable mode with test dependencies.

```bash
python -m pip install -e .[viz]
```
Install optional plotting dependencies.

```bash
pytest -q
```
Run the full test suite (configured via `tool.pytest.ini_options`).

```bash
python -m wop.cli --example box --x0 "3 0 0" --n 50000
```
Run a reproducible CLI example and inspect estimator diagnostics.

## Coding Style & Naming Conventions
Target Python `>=3.10`. Follow PEP 8 with 4-space indentation and explicit type hints (including NumPy types where practical). Use:
- `snake_case` for functions/variables/modules
- `PascalCase` for classes/dataclasses
- clear, validated numeric inputs (shape/type checks for vectors and parameters)

Keep algorithmic logic in `wop/` modules; keep CLI/parsing concerns in `cli.py` and `argparse_utils.py`.

## Testing Guidelines
Use `pytest` and place tests under `external_wop/tests/` with names `test_*.py` and functions `test_*`. Favor deterministic tests with fixed seeds (for example, `np.random.default_rng(314159)`). For numerical methods, assert tolerances and invariants (sample counts, variance bounds, truncation behavior), not exact floating-point equality.

## Commit & Pull Request Guidelines
Git history is not included in this workspace snapshot, so no project-specific commit pattern can be inferred here. Use Conventional Commit style by default (for example, `feat: add r_max guard in trajectory tracing`, `fix: tighten plane-hit tolerance check`).

PRs should include:
- concise problem/solution summary
- linked issue (if available)
- test evidence (`pytest -q` output)
- notes on numerical impact and, for visualization changes, a screenshot or plot.
