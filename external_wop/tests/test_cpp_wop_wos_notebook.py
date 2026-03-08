from __future__ import annotations

import json
from pathlib import Path


def test_cpp_only_wop_wos_notebook_has_required_experiment_setup() -> None:
    workspace_root = Path(__file__).resolve().parents[2]
    candidates = [
        workspace_root / "wop-vs-wos-cpp-unit-cube.ipynb",
        workspace_root / "output" / "jupyter-notebook" / "wop-vs-wos-cpp-unit-cube.ipynb",
    ]
    notebook_path = next((path for path in candidates if path.exists()), None)
    assert notebook_path is not None, f"Notebook not found in any expected location: {candidates}"

    nb = json.loads(notebook_path.read_text(encoding="utf-8"))
    code_text = "\n".join("".join(cell.get("source", [])) for cell in nb.get("cells", []) if cell.get("cell_type") == "code")

    assert "N_VALUES = [100, 1000, 10000, 100000, 1000000]" in code_text
    assert 'METHODS = ["wop", "wos"]' in code_text

    assert "from wop.wop import estimate_wop" not in code_text
    assert "from wop.wos_box import estimate_wos_box" not in code_text
