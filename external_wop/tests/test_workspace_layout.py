from __future__ import annotations

import json
from pathlib import Path


def test_cpp_project_is_separated_to_workspace_root() -> None:
    workspace_root = Path(__file__).resolve().parents[2]
    cpp_root = workspace_root / "external_wop_cpp"

    assert cpp_root.exists()
    assert cpp_root.is_dir()
    assert (cpp_root / "CMakeLists.txt").exists()


def test_comparison_notebook_exists_for_cpp_python_runs() -> None:
    workspace_root = Path(__file__).resolve().parents[2]
    notebook_path = workspace_root / "compare_wop_cpp_python.ipynb"

    assert notebook_path.exists()

    nb = json.loads(notebook_path.read_text(encoding="utf-8"))
    cells = nb.get("cells", [])
    source_text = "\n".join(
        "".join(cell.get("source", []))
        for cell in cells
        if isinstance(cell, dict)
    )
    assert "external_wop_cpp" in source_text
    assert "external_wop" in source_text
    assert "wop_cli" in source_text


def test_comparison_notebook_has_cmake_discovery_logic() -> None:
    workspace_root = Path(__file__).resolve().parents[2]
    notebook_path = workspace_root / "compare_wop_cpp_python.ipynb"

    nb = json.loads(notebook_path.read_text(encoding="utf-8"))
    cells = nb.get("cells", [])
    source_text = "\n".join(
        "".join(cell.get("source", []))
        for cell in cells
        if isinstance(cell, dict)
    )
    assert "shutil.which(\"cmake\")" in source_text
    assert "CMake is not available in PATH" in source_text
