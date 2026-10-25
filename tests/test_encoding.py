"""Text files are read and written as UTF-8 everywhere: figures, pages and layout files carry α, β, η, ′ and −,
which the Windows default code page (cp1252) cannot encode."""

import ast
from pathlib import Path

ROOT = Path(__file__).parents[1]


def _unspecified(path: Path) -> list[str]:
    out = []
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if not isinstance(node, ast.Call) or any(k.arg == "encoding" for k in node.keywords):
            continue
        f = node.func
        if isinstance(f, ast.Attribute) and f.attr in ("read_text", "write_text"):
            out.append(f"{path.relative_to(ROOT)}:{node.lineno}")
        elif isinstance(f, ast.Name) and f.id == "open":
            mode = (
                node.args[1] if len(node.args) > 1 else next((k.value for k in node.keywords if k.arg == "mode"), None)
            )
            if not (isinstance(mode, ast.Constant) and "b" in str(mode.value)):
                out.append(f"{path.relative_to(ROOT)}:{node.lineno}")
    return out


def test_library_and_scripts_name_the_text_encoding():
    files = sorted([*(ROOT / "src").rglob("*.py"), *(ROOT / "scripts").rglob("*.py")])
    missing = [hit for f in files for hit in _unspecified(f)]
    assert not missing, "add encoding='utf-8' at " + ", ".join(missing)
