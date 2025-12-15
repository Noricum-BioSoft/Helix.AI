"""
Tool inventory: a single source of truth for "what tools do you have?"

This module builds a deterministic inventory from:
- Central MCP/tool schema registry (backend/tool_schemas.py)
- Static discovery of @tool-decorated python functions (AST scan; no imports)
- Environment capabilities (local CLI tools discovered via `shutil.which`)
- Expected EC2 toolchain (from scripts/aws setup + EC2 user-data in ec2_executor.py)
"""

from __future__ import annotations

import ast
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from backend.tool_schemas import list_tool_schemas


# Keep in sync with:
# - backend/ec2_executor.py user-data conda install line
# - scripts/aws/setup-bioinformatics-ec2.sh
# - scripts/aws/complete-ec2-setup.sh
EC2_EXPECTED_CLI_TOOLS: List[str] = [
    "bbmerge.sh",
    "flash",
    "pear",
    "samtools",
    "bcftools",
    "bwa",
    "bowtie2",
    "fastqc",
]


DEFAULT_LOCAL_CLI_CANDIDATES: List[str] = [
    # Bioinformatics tools
    *EC2_EXPECTED_CLI_TOOLS,
    # Common runtime/tooling
    "Rscript",
    "python3",
    "conda",
    "aws",
]


@dataclass(frozen=True)
class DecoratedTool:
    name: str
    file: str
    doc: str


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _safe_docstring(node: ast.AST) -> str:
    try:
        return (ast.get_docstring(node) or "").strip()
    except Exception:
        return ""


def _is_tool_decorator(dec: ast.AST) -> bool:
    # Accept:
    # - @tool
    # - @tool(...)
    if isinstance(dec, ast.Name) and dec.id == "tool":
        return True
    if isinstance(dec, ast.Call) and isinstance(dec.func, ast.Name) and dec.func.id == "tool":
        return True
    return False


def discover_decorated_tools(
    search_dirs: Optional[Sequence[Path]] = None,
) -> List[DecoratedTool]:
    """
    Find python functions decorated with @tool without importing modules.
    This keeps discovery safe (no side effects from module imports).
    """
    root = _repo_root()
    if search_dirs is None:
        search_dirs = [
            root / "backend",
            root / "tools",
        ]

    out: List[DecoratedTool] = []
    for d in search_dirs:
        if not d.exists():
            continue
        for py_file in d.rglob("*.py"):
            # Skip huge virtualenv-like folders if any ever get checked in
            if any(part in ("node_modules", ".venv", "venv", "__pycache__") for part in py_file.parts):
                continue
            try:
                src = py_file.read_text(encoding="utf-8")
            except Exception:
                continue
            try:
                tree = ast.parse(src)
            except SyntaxError:
                continue

            for node in tree.body:
                if isinstance(node, ast.FunctionDef) and node.decorator_list:
                    if any(_is_tool_decorator(dec) for dec in node.decorator_list):
                        out.append(
                            DecoratedTool(
                                name=node.name,
                                file=str(py_file.relative_to(root)),
                                doc=_safe_docstring(node),
                            )
                        )
    # Deterministic ordering
    out.sort(key=lambda t: (t.name, t.file))
    return out


def detect_cli_tools(candidates: Sequence[str]) -> List[Dict[str, str]]:
    found: List[Dict[str, str]] = []
    for name in candidates:
        name = name.strip()
        if not name:
            continue
        path = shutil.which(name)
        if path:
            found.append({"name": name, "path": path})
    # Deduplicate by name deterministically
    uniq: Dict[str, Dict[str, str]] = {}
    for item in found:
        uniq[item["name"]] = item
    return [uniq[k] for k in sorted(uniq.keys())]


def build_toolbox_inventory() -> Dict[str, Any]:
    """
    Main inventory payload used for:
    - answering "what tools do you have?"
    - diagnostics / debugging tool availability
    """
    schemas = list_tool_schemas()
    schema_names = sorted({s.get("name", "") for s in schemas if isinstance(s, dict) and s.get("name")})

    decorated = discover_decorated_tools()

    env_candidates_raw = os.getenv("HELIX_CLI_TOOL_CANDIDATES", "")
    if env_candidates_raw.strip():
        candidates = [p.strip() for p in env_candidates_raw.split(",") if p.strip()]
    else:
        candidates = DEFAULT_LOCAL_CLI_CANDIDATES

    local_cli = detect_cli_tools(candidates)

    use_ec2 = os.getenv("HELIX_USE_EC2", "false").lower().strip() == "true"
    ec2_configured = bool(os.getenv("HELIX_EC2_INSTANCE_ID") and os.getenv("HELIX_EC2_KEY_FILE"))

    return {
        "tool_schemas": {
            "count": len(schemas),
            "tools": schema_names,
        },
        "python_decorated_tools": {
            "count": len(decorated),
            "tools": [
                {"name": t.name, "file": t.file, "description": t.doc}
                for t in decorated
            ],
        },
        "local_cli_tools": {
            "candidates": candidates,
            "found": local_cli,
            "found_count": len(local_cli),
        },
        "ec2": {
            "enabled": use_ec2,
            "configured": ec2_configured,
            "expected_cli_tools": EC2_EXPECTED_CLI_TOOLS,
            "note": (
                "EC2-installed tools are listed as 'expected'. "
                "To verify actual EC2 availability at runtime, the inventory would need SSH access "
                "(HELIX_EC2_INSTANCE_ID + HELIX_EC2_KEY_FILE)."
            ),
        },
        "dynamic_tool_generation": {
            "available": True,
            "implementation": "backend/tool_generator_agent.py",
            "note": "When no pre-existing tool exists, Helix can generate and run a one-off tool (optionally on EC2/EMR).",
        },
    }


def format_toolbox_inventory_markdown(inv: Dict[str, Any]) -> str:
    """
    Human-readable markdown summary suitable for chat answers.
    """
    schemas = inv.get("tool_schemas", {})
    schema_tools = schemas.get("tools", [])

    decorated = inv.get("python_decorated_tools", {})
    decorated_tools = decorated.get("tools", [])

    local_cli = inv.get("local_cli_tools", {})
    local_found = local_cli.get("found", [])

    ec2 = inv.get("ec2", {})

    lines: List[str] = []
    lines.append("### Toolbox inventory")
    lines.append("")

    lines.append(f"- **Registered Helix tools (MCP schemas)**: {schemas.get('count', 0)}")
    for name in schema_tools:
        lines.append(f"  - `{name}`")
    lines.append("")

    lines.append(f"- **Python `@tool` functions discovered in repo**: {decorated.get('count', 0)}")
    for t in decorated_tools:
        nm = t.get("name", "")
        fp = t.get("file", "")
        lines.append(f"  - `{nm}` ({fp})")
    lines.append("")

    lines.append(f"- **Local CLI tools detected on this machine**: {local_cli.get('found_count', 0)}")
    for item in local_found:
        lines.append(f"  - `{item.get('name')}` → `{item.get('path')}`")
    lines.append("")

    lines.append("- **EC2 bioinformatics toolchain**:")
    lines.append(f"  - **HELIX_USE_EC2 enabled**: {bool(ec2.get('enabled'))}")
    lines.append(f"  - **EC2 configured (instance+key)**: {bool(ec2.get('configured'))}")
    lines.append("  - **Expected tools on the Helix EC2 image/setup**:")
    for name in ec2.get("expected_cli_tools", []):
        lines.append(f"    - `{name}`")
    lines.append("")

    lines.append("- **Dynamic tool generation**:")
    lines.append("  - If a capability isn’t in the toolbox, Helix can generate and execute a custom tool (via `backend/tool_generator_agent.py`).")
    return "\n".join(lines)


