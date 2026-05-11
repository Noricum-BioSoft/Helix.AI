"""
Sandboxed Python/pandas executor for LLM-generated data analysis code.

The execution namespace is intentionally minimal:
  - ``df``   — the user's DataFrame (a defensive copy)
  - ``pd``   — pandas
  - ``np``   — numpy
  - ``plt``  — matplotlib.pyplot (figures captured, not shown)
  - ``__builtins__`` — restricted to safe built-ins only

Generated code must store its final answer in a variable called ``result``.
"""
from __future__ import annotations

import ast
import builtins
import logging
import threading
from typing import Any, Dict

logger = logging.getLogger(__name__)


class _ImportAndCallRemover(ast.NodeTransformer):
    """AST transformer that removes import statements AND direct __import__ calls."""

    removed: int = 0

    def visit_Import(self, node: ast.Import) -> None:  # type: ignore[override]
        self.removed += 1
        return None  # remove from AST

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:  # type: ignore[override]
        self.removed += 1
        return None  # remove from AST

    def visit_Call(self, node: ast.Call) -> ast.AST:
        """Replace ``__import__(...)`` calls with None so the surrounding
        expression stays syntactically valid but produces no side-effect."""
        if (
            isinstance(node.func, ast.Name) and node.func.id == "__import__"
        ) or (
            isinstance(node.func, ast.Attribute) and node.func.attr == "__import__"
        ):
            self.removed += 1
            return ast.Constant(value=None)  # harmless replacement
        return self.generic_visit(node)


def _compile_without_imports(code: str) -> tuple[Any, int]:
    """Parse *code*, strip all import statements and ``__import__`` calls at
    the AST level, and compile the result to a code object.

    Returns ``(code_object, n_removed)``.  Falls back to compiling the raw
    string when the code has a ``SyntaxError`` (in that case n_removed = 0
    and execution will fail immediately with the real SyntaxError, which is
    more useful than a misleading ImportError).
    """
    try:
        tree = ast.parse(code, "<tabular_qa>", "exec")
    except SyntaxError:
        return compile(code, "<tabular_qa>", "exec"), 0

    remover = _ImportAndCallRemover()
    new_tree = remover.visit(tree)
    ast.fix_missing_locations(new_tree)

    if remover.removed:
        logger.warning(
            "[executor] stripped %d import/call node(s) from generated code before execution",
            remover.removed,
        )

    return compile(new_tree, "<tabular_qa>", "exec"), remover.removed

# Built-ins that are safe inside the sandbox
_SAFE_BUILTIN_NAMES = {
    "abs", "all", "any", "bool", "dict", "dir", "divmod", "enumerate",
    "filter", "float", "format", "frozenset", "getattr", "hasattr",
    "hash", "int", "isinstance", "issubclass", "iter", "len", "list",
    "map", "max", "min", "next", "print", "range", "repr", "reversed",
    "round", "set", "slice", "sorted", "str", "sum", "tuple", "type",
    "zip", "None", "True", "False",
}

_SAFE_BUILTINS = {k: getattr(builtins, k) for k in _SAFE_BUILTIN_NAMES if hasattr(builtins, k)}

# Maximum bytes of text/repr output returned to the LLM
_MAX_REPR_BYTES = 16_000


class ExecutionError(Exception):
    """Raised when sandboxed code raises an exception."""


class TimeoutError(Exception):
    """Raised when sandboxed execution exceeds the time limit."""


def _serialise_result(result: Any) -> Dict[str, Any]:
    """Convert the ``result`` variable to a JSON-friendly dict."""
    import pandas as pd
    import numpy as np

    if result is None:
        return {"type": "null", "value": None}

    if isinstance(result, pd.DataFrame):
        # Truncate large tables
        if len(result) > 500:
            truncated = True
            result = result.head(500)
        else:
            truncated = False
        return {
            "type": "dataframe",
            "value": result.to_dict(orient="records"),
            "columns": list(result.columns),
            "n_rows": len(result),
            "truncated": truncated,
        }

    if isinstance(result, pd.Series):
        return {
            "type": "series",
            "name": result.name,
            "value": result.head(200).to_dict(),
        }

    if isinstance(result, (np.integer, np.floating)):
        return {"type": "scalar", "value": result.item()}

    if isinstance(result, (int, float, bool, str)):
        return {"type": "scalar", "value": result}

    if isinstance(result, dict):
        return {"type": "dict", "value": result}

    if isinstance(result, (list, tuple)):
        return {"type": "list", "value": list(result)[:200]}

    # Fallback: repr (truncated)
    raw = repr(result)
    if len(raw) > _MAX_REPR_BYTES:
        raw = raw[:_MAX_REPR_BYTES] + "…"
    return {"type": "repr", "value": raw}


def execute_code(
    code: str,
    df: Any,  # pandas DataFrame
    *,
    timeout_s: int = 30,
) -> Dict[str, Any]:
    """
    Execute *code* in a restricted namespace containing *df*.

    Returns
    -------
    dict with keys:
        success  bool
        result   serialised result dict (if success)
        error    str (if not success)
        stdout   str  (captured print output)
        code     str  (the code that was executed)
    """
    import pandas as pd
    import numpy as np

    # Matplotlib: non-interactive backend so plt.show() is a no-op
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    captured_output: list[str] = []

    # Redirect print() inside sandbox.
    # Must be defined before the namespace so the closure captures the local list.
    def _safe_print(*args: Any, **kwargs: Any) -> None:
        sep = kwargs.get("sep", " ")
        captured_output.append(sep.join(str(a) for a in args))

    # Make a fresh copy of the safe-builtins dict for every call.
    # We must NOT mutate the module-level _SAFE_BUILTINS singleton: that would
    # permanently replace the real `print` with a closure over a previous
    # call's captured_output list, losing output on subsequent calls and
    # creating a race condition under concurrent execution.
    call_builtins = {**_SAFE_BUILTINS, "print": _safe_print}

    # Pre-load scientific libraries so LLM-generated code can use them
    # without needing import statements (which the sandbox blocks).
    try:
        import scipy.stats as _scipy_stats
        import scipy.cluster.hierarchy as _scipy_hierarchy
        import scipy.spatial.distance as _scipy_distance
        import scipy as _scipy
        _scipy.stats = _scipy_stats
    except ImportError:
        _scipy_stats = None
        _scipy = None

    try:
        import seaborn as _sns
    except ImportError:
        _sns = None

    namespace: Dict[str, Any] = {
        "__builtins__": call_builtins,
        "df": df.copy(),
        "pd": pd,
        "np": np,
        "plt": plt,
        **({"scipy": _scipy, "stats": _scipy_stats} if _scipy_stats else {}),
        **({"sns": _sns} if _sns else {}),
    }

    result_holder: Dict[str, Any] = {}
    exc_holder: list[Exception] = []

    def _run() -> None:
        try:
            code_obj, _ = _compile_without_imports(code)
            exec(code_obj, namespace)  # noqa: S102
            result_holder["value"] = namespace.get("result")
        except Exception as e:  # noqa: BLE001
            exc_holder.append(e)

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()
    thread.join(timeout=timeout_s)

    stdout = "\n".join(captured_output)

    if thread.is_alive():
        return {
            "success": False,
            "error": f"Execution timed out after {timeout_s}s.",
            "stdout": stdout,
            "code": code,
        }

    if exc_holder:
        raw_err = f"{type(exc_holder[0]).__name__}: {exc_holder[0]}"
        # `import foo` inside the sandbox causes `KeyError: '__import__'` because
        # __import__ is absent from the restricted builtins dict.  That message
        # is opaque to an LLM retry — translate it into an actionable instruction.
        if isinstance(exc_holder[0], KeyError) and "__import__" in str(exc_holder[0]):
            raw_err = (
                "SandboxImportError: import statements are not allowed in this sandbox. "
                "Do NOT write any import or from...import lines. "
                "The following names are already pre-bound and ready to use: "
                "df, pd, np, plt, sns, stats (scipy.stats), scipy. "
                "Remove all import statements from your code and use these names directly."
            )
        return {
            "success": False,
            "error": raw_err,
            "stdout": stdout,
            "code": code,
        }

    return {
        "success": True,
        "result": _serialise_result(result_holder.get("value")),
        "stdout": stdout,
        "code": code,
    }
