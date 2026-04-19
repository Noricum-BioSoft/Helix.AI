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

import builtins
import threading
from typing import Any, Dict

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

    namespace: Dict[str, Any] = {
        "__builtins__": call_builtins,
        "df": df.copy(),
        "pd": pd,
        "np": np,
        "plt": plt,
    }

    result_holder: Dict[str, Any] = {}
    exc_holder: list[Exception] = []

    def _run() -> None:
        try:
            exec(compile(code, "<tabular_qa>", "exec"), namespace)  # noqa: S102
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
        return {
            "success": False,
            "error": f"{type(exc_holder[0]).__name__}: {exc_holder[0]}",
            "stdout": stdout,
            "code": code,
        }

    return {
        "success": True,
        "result": _serialise_result(result_holder.get("value")),
        "stdout": stdout,
        "code": code,
    }
