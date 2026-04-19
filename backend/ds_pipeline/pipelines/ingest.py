"""
Data ingestion: load tabular data (CSV / TSV / Excel) into DuckDB + pandas.

Supported formats
-----------------
* .csv  — comma-separated
* .tsv  — tab-separated
* .xlsx / .xls — Excel workbook (sheet selection supported)
* .txt  — treated as CSV; delimiter auto-detected

The primary entry point is ``ingest_tabular``.  ``ingest_csv`` is kept as a
backwards-compatibility shim that delegates to ``ingest_tabular``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

TABULAR_EXTENSIONS = {".csv", ".tsv", ".xlsx", ".xls", ".txt"}


def list_sheets(path: str | Path) -> List[str]:
    """Return sheet names for an Excel workbook, or ``["default"]`` for flat files."""
    import pandas as pd

    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        xl = pd.ExcelFile(path)
        return xl.sheet_names
    return ["default"]


def ingest_tabular(
    path: str | Path,
    *,
    sheet: Optional[str] = None,
    delimiter: Optional[str] = None,
) -> Tuple[Any, Any, Dict[str, Any]]:
    """
    Load tabular data from *path* into DuckDB + pandas.

    Parameters
    ----------
    path :
        Filesystem path to the data file.  Must exist.
    sheet :
        For Excel workbooks, the sheet name or 0-based index to load.
        Defaults to the first sheet when ``None``.
    delimiter :
        Override the column separator for CSV/TSV/TXT files.
        Inferred from extension when ``None`` (``,`` for .csv/.txt, ``\\t`` for .tsv).

    Returns
    -------
    conn : duckdb.DuckDBPyConnection
        In-memory DuckDB connection with the DataFrame registered as ``"data"``.
    df : pandas.DataFrame
        The loaded dataset.
    metadata : dict
        Shape, column types, nullability, 5-row sample, source sheet/format info.
    """
    import duckdb
    import pandas as pd

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")

    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xls"}:
        xl = pd.ExcelFile(path)
        available_sheets = xl.sheet_names
        if sheet is None:
            sheet_to_load = available_sheets[0]
        elif isinstance(sheet, int):
            sheet_to_load = available_sheets[sheet]
        elif sheet in available_sheets:
            sheet_to_load = sheet
        else:
            raise ValueError(
                f"Sheet '{sheet}' not found in '{path.name}'. "
                f"Available sheets: {available_sheets}"
            )
        df = pd.read_excel(path, sheet_name=sheet_to_load)
        source_format = "excel"
        source_sheet = sheet_to_load
    else:
        if delimiter is None:
            delimiter = "\t" if suffix == ".tsv" else ","
        df = pd.read_csv(path, sep=delimiter)
        source_format = "tsv" if suffix == ".tsv" else "csv"
        source_sheet = None

    conn = duckdb.connect(database=":memory:")
    conn.register("data", df)

    col_info = []
    for col in df.columns:
        col_info.append({
            "name": col,
            "dtype": str(df[col].dtype),
            "n_missing": int(df[col].isnull().sum()),
            "pct_missing": round(float(df[col].isnull().mean()), 4),
            "n_unique": int(df[col].nunique()),
        })

    metadata: Dict[str, Any] = {
        "path": str(path),
        "source_format": source_format,
        "source_sheet": source_sheet,
        "available_sheets": list_sheets(path) if source_format == "excel" else None,
        "n_rows": len(df),
        "n_cols": len(df.columns),
        "columns": col_info,
        "sample": df.head(5).to_dict(orient="records"),
        "size_bytes": path.stat().st_size,
    }
    return conn, df, metadata


def ingest_csv(path: str | Path) -> Tuple[Any, Any, Dict[str, Any]]:
    """Backwards-compatible shim — delegates to :func:`ingest_tabular`."""
    return ingest_tabular(path)
