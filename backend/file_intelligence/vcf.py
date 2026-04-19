"""Upload-time profiler for VCF / BCF variant files."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


def profile_vcf(path: str | Path) -> Dict[str, Any]:
    path = Path(path)

    try:
        return _profile_with_pandas(path)
    except Exception as exc:
        return {
            "format": "vcf",
            "family": "variant",
            "n_records": None,
            "summary": {},
            "schema": {},
            "sample": [],
            "available_sheets": None,
            "raw_metadata": {},
            "profiler_error": str(exc),
        }


def _profile_with_pandas(path: Path) -> Dict[str, Any]:
    """Parse VCF using pandas (no pysam dependency required)."""
    import gzip
    import pandas as pd

    suffix = path.suffix.lower()
    is_gz = suffix == ".gz" or str(path).endswith(".vcf.gz") or str(path).endswith(".bcf.gz")
    opener = gzip.open if is_gz else open

    header_lines: list[str] = []
    col_header: list[str] = []
    data_lines: list[str] = []
    max_data_preview = 200  # parse up to 200 rows for stats

    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line.rstrip("\n"))
            elif line.startswith("#CHROM"):
                col_header = line.lstrip("#").rstrip("\n").split("\t")
            elif len(data_lines) < max_data_preview:
                data_lines.append(line.rstrip("\n"))

    n_meta = len(header_lines)
    n_samples = max(0, len(col_header) - 9) if col_header else 0

    # Parse preview rows into a DataFrame
    df: Any = None
    if col_header and data_lines:
        import io
        text = "\t".join(col_header) + "\n" + "\n".join(data_lines)
        df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)

    chrom_dist: Dict[str, int] = {}
    filter_dist: Dict[str, int] = {}
    if df is not None and not df.empty:
        if "CHROM" in df.columns:
            chrom_dist = df["CHROM"].value_counts().head(20).to_dict()
        if "FILTER" in df.columns:
            filter_dist = df["FILTER"].value_counts().head(10).to_dict()

    # Count total variants (full scan for small files; estimate for large)
    n_variants = len(data_lines)  # at least what we parsed
    if len(data_lines) < max_data_preview:
        pass  # already got all of them
    else:
        # Continue counting without storing
        with opener(path, "rt", errors="replace") as fh:
            total = 0
            for line in fh:
                if not line.startswith("#"):
                    total += 1
        n_variants = total

    sample_rows = []
    if df is not None:
        for row in df.head(5).to_dict(orient="records"):
            sample_rows.append({k: v for k, v in row.items() if k in
                                 ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]})

    # Extract INFO field keys from header
    info_keys = [
        line.split("ID=")[1].split(",")[0]
        for line in header_lines
        if line.startswith("##INFO=<ID=")
    ]

    return {
        "format": "vcf",
        "family": "variant",
        "n_records": n_variants,
        "summary": {
            "n_variants": n_variants,
            "n_samples": n_samples,
            "n_header_lines": n_meta,
            "chrom_distribution": chrom_dist,
            "filter_distribution": filter_dist,
            "compressed": is_gz,
        },
        "schema": {
            "columns": col_header,
            "sample_names": col_header[9:] if len(col_header) > 9 else [],
            "info_fields": info_keys[:30],
        },
        "sample": sample_rows,
        "available_sheets": None,
        "raw_metadata": {"meta_lines": header_lines[:10]},
        "profiler_error": None,
    }
