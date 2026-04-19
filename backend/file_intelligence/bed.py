"""Upload-time profiler for BED, GFF, GTF, GFF3 genomic interval files."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


_BED_STANDARD_COLS = ["chrom", "chromStart", "chromEnd", "name",
                      "score", "strand", "thickStart", "thickEnd",
                      "itemRgb", "blockCount", "blockSizes", "blockStarts"]

_GFF_COLS = ["seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes"]


def profile_bed(path: str | Path) -> Dict[str, Any]:
    path = Path(path)
    suffix = path.suffix.lower()
    fmt = _detect_format(path)

    try:
        if fmt in {"gff", "gff3", "gtf"}:
            return _profile_gff(path, fmt)
        return _profile_bed(path, fmt)
    except Exception as exc:
        return {
            "format": fmt,
            "family": "genomic_interval",
            "n_records": None,
            "summary": {},
            "schema": {},
            "sample": [],
            "available_sheets": None,
            "raw_metadata": {},
            "profiler_error": str(exc),
        }


def _detect_format(path: Path) -> str:
    suffix = path.suffix.lower()
    name = path.name.lower()
    if suffix in {".gff", ".gff3"} or ".gff" in name:
        return "gff3"
    if suffix == ".gtf" or ".gtf" in name:
        return "gtf"
    if "narrowpeak" in name or "broadpeak" in name:
        return "peak"
    return "bed"


def _profile_bed(path: Path, fmt: str) -> Dict[str, Any]:
    import pandas as pd

    # Read skipping track/browser header lines
    rows: list[str] = []
    with open(path, errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(("track", "browser", "#")):
                continue
            rows.append(s)

    if not rows:
        raise ValueError("BED file contains no data rows.")

    n_cols = len(rows[0].split("\t"))
    col_names = _BED_STANDARD_COLS[:n_cols] if n_cols <= 12 else [f"col{i}" for i in range(n_cols)]

    import io
    df = pd.read_csv(io.StringIO("\n".join(rows)), sep="\t", header=None,
                     names=col_names, dtype=str)

    chrom_dist: Dict[str, int] = {}
    if "chrom" in df.columns:
        chrom_dist = df["chrom"].value_counts().head(20).to_dict()

    interval_widths: list[int] = []
    if "chromStart" in df.columns and "chromEnd" in df.columns:
        df["_width"] = pd.to_numeric(df["chromEnd"], errors="coerce") - \
                       pd.to_numeric(df["chromStart"], errors="coerce")
        widths = df["_width"].dropna()
        if not widths.empty:
            interval_widths = [int(widths.min()), int(widths.max()), round(float(widths.mean()), 1)]

    strand_dist: Dict[str, int] = {}
    if "strand" in df.columns:
        strand_dist = df["strand"].value_counts().to_dict()

    sample = df.drop(columns=["_width"], errors="ignore").head(5).to_dict(orient="records")

    return {
        "format": fmt,
        "family": "genomic_interval",
        "n_records": len(df),
        "summary": {
            "n_intervals": len(df),
            "n_cols": n_cols,
            "chrom_distribution": chrom_dist,
            "strand_distribution": strand_dist,
            "interval_width_min_max_mean": interval_widths or None,
        },
        "schema": {"columns": col_names},
        "sample": sample,
        "available_sheets": None,
        "raw_metadata": {},
        "profiler_error": None,
    }


def _profile_gff(path: Path, fmt: str) -> Dict[str, Any]:
    import pandas as pd
    import io

    rows: list[str] = []
    directive_lines: list[str] = []
    with open(path, errors="replace") as fh:
        for line in fh:
            s = line.rstrip("\n")
            if s.startswith("##"):
                directive_lines.append(s)
            elif not s.startswith("#") and s.strip():
                rows.append(s)

    if not rows:
        raise ValueError(f"{fmt.upper()} file contains no data rows.")

    df = pd.read_csv(io.StringIO("\n".join(rows)), sep="\t", header=None,
                     names=_GFF_COLS, dtype=str)

    feature_dist: Dict[str, int] = {}
    if "feature" in df.columns:
        feature_dist = df["feature"].value_counts().head(15).to_dict()

    seqname_dist: Dict[str, int] = {}
    if "seqname" in df.columns:
        seqname_dist = df["seqname"].value_counts().head(20).to_dict()

    return {
        "format": fmt,
        "family": "genomic_interval",
        "n_records": len(df),
        "summary": {
            "n_features": len(df),
            "feature_types": feature_dist,
            "seqname_distribution": seqname_dist,
            "n_directives": len(directive_lines),
        },
        "schema": {"columns": _GFF_COLS},
        "sample": df.head(5).to_dict(orient="records"),
        "available_sheets": None,
        "raw_metadata": {"directives": directive_lines[:5]},
        "profiler_error": None,
    }
