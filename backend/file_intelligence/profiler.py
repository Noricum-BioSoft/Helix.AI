"""
Dispatch upload-time profiling by file format.

Each profiler returns a ``FileProfile`` dict with a consistent top-level shape:
    format          str   — canonical format name
    family          str   — capability family (tabular / sequence / alignment / …)
    n_records       int   — rows / reads / variants / cells (format-dependent)
    summary         dict  — format-specific top-level stats
    schema          dict  — column/field definitions where applicable
    sample          list  — first few records as plain dicts / strings
    available_sheets list — Excel sheet names (tabular only)
    raw_metadata    dict  — any additional format-specific metadata
    profiler_error  str   — set only when profiling partially failed (non-fatal)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional


# ---------------------------------------------------------------------------
# Extension → (family, profiler_fn) registry
# ---------------------------------------------------------------------------

SUPPORTED_EXTENSIONS: Dict[str, str] = {
    # Tabular
    ".csv":  "tabular",
    ".tsv":  "tabular",
    ".xlsx": "tabular",
    ".xls":  "tabular",
    ".txt":  "tabular",
    # Sequence
    ".fasta": "sequence",
    ".fa":    "sequence",
    ".fas":   "sequence",
    ".fastq": "sequence",
    ".fq":    "sequence",
    # Variants
    ".vcf":  "variant",
    ".bcf":  "variant",
    # Genomic intervals / annotation
    ".bed":  "genomic_interval",
    ".gff":  "genomic_interval",
    ".gff3": "genomic_interval",
    ".gtf":  "genomic_interval",
    # Single-cell
    ".h5ad": "single_cell",
    ".loom": "single_cell",
    ".h5":   "single_cell",
    # Alignment
    ".sam":  "alignment",
    ".bam":  "alignment",
    ".cram": "alignment",
    # Generic compressed (resolve by inner extension)
    ".gz":   "compressed",
}


def _resolve_gz_inner_suffix(path: Path) -> str:
    """Return the inner extension of a .gz file (e.g. '.fastq' for 'reads.fastq.gz')."""
    stem = path.stem  # e.g. 'reads.fastq' from 'reads.fastq.gz'
    return Path(stem).suffix.lower()  # '.fastq'


# Families whose profilers natively read gzip-compressed files.
_GZ_NATIVE_FAMILIES = {"sequence"}


def profile_file(path: str | Path, *, sheet: Optional[str] = None) -> Dict[str, Any]:
    """
    Profile *path* at upload time and return a ``FileProfile`` dict.

    Parameters
    ----------
    path :
        Filesystem path to the uploaded file.
    sheet :
        For multi-sheet files (Excel), the sheet to profile first.
        ``None`` → profile the first sheet.

    Returns
    -------
    dict  Always present keys: ``format``, ``family``, ``n_records``, ``summary``,
          ``schema``, ``sample``.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    # For .gz files, use the inner extension only for family dispatch.
    # The original compressed path is passed to profilers that handle .gz
    # natively (sequence).  For all other families, .gz is not yet supported
    # because their underlying parsers (pandas, pysam, h5py…) require the
    # decompressed file; we return a clear error rather than silently failing.
    if suffix == ".gz":
        inner_suffix = _resolve_gz_inner_suffix(path)
        family = SUPPORTED_EXTENSIONS.get(inner_suffix, "unknown")
        if family not in _GZ_NATIVE_FAMILIES:
            return {
                "format": f"{inner_suffix.lstrip('.')}.gz",
                "family": family if family != "unknown" else "unknown",
                "n_records": None,
                "summary": {},
                "schema": {},
                "sample": [],
                "available_sheets": None,
                "raw_metadata": {},
                "profiler_error": (
                    f"gzip-compressed {inner_suffix.lstrip('.')} files are not yet "
                    "supported for upload-time profiling. Please decompress the file "
                    "before uploading."
                ),
            }
        dispatch_suffix = inner_suffix
    else:
        dispatch_suffix = suffix
    family = SUPPORTED_EXTENSIONS.get(dispatch_suffix, "unknown")

    if family == "tabular":
        from backend.file_intelligence.tabular import profile_tabular
        return profile_tabular(path, sheet=sheet)

    if family == "sequence":
        from backend.file_intelligence.sequence import profile_sequence
        return profile_sequence(path)

    if family == "variant":
        from backend.file_intelligence.vcf import profile_vcf
        return profile_vcf(path)

    if family == "genomic_interval":
        from backend.file_intelligence.bed import profile_bed
        return profile_bed(path)

    if family == "single_cell":
        from backend.file_intelligence.singlecell import profile_singlecell
        return profile_singlecell(path)

    if family == "alignment":
        from backend.file_intelligence.alignment import profile_alignment
        return profile_alignment(path)

    # Fallback for unknown / compressed-only
    return {
        "format": suffix.lstrip(".") or "unknown",
        "family": family,
        "n_records": None,
        "summary": {},
        "schema": {},
        "sample": [],
        "available_sheets": None,
        "raw_metadata": {},
        "profiler_error": f"No profiler registered for extension '{suffix}'.",
    }
