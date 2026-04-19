"""Upload-time profiler for SAM / BAM / CRAM alignment files."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


def profile_alignment(path: str | Path) -> Dict[str, Any]:
    path = Path(path)
    suffix = path.suffix.lower()

    try:
        import pysam  # optional dependency

        mode = "rb" if suffix in {".bam", ".cram"} else "r"
        with pysam.AlignmentFile(str(path), mode, check_sq=False) as af:
            header = af.header.to_dict()
            references = list(af.references)[:20]
            lengths = list(af.lengths)[:20]

            # Read a sample of alignments
            sample_reads: list[Dict[str, Any]] = []
            n_reads = 0
            n_mapped = 0
            n_unmapped = 0
            for read in af.fetch(until_eof=True):
                n_reads += 1
                if read.is_unmapped:
                    n_unmapped += 1
                else:
                    n_mapped += 1
                if len(sample_reads) < 5:
                    sample_reads.append({
                        "query_name": read.query_name,
                        "reference_name": read.reference_name,
                        "reference_start": read.reference_start,
                        "mapping_quality": read.mapping_quality,
                        "query_length": read.query_length,
                        "is_paired": read.is_paired,
                    })
                if n_reads >= 10000:
                    break  # cap scan for large files

        return {
            "format": suffix.lstrip("."),
            "family": "alignment",
            "n_records": n_reads,
            "summary": {
                "n_reads_sampled": n_reads,
                "n_mapped": n_mapped,
                "n_unmapped": n_unmapped,
                "pct_mapped": round(100 * n_mapped / n_reads, 2) if n_reads else None,
                "references": references,
                "reference_lengths": lengths,
            },
            "schema": {"fields": ["query_name", "flag", "reference_name", "pos",
                                   "mapq", "cigar", "seq", "qual"]},
            "sample": sample_reads,
            "available_sheets": None,
            "raw_metadata": {
                "header_SQ_count": len(header.get("SQ", [])),
                "program_records": [p.get("ID") for p in header.get("PG", [])][:5],
            },
            "profiler_error": None,
        }

    except ImportError:
        return {
            "format": suffix.lstrip("."),
            "family": "alignment",
            "n_records": None,
            "summary": {"note": "pysam not installed; alignment profiling unavailable."},
            "schema": {},
            "sample": [],
            "available_sheets": None,
            "raw_metadata": {},
            "profiler_error": "pysam not available",
        }
    except Exception as exc:
        return {
            "format": suffix.lstrip("."),
            "family": "alignment",
            "n_records": None,
            "summary": {},
            "schema": {},
            "sample": [],
            "available_sheets": None,
            "raw_metadata": {},
            "profiler_error": str(exc),
        }
