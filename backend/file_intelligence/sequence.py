"""Upload-time profiler for FASTA / FASTQ files."""
from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Any, Dict


def profile_sequence(path: str | Path) -> Dict[str, Any]:
    path = Path(path)
    suffix = path.suffix.lower()
    is_gz = suffix == ".gz"
    inner = Path(path.stem).suffix.lower() if is_gz else suffix
    is_fastq = inner in {".fastq", ".fq"}

    try:
        opener = gzip.open if is_gz else open
        mode = "rt"

        n_records = 0
        seq_lengths: list[int] = []
        sample_records: list[Dict[str, Any]] = []
        sample_limit = 5

        with opener(path, mode, errors="replace") as fh:
            if is_fastq:
                while True:
                    header = fh.readline()
                    if not header:
                        break
                    seq = fh.readline().rstrip("\n")
                    fh.readline()  # +
                    qual = fh.readline().rstrip("\n")
                    n_records += 1
                    seq_lengths.append(len(seq))
                    if len(sample_records) < sample_limit:
                        sample_records.append({
                            "header": header.strip(),
                            "seq_len": len(seq),
                            "qual_min": min(ord(c) - 33 for c in qual) if qual else None,
                            "qual_max": max(ord(c) - 33 for c in qual) if qual else None,
                        })
            else:
                # FASTA
                current_id = None
                current_len = 0
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        if current_id is not None:
                            seq_lengths.append(current_len)
                            if len(sample_records) < sample_limit:
                                sample_records.append({"id": current_id, "seq_len": current_len})
                        current_id = line[1:].split()[0]
                        current_len = 0
                        n_records += 1
                    else:
                        current_len += len(re.sub(r"\s", "", line))
                if current_id is not None:
                    seq_lengths.append(current_len)
                    if len(sample_records) < sample_limit:
                        sample_records.append({"id": current_id, "seq_len": current_len})

        mean_len = sum(seq_lengths) / len(seq_lengths) if seq_lengths else 0
        min_len = min(seq_lengths) if seq_lengths else 0
        max_len = max(seq_lengths) if seq_lengths else 0

        return {
            "format": "fastq" if is_fastq else "fasta",
            "family": "sequence",
            "n_records": n_records,
            "summary": {
                "n_sequences": n_records,
                "mean_length": round(mean_len, 1),
                "min_length": min_len,
                "max_length": max_len,
                "compressed": is_gz,
            },
            "schema": {"fields": ["id", "sequence"] + (["quality"] if is_fastq else [])},
            "sample": sample_records,
            "available_sheets": None,
            "raw_metadata": {},
            "profiler_error": None,
        }

    except Exception as exc:
        return _error_profile("fastq" if is_fastq else "fasta", "sequence", exc)


def _error_profile(fmt: str, family: str, exc: Exception) -> Dict[str, Any]:
    return {
        "format": fmt, "family": family, "n_records": None,
        "summary": {}, "schema": {}, "sample": [],
        "available_sheets": None, "raw_metadata": {},
        "profiler_error": str(exc),
    }
