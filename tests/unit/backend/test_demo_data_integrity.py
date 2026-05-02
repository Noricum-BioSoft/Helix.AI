"""
test_demo_data_integrity.py

Fast, dependency-free checks that catch the two classes of demo breakage
we have repeatedly encountered:

  1. S3 data referenced in a followUpPrompt does not exist.
  2. A followUpPrompt contains accession numbers / non-FASTA text that
     would be misrouted or passed directly to a FASTA-parsing tool.
  3. A followUpPrompt is empty / None when a demo card shows "Load & run".

These run offline (boto3 S3 head_object calls are mocked) and are
intended to be part of the standard unit-test suite so they catch
regressions before any deployment step.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import List, Tuple
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Pull demo scenario data from the TypeScript source so there is a single
# source of truth — we parse it with a lightweight regex rather than running
# a full TS compiler.
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEMO_FILE = PROJECT_ROOT / "frontend" / "src" / "data" / "demoScenarios.ts"


def _extract_demo_fields(ts_source: str) -> List[dict]:
    """
    Extract (id, tool, expectedBehavior, followUpPrompt) from the TS file
    using regex so we don't need a JS runtime.
    """
    # Split on top-level object literals inside the demoScenarios array
    block_pattern = re.compile(
        r"\{\s*id:\s*'([^']+)'.*?(?=\{\s*id:|];)", re.DOTALL
    )
    demos = []
    for m in block_pattern.finditer(ts_source):
        block = m.group(0)
        demo: dict = {"id": m.group(1)}

        tool_m = re.search(r"tool:\s*'([^']+)'", block)
        demo["tool"] = tool_m.group(1) if tool_m else None

        behavior_m = re.search(r"expectedBehavior:\s*'([^']+)'", block)
        demo["expectedBehavior"] = behavior_m.group(1) if behavior_m else None

        # Extract backtick-delimited followUpPrompt (may be multiline)
        fup_m = re.search(r"followUpPrompt:\s*`(.*?)`", block, re.DOTALL)
        demo["followUpPrompt"] = fup_m.group(1).strip() if fup_m else None

        demos.append(demo)
    return demos


@pytest.fixture(scope="module")
def demos() -> List[dict]:
    src = DEMO_FILE.read_text(encoding="utf-8")
    result = _extract_demo_fields(src)
    assert result, "Could not parse any demo scenarios from demoScenarios.ts"
    return result


@pytest.fixture(scope="module")
def demos_with_followup(demos) -> List[dict]:
    return [d for d in demos if d.get("followUpPrompt")]


# ---------------------------------------------------------------------------
# 1. S3 data integrity
# ---------------------------------------------------------------------------

_S3_URI_RE = re.compile(r"s3://([^/\s]+)/([^\s,`]+)")


def _extract_s3_uris(text: str) -> List[Tuple[str, str]]:
    """Return list of (bucket, key) from all s3:// URIs in text."""
    return [(m.group(1), m.group(2)) for m in _S3_URI_RE.finditer(text)]


@pytest.mark.s3
class TestS3DataExists:
    """Every S3 URI embedded in a followUpPrompt must exist in S3."""

    def _check_object(self, bucket: str, key: str) -> None:
        """Raise AssertionError with a helpful message if the object is missing."""
        import boto3
        from botocore.exceptions import ClientError

        s3 = boto3.client("s3")
        try:
            s3.head_object(Bucket=bucket, Key=key)
        except ClientError as exc:
            code = exc.response["Error"]["Code"]
            pytest.fail(
                f"S3 object s3://{bucket}/{key} does not exist (HTTP {code}).\n"
                f"Fix: upload the demo data file or update the followUpPrompt."
            )

    def test_all_followup_s3_uris_exist(self, demos_with_followup):
        """All S3 URIs referenced in demo followUpPrompts must be reachable."""
        missing = []
        for demo in demos_with_followup:
            uris = _extract_s3_uris(demo["followUpPrompt"])
            for bucket, key in uris:
                try:
                    import boto3
                    from botocore.exceptions import ClientError
                    boto3.client("s3").head_object(Bucket=bucket, Key=key)
                except Exception as exc:
                    missing.append(
                        f"  [{demo['id']}] s3://{bucket}/{key} — {exc}"
                    )
        if missing:
            pytest.fail(
                "The following S3 data files referenced in demo followUpPrompts "
                "do not exist:\n" + "\n".join(missing)
            )


# ---------------------------------------------------------------------------
# 2. followUpPrompt format sanity
# ---------------------------------------------------------------------------

# Accession patterns that should NOT appear directly in followUpPrompts as
# the sole data source (they look like data but aren't fetchable by most tools)
_ACCESSION_ONLY_RE = re.compile(
    r"^\s*\w[\w\-]+:\s*[A-Z]{1,2}_?\d{5,9}\.\d\b",   # e.g.  Wuhan-Hu-1: MN908947.3
    re.MULTILINE,
)
_FASTA_RE = re.compile(r"^>", re.MULTILINE)
_S3_PATH_RE = re.compile(r"s3://")


class TestFollowUpPromptFormat:
    """Guard against the accession-numbers-without-S3 class of bug."""

    def test_no_followup_is_none_when_load_and_run_exists(self, demos):
        """Demos whose initial prompt produces needs_inputs must have a followUpPrompt."""
        for demo in demos:
            if demo.get("expectedBehavior") == "needs_inputs":
                assert demo.get("followUpPrompt"), (
                    f"Demo '{demo['id']}' has expectedBehavior='needs_inputs' but "
                    f"no followUpPrompt — the 'Load & run' button will do nothing."
                )

    def test_followup_prompts_not_empty(self, demos_with_followup):
        for demo in demos_with_followup:
            assert demo["followUpPrompt"].strip(), (
                f"Demo '{demo['id']}' has an empty followUpPrompt."
            )

    def test_followup_not_accessions_only(self, demos_with_followup):
        """
        A followUpPrompt must NOT consist solely of accession-number lines
        (e.g. 'Wuhan-Hu-1: MN908947.3') without an S3 URI or inline FASTA.

        This is the pattern that caused the phylogenetics 'No sequences provided'
        and 'At least 2 sequences required' errors.
        """
        for demo in demos_with_followup:
            prompt = demo["followUpPrompt"]
            has_s3 = bool(_S3_PATH_RE.search(prompt))
            has_fasta = bool(_FASTA_RE.search(prompt))
            accession_lines = _ACCESSION_ONLY_RE.findall(prompt)

            if accession_lines and not has_s3 and not has_fasta:
                pytest.fail(
                    f"Demo '{demo['id']}' followUpPrompt contains accession numbers "
                    f"({accession_lines[:2]}...) but no S3 URI or inline FASTA.\n"
                    f"Tools will receive these as sequences and fail with "
                    f"'No sequences provided' or 'At least 2 sequences required'.\n"
                    f"Fix: add 'sequences: s3://...' pointing to actual FASTA data, "
                    f"or embed inline FASTA in the prompt."
                )


# ---------------------------------------------------------------------------
# 3. Tool input resolution unit tests
# ---------------------------------------------------------------------------

class TestFastaInputResolution:
    """Unit tests for backend/phylogenetic_tree._resolve_fasta_input."""

    @pytest.fixture(autouse=True)
    def _import(self):
        sys.path.insert(0, str(PROJECT_ROOT))
        from backend.phylogenetic_tree import _resolve_fasta_input, SARS_COV2_SPIKE_FASTA
        self._resolve = _resolve_fasta_input
        self._builtin = SARS_COV2_SPIKE_FASTA

    def test_plain_fasta_returned_unchanged(self):
        fasta = ">seq1\nATCG\n>seq2\nGCTA"
        assert self._resolve(fasta) == fasta

    def test_empty_string_returned_unchanged(self):
        assert self._resolve("") == ""

    def test_s3_uri_fetches_content(self):
        fake_body = b">Wuhan\nMFVF\n>Alpha\nMFVF"
        mock_s3 = MagicMock()
        mock_s3.get_object.return_value = {"Body": MagicMock(read=lambda: fake_body)}
        with patch("boto3.client", return_value=mock_s3):
            result = self._resolve("s3://my-bucket/path/to/seqs.fasta")
        assert result == fake_body.decode("utf-8")
        mock_s3.get_object.assert_called_once_with(
            Bucket="my-bucket", Key="path/to/seqs.fasta"
        )

    def test_s3_fetch_failure_returns_original(self):
        """If S3 fetch fails, fall back gracefully (don't crash)."""
        with patch("boto3.client", side_effect=Exception("no credentials")):
            result = self._resolve("s3://bucket/key.fasta")
        assert result == "s3://bucket/key.fasta"

    def test_local_file_path_read(self, tmp_path):
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">s1\nACGT\n>s2\nTGCA")
        result = self._resolve(str(fasta_file))
        assert ">s1" in result

    def test_nonexistent_file_returns_original(self):
        result = self._resolve("/nonexistent/path.fasta")
        assert result == "/nonexistent/path.fasta"

    def test_run_phylogenetic_tree_raw_falls_back_to_builtin_on_empty(self):
        """Empty input → built-in SARS-CoV-2 dataset used, result has status success."""
        from backend.phylogenetic_tree import run_phylogenetic_tree_raw
        result = run_phylogenetic_tree_raw("")
        assert result.get("status") != "error", (
            f"Empty input should fall back to built-in dataset, got: {result}"
        )
        assert "newick" in result or "tree_newick" in result or "text" in result
