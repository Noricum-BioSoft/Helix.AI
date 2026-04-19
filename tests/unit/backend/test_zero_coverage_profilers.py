"""
Targeted coverage tests for the two 0%-coverage profilers:
  - backend/file_intelligence/singlecell.py
  - backend/file_intelligence/alignment.py

Strategy
--------
* H5AD / HDF5   — anndata + h5py are available; use real in-memory files.
* Loom          — loompy NOT installed; mock the module so all code paths run.
* SAM/BAM       — pysam NOT installed; test the graceful ImportError fallback
                  AND the full read-path via a mock.
"""
from __future__ import annotations

import sys
import types
from pathlib import Path
from unittest.mock import MagicMock, patch, call

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helpers to build synthetic fixture files
# ---------------------------------------------------------------------------

def _make_h5ad(tmp_path: Path) -> Path:
    """Create a minimal AnnData H5AD file."""
    import anndata as ad
    import pandas as pd

    n_cells, n_genes = 20, 50
    X = np.random.default_rng(0).random((n_cells, n_genes)).astype("float32")
    obs = pd.DataFrame(
        {"cell_type": np.random.choice(["T", "B", "NK"], n_cells)},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(
        {"gene_name": [f"GENE_{j}" for j in range(n_genes)]},
        index=[f"gene_{j}" for j in range(n_genes)],
    )
    adata = ad.AnnData(X=X, obs=obs, var=var)

    p = tmp_path / "test.h5ad"
    adata.write_h5ad(p)
    return p


def _make_h5(tmp_path: Path) -> Path:
    """Create a generic HDF5 file (NOT valid AnnData)."""
    import h5py

    p = tmp_path / "test.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("expression", data=np.zeros((10, 5)))
        grp = f.create_group("metadata")
        grp.create_dataset("sample_ids", data=["s1", "s2"])
    return p


# ---------------------------------------------------------------------------
# Tests: singlecell.py — H5AD path (real anndata)
# ---------------------------------------------------------------------------

class TestSingleCellH5AD:
    def test_profile_h5ad_basic_shape(self, tmp_path: Path):
        from backend.file_intelligence.singlecell import profile_singlecell

        h5ad = _make_h5ad(tmp_path)
        result = profile_singlecell(h5ad)

        assert result["format"] == "h5ad"
        assert result["family"] == "single_cell"
        assert result["n_records"] == 20  # n_cells
        assert result["profiler_error"] is None

    def test_profile_h5ad_summary_fields(self, tmp_path: Path):
        from backend.file_intelligence.singlecell import profile_singlecell

        result = profile_singlecell(_make_h5ad(tmp_path))
        summary = result["summary"]

        assert summary["n_cells"] == 20
        assert summary["n_genes"] == 50
        assert isinstance(summary["obs_metadata_keys"], list)
        assert "cell_type" in summary["obs_metadata_keys"]

    def test_profile_h5ad_schema_keys(self, tmp_path: Path):
        from backend.file_intelligence.singlecell import profile_singlecell

        result = profile_singlecell(_make_h5ad(tmp_path))
        assert "obs_columns" in result["schema"]
        assert "var_columns" in result["schema"]

    def test_profile_h5ad_sample_cells(self, tmp_path: Path):
        from backend.file_intelligence.singlecell import profile_singlecell

        result = profile_singlecell(_make_h5ad(tmp_path))
        assert len(result["sample"]) > 0

    def test_profile_h5ad_embeddings_empty_when_none(self, tmp_path: Path):
        from backend.file_intelligence.singlecell import profile_singlecell

        result = profile_singlecell(_make_h5ad(tmp_path))
        assert result["summary"]["has_umap"] is False
        assert result["summary"]["has_pca"] is False

    def test_profile_h5ad_dispatched_via_profiler(self, tmp_path: Path):
        """Verify profile_file() dispatches to singlecell profiler for .h5ad."""
        from backend.file_intelligence.profiler import profile_file

        result = profile_file(str(_make_h5ad(tmp_path)))
        assert result["family"] == "single_cell"
        assert result["n_records"] == 20


# ---------------------------------------------------------------------------
# Tests: singlecell.py — HDF5 generic fallback (real h5py)
# ---------------------------------------------------------------------------

class TestSingleCellHDF5Generic:
    def test_profile_generic_hdf5(self, tmp_path: Path):
        from backend.file_intelligence.singlecell import profile_singlecell

        h5 = _make_h5(tmp_path)
        result = profile_singlecell(h5)

        # Falls through h5ad attempt → generic HDF5 path
        assert result["family"] == "single_cell"
        assert result["format"] in ("hdf5", "h5")
        # HDF5 keys should be listed
        assert "hdf5_keys" in result["summary"]
        assert any("expression" in k for k in result["summary"]["hdf5_keys"])

    def test_profile_generic_hdf5_dispatched_via_profiler(self, tmp_path: Path):
        from backend.file_intelligence.profiler import profile_file

        result = profile_file(str(_make_h5(tmp_path)))
        assert result["family"] == "single_cell"


# ---------------------------------------------------------------------------
# Tests: singlecell.py — Loom path (loompy mocked)
# ---------------------------------------------------------------------------

class TestSingleCellLoom:
    def _make_fake_loompy(self):
        """Build a minimal fake loompy module."""
        loompy = types.ModuleType("loompy")

        class FakeDS:
            shape = (100, 30)  # (n_genes, n_cells) in loom convention
            col_attrs = {"CellID": None, "cell_type": None}
            row_attrs = {"Gene": None, "Accession": None}
            layers = {"": None, "spliced": None}

            def __enter__(self):
                return self

            def __exit__(self, *_):
                pass

        def connect(path, mode="r"):
            return FakeDS()

        loompy.connect = connect
        return loompy

    def test_profile_loom_basic(self, tmp_path: Path):
        from backend.file_intelligence import singlecell

        fake_loom_path = tmp_path / "test.loom"
        fake_loom_path.touch()

        fake_loompy = self._make_fake_loompy()
        with patch.dict(sys.modules, {"loompy": fake_loompy}):
            result = singlecell._profile_loom(fake_loom_path)

        assert result["format"] == "loom"
        assert result["family"] == "single_cell"
        assert result["n_records"] == 30   # n_cells
        assert result["summary"]["n_genes"] == 100
        assert "CellID" in result["summary"]["cell_attributes"]
        assert "Gene" in result["summary"]["gene_attributes"]
        assert "spliced" in result["summary"]["layers"]
        assert result["profiler_error"] is None

    def test_profile_loom_dispatched_via_profile_singlecell(self, tmp_path: Path):
        from backend.file_intelligence import singlecell

        fake_loom_path = tmp_path / "test.loom"
        fake_loom_path.touch()

        fake_loompy = self._make_fake_loompy()
        with patch.dict(sys.modules, {"loompy": fake_loompy}):
            result = singlecell.profile_singlecell(fake_loom_path)

        assert result["family"] == "single_cell"
        assert result["n_records"] == 30


# ---------------------------------------------------------------------------
# Tests: alignment.py — ImportError fallback (no pysam)
# ---------------------------------------------------------------------------

class TestAlignmentNoPysam:
    def test_bam_no_pysam_returns_graceful_fallback(self, tmp_path: Path):
        from backend.file_intelligence.alignment import profile_alignment

        fake_bam = tmp_path / "test.bam"
        fake_bam.write_bytes(b"\x00" * 16)  # dummy content

        # pysam is not installed in this environment
        result = profile_alignment(fake_bam)

        assert result["family"] == "alignment"
        assert result["format"] == "bam"
        assert result["profiler_error"] == "pysam not available"
        assert "pysam" in result["summary"]["note"].lower()

    def test_sam_no_pysam_returns_graceful_fallback(self, tmp_path: Path):
        from backend.file_intelligence.alignment import profile_alignment

        fake_sam = tmp_path / "test.sam"
        fake_sam.write_text("@HD\tVN:1.6\n")

        result = profile_alignment(fake_sam)

        assert result["family"] == "alignment"
        assert result["format"] == "sam"
        assert result["profiler_error"] == "pysam not available"


# ---------------------------------------------------------------------------
# Tests: alignment.py — Full read path (pysam mocked)
# ---------------------------------------------------------------------------

class TestAlignmentWithPysam:
    def _make_fake_pysam(self, n_reads: int = 150, n_mapped: int = 130):
        """Build a minimal fake pysam module."""
        pysam = types.ModuleType("pysam")

        class FakeRead:
            def __init__(self, i: int):
                self.query_name = f"read_{i}"
                self.reference_name = "chr1"
                self.reference_start = i * 100
                self.mapping_quality = 60
                self.query_length = 150
                self.is_paired = True
                self.is_unmapped = i >= n_mapped

        class FakeHeader:
            def to_dict(self):
                return {"SQ": [{"SN": "chr1", "LN": 248956422}], "PG": [{"ID": "bwa"}]}

        class FakeAlignmentFile:
            references = ("chr1", "chr2")
            lengths = (248956422, 242193529)
            header = FakeHeader()

            def __init__(self, path, mode, check_sq=True):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *_):
                pass

            def fetch(self, until_eof=False):
                for i in range(n_reads):
                    yield FakeRead(i)

        pysam.AlignmentFile = FakeAlignmentFile
        return pysam

    def test_bam_full_profile(self, tmp_path: Path):
        from backend.file_intelligence import alignment

        fake_bam = tmp_path / "test.bam"
        fake_bam.write_bytes(b"\x00" * 16)

        fake_pysam = self._make_fake_pysam(n_reads=150, n_mapped=130)
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = alignment.profile_alignment(fake_bam)

        assert result["format"] == "bam"
        assert result["family"] == "alignment"
        assert result["n_records"] == 150
        assert result["summary"]["n_mapped"] == 130
        assert result["summary"]["n_unmapped"] == 20
        assert result["summary"]["pct_mapped"] == pytest.approx(86.67, abs=0.1)
        assert "chr1" in result["summary"]["references"]
        assert result["profiler_error"] is None

    def test_sam_full_profile(self, tmp_path: Path):
        from backend.file_intelligence import alignment

        fake_sam = tmp_path / "test.sam"
        fake_sam.write_text("@HD\tVN:1.6\n")

        fake_pysam = self._make_fake_pysam(n_reads=50, n_mapped=48)
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = alignment.profile_alignment(fake_sam)

        assert result["format"] == "sam"
        assert result["n_records"] == 50
        assert result["summary"]["n_mapped"] == 48

    def test_sample_reads_capped_at_5(self, tmp_path: Path):
        from backend.file_intelligence import alignment

        fake_bam = tmp_path / "test.bam"
        fake_bam.write_bytes(b"\x00" * 16)

        fake_pysam = self._make_fake_pysam(n_reads=100)
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = alignment.profile_alignment(fake_bam)

        assert len(result["sample"]) == 5  # always capped

    def test_read_count_capped_at_10000(self, tmp_path: Path):
        from backend.file_intelligence import alignment

        fake_bam = tmp_path / "test.bam"
        fake_bam.write_bytes(b"\x00" * 16)

        fake_pysam = self._make_fake_pysam(n_reads=15000, n_mapped=14000)
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = alignment.profile_alignment(fake_bam)

        assert result["n_records"] == 10000  # hard cap

    def test_header_metadata_captured(self, tmp_path: Path):
        from backend.file_intelligence import alignment

        fake_bam = tmp_path / "test.bam"
        fake_bam.write_bytes(b"\x00" * 16)

        fake_pysam = self._make_fake_pysam()
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = alignment.profile_alignment(fake_bam)

        assert result["raw_metadata"]["header_SQ_count"] == 1
        assert "bwa" in result["raw_metadata"]["program_records"]

    def test_dispatched_via_profiler(self, tmp_path: Path):
        from backend.file_intelligence import profiler, alignment

        fake_bam = tmp_path / "test.bam"
        fake_bam.write_bytes(b"\x00" * 16)

        fake_pysam = self._make_fake_pysam()
        with patch.dict(sys.modules, {"pysam": fake_pysam}):
            result = profiler.profile_file(str(fake_bam))

        assert result["family"] == "alignment"
