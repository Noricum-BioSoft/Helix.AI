"""Upload-time profiler for single-cell formats: H5AD (AnnData), Loom, HDF5."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


def profile_singlecell(path: str | Path) -> Dict[str, Any]:
    path = Path(path)
    suffix = path.suffix.lower()

    try:
        if suffix == ".h5ad":
            return _profile_h5ad(path)
        if suffix == ".loom":
            return _profile_loom(path)
        # Generic HDF5 — attempt h5ad first, fall back to raw key listing
        try:
            return _profile_h5ad(path)
        except Exception:
            return _profile_hdf5_generic(path)
    except Exception as exc:
        return {
            "format": suffix.lstrip("."),
            "family": "single_cell",
            "n_records": None,
            "summary": {},
            "schema": {},
            "sample": [],
            "available_sheets": None,
            "raw_metadata": {},
            "profiler_error": str(exc),
        }


def _profile_h5ad(path: Path) -> Dict[str, Any]:
    import anndata as ad  # optional dependency

    adata = ad.read_h5ad(path, backed="r")  # memory-mapped, don't load full matrix

    n_cells = adata.n_obs
    n_genes = adata.n_vars

    obs_keys = list(adata.obs.columns)
    var_keys = list(adata.var.columns)
    obsm_keys = list(adata.obsm.keys())
    obsp_keys = list(adata.obsp.keys())
    uns_keys = list(adata.uns.keys())

    # Sample: first 5 cells with obs metadata
    sample_cells = adata.obs.head(5).reset_index().to_dict(orient="records")

    # Top expressed genes heuristic (from var if mean available)
    top_genes: list[str] = []
    if "highly_variable" in var_keys:
        hv = adata.var[adata.var["highly_variable"]].index.tolist()
        top_genes = hv[:10]
    elif n_genes <= 50000:
        top_genes = adata.var_names[:10].tolist()

    adata.file.close()

    return {
        "format": "h5ad",
        "family": "single_cell",
        "n_records": n_cells,
        "summary": {
            "n_cells": n_cells,
            "n_genes": n_genes,
            "obs_metadata_keys": obs_keys,
            "var_metadata_keys": var_keys,
            "embeddings": obsm_keys,
            "graphs": obsp_keys,
            "uns_keys": uns_keys,
            "has_umap": "X_umap" in obsm_keys,
            "has_pca": "X_pca" in obsm_keys,
        },
        "schema": {
            "obs_columns": obs_keys,
            "var_columns": var_keys,
        },
        "sample": sample_cells,
        "available_sheets": None,
        "raw_metadata": {"top_genes_preview": top_genes},
        "profiler_error": None,
    }


def _profile_loom(path: Path) -> Dict[str, Any]:
    import loompy  # optional dependency

    with loompy.connect(str(path), mode="r") as ds:
        n_genes, n_cells = ds.shape
        col_attrs = list(ds.col_attrs.keys())
        row_attrs = list(ds.row_attrs.keys())
        layers = list(ds.layers.keys())

    return {
        "format": "loom",
        "family": "single_cell",
        "n_records": n_cells,
        "summary": {
            "n_cells": n_cells,
            "n_genes": n_genes,
            "cell_attributes": col_attrs,
            "gene_attributes": row_attrs,
            "layers": layers,
        },
        "schema": {"cell_attributes": col_attrs, "gene_attributes": row_attrs},
        "sample": [],
        "available_sheets": None,
        "raw_metadata": {},
        "profiler_error": None,
    }


def _profile_hdf5_generic(path: Path) -> Dict[str, Any]:
    import h5py

    keys: list[str] = []
    with h5py.File(path, "r") as f:
        f.visit(keys.append)

    return {
        "format": "hdf5",
        "family": "single_cell",
        "n_records": None,
        "summary": {"hdf5_keys": keys[:30]},
        "schema": {},
        "sample": [],
        "available_sheets": None,
        "raw_metadata": {},
        "profiler_error": None,
    }
