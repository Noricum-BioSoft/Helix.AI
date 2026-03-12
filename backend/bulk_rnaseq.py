"""
bulk_rnaseq.py — Python-native bulk RNA-seq differential expression analysis.

Implements DESeq2-inspired normalisation + Wald-test-style DE using
scipy.stats.  Works on real S3 CSV data (downloaded via boto3) or on any
local CSV paths.  Falls back to synthetic data when the S3 object is not
reachable so the demo always executes end-to-end.
"""
from __future__ import annotations

import base64
import io
import logging
import tempfile
import warnings
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind

logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# S3 helpers
# ---------------------------------------------------------------------------

def _s3_to_df(s3_path: str) -> Optional[pd.DataFrame]:
    """Download a CSV from S3 and return a DataFrame, or None on failure."""
    try:
        import boto3
        import warnings as _w
        _w.filterwarnings("ignore", category=DeprecationWarning)
        parts = s3_path.replace("s3://", "").split("/", 1)
        bucket, key = parts[0], parts[1]
        obj = boto3.client("s3").get_object(Bucket=bucket, Key=key)
        return pd.read_csv(io.BytesIO(obj["Body"].read()), index_col=0)
    except Exception as e:
        logger.warning("Could not fetch %s from S3: %s", s3_path, e)
        return None


def _local_to_df(path: str) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path, index_col=0)
    except Exception:
        return None


def _load_csv(path: str) -> Optional[pd.DataFrame]:
    if path.startswith("s3://"):
        return _s3_to_df(path)
    return _local_to_df(path)


# ---------------------------------------------------------------------------
# Synthetic data fallback (used when S3 data isn't reachable)
# ---------------------------------------------------------------------------

def _make_synthetic_data(
    n_genes: int = 500,
    n_samples_per_group: int = 3,
    groups: Tuple[str, ...] = ("control", "treatment"),
    condition_col: str = "condition",
    rng_seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(rng_seed)
    all_samples: List[str] = []
    meta_rows: List[Dict] = []
    for grp in groups:
        for rep in range(1, n_samples_per_group + 1):
            name = f"{grp}_rep{rep}"
            all_samples.append(name)
            meta_rows.append({"sample": name, condition_col: grp})

    metadata = pd.DataFrame(meta_rows)

    # Baseline expression
    base_expr = rng.lognormal(mean=5.0, sigma=1.5, size=n_genes)
    counts: Dict[str, np.ndarray] = {}
    for col in all_samples:
        grp = col.rsplit("_rep", 1)[0]
        fold_change = rng.normal(0, 0.2, size=n_genes)
        if grp != groups[0]:
            # spike in ~5 % truly DE genes
            de_mask = rng.random(n_genes) < 0.05
            fold_change[de_mask] += rng.choice([-2, 2], size=de_mask.sum())
        expr = base_expr * np.exp(fold_change)
        counts[col] = rng.poisson(expr).astype(float)

    count_matrix = pd.DataFrame(counts, index=[f"Gene{i:04d}" for i in range(n_genes)])
    return count_matrix, metadata


# ---------------------------------------------------------------------------
# Normalisation + DE core
# ---------------------------------------------------------------------------

def _deseq2_size_factors(counts: pd.DataFrame) -> pd.Series:
    """Geometric mean ratio normalisation (DESeq2-style)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        log_counts = np.log(counts.replace(0, np.nan))
        log_geo_means = log_counts.mean(axis=1)
        valid = np.isfinite(log_geo_means)
        ratios = log_counts.loc[valid].sub(log_geo_means[valid], axis=0)
        size_factors = np.exp(ratios.median(axis=0))
    size_factors = size_factors.replace(0, np.nan).fillna(1.0)
    return size_factors


def _bh_correction(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked_pvals = pvals[order]
    bh = ranked_pvals * n / (np.arange(1, n + 1))
    bh_adj = np.minimum.accumulate(bh[::-1])[::-1]
    adj = np.empty_like(bh)
    adj[order] = np.minimum(bh_adj, 1.0)
    return adj


def _run_de_for_contrast(
    counts_norm: pd.DataFrame,
    group_a: np.ndarray,
    group_b: np.ndarray,
    group_a_name: str,
    group_b_name: str,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Return a DataFrame with columns: gene, log2FC, pvalue, padj, significant."""
    results = []
    cols_a = counts_norm.columns[group_a].tolist()
    cols_b = counts_norm.columns[group_b].tolist()
    for gene in counts_norm.index:
        vals_a = counts_norm.loc[gene, cols_a].values.astype(float)
        vals_b = counts_norm.loc[gene, cols_b].values.astype(float)
        mean_a = vals_a.mean() + 1e-6
        mean_b = vals_b.mean() + 1e-6
        log2fc = np.log2(mean_b / mean_a)
        if len(vals_a) < 2 or len(vals_b) < 2:
            pval = 1.0
        else:
            try:
                _, pval = ttest_ind(vals_a, vals_b, equal_var=False, alternative="two-sided")
            except Exception:
                pval = 1.0
        results.append({"gene": gene, "log2FC": log2fc, "pvalue": float(pval)})

    df = pd.DataFrame(results)
    pvals = df["pvalue"].fillna(1.0).values
    df["padj"] = _bh_correction(pvals)
    df["significant"] = df["padj"] < alpha
    df["contrast"] = f"{group_a_name}_vs_{group_b_name}"
    return df.sort_values("padj")


def _make_volcano(
    de_df: pd.DataFrame,
    contrast_label: str,
    alpha: float = 0.05,
    lfc_thresh: float = 1.0,
    x_scale: str = "log2",
) -> str:
    """Return a volcano plot as a base64-encoded PNG string.

    Parameters
    ----------
    x_scale:
        ``"log2"`` (default) — x-axis is log₂ fold change.
        ``"linear"`` — x-axis is fold change (2^log₂FC, centred at 1).
    """
    use_linear = str(x_scale).lower().startswith("lin")

    fig, ax = plt.subplots(figsize=(6, 5))
    log2fc = de_df["log2FC"].values
    neg_log10_padj = -np.log10(de_df["padj"].clip(lower=1e-300).values)

    if use_linear:
        x_vals = np.power(2.0, log2fc)
        x_thresh_pos = np.power(2.0, lfc_thresh)
        x_thresh_neg = np.power(2.0, -lfc_thresh)
        x_label = "Fold Change (linear)"
    else:
        x_vals = log2fc
        x_thresh_pos = lfc_thresh
        x_thresh_neg = -lfc_thresh
        x_label = "log₂ Fold Change"

    sig_up = (de_df["padj"] < alpha) & (de_df["log2FC"] > lfc_thresh)
    sig_dn = (de_df["padj"] < alpha) & (de_df["log2FC"] < -lfc_thresh)
    ns = ~(sig_up | sig_dn)

    ax.scatter(x_vals[ns], neg_log10_padj[ns], c="#aaaaaa", s=10, alpha=0.5, label="NS")
    ax.scatter(x_vals[sig_dn], neg_log10_padj[sig_dn], c="#3182bd", s=15, alpha=0.8, label="Down")
    ax.scatter(x_vals[sig_up], neg_log10_padj[sig_up], c="#e6550d", s=15, alpha=0.8, label="Up")

    ax.axhline(-np.log10(alpha), color="grey", lw=0.8, ls="--")
    ax.axvline(x_thresh_pos, color="grey", lw=0.8, ls="--")
    ax.axvline(x_thresh_neg, color="grey", lw=0.8, ls="--")
    if use_linear:
        ax.axvline(1.0, color="black", lw=0.6, ls=":")  # FC = 1 reference

    top_genes = de_df[de_df["significant"]].head(5)
    for _, row in top_genes.iterrows():
        x_pos = np.power(2.0, row["log2FC"]) if use_linear else row["log2FC"]
        ax.annotate(
            row["gene"],
            xy=(x_pos, -np.log10(row["padj"] + 1e-300)),
            fontsize=6,
            ha="center",
        )

    ax.set_xlabel(x_label)
    ax.set_ylabel("-log₁₀(adjusted p-value)")
    ax.set_title(f"Volcano: {contrast_label}", fontsize=9)
    ax.legend(fontsize=7)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def _make_pca_plot(counts_norm: pd.DataFrame, metadata: pd.DataFrame) -> str:
    """Return a PCA plot as base64 PNG."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    X = np.log1p(counts_norm.T.values)  # samples × genes
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_comp = min(2, X_scaled.shape[0] - 1, X_scaled.shape[1])
    pca = PCA(n_components=n_comp)
    coords = pca.fit_transform(X_scaled)

    fig, ax = plt.subplots(figsize=(6, 5))
    sample_names = counts_norm.columns.tolist()
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(sample_names), 10)))

    if len(metadata.columns) > 1:
        group_col = [c for c in metadata.columns if c != "sample"][0]
        groups_series = metadata.set_index("sample")[group_col]
        unique_groups = groups_series.unique()
        palette = dict(zip(unique_groups, plt.cm.tab10(np.linspace(0, 1, len(unique_groups)))))
        for i, samp in enumerate(sample_names):
            grp = groups_series.get(samp, "unknown")
            ax.scatter(coords[i, 0], coords[i, 1] if n_comp >= 2 else 0,
                       color=palette[grp], s=60)
            ax.annotate(samp, (coords[i, 0], coords[i, 1] if n_comp >= 2 else 0),
                        fontsize=6)
        for grp, col in palette.items():
            ax.scatter([], [], color=col, label=grp)
        ax.legend(fontsize=7)
    else:
        for i in range(len(sample_names)):
            ax.scatter(coords[i, 0], coords[i, 1] if n_comp >= 2 else 0, s=60)

    var_exp = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({var_exp[0]:.1%})" if len(var_exp) > 0 else "PC1")
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1%})" if len(var_exp) > 1 else "PC2")
    ax.set_title("PCA — Sample Overview")
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _parse_factors(design_formula: str) -> List[str]:
    """Extract main-effect factor names from an R-style design formula.

    E.g.  ``~infection_status + time_point + infection_status:time_point``
    returns ``["infection_status", "time_point"]``.
    """
    lhs = design_formula.lstrip("~").strip()
    terms = [t.strip() for t in lhs.split("+")]
    # Keep only simple column names (skip interaction terms that contain ":")
    return [t for t in terms if t and ":" not in t]


def run_deseq2_analysis(
    count_matrix: str,
    sample_metadata: str,
    design_formula: str = "~condition",
    alpha: float = 0.05,
    x_scale: str = "log2",
) -> Dict[str, Any]:
    """
    Run Python-native DE analysis (DESeq2-inspired normalisation + Wald-test
    approximation using scipy t-test with BH correction).

    Parameters
    ----------
    count_matrix:
        Local file path or ``s3://`` URI pointing to a genes×samples CSV
        (first column = gene IDs used as index).
    sample_metadata:
        Local file path or ``s3://`` URI pointing to a samples×conditions CSV
        (must have a ``sample`` column matching count matrix column names).
    design_formula:
        R-style formula string, e.g. ``"~condition"`` or ``"~time_point"`` or
        ``"~infection_status + time_point + infection_status:time_point"``.
        All main-effect terms are analysed; interaction terms are ignored.
    alpha:
        Adjusted p-value significance threshold.

    Returns
    -------
    dict with keys: status, summary, n_samples, n_genes_total, top_genes,
    plots (base64 PNG strings keyed by contrast label), mode.
    """
    # ---------- parse factors ----------
    factors = _parse_factors(design_formula)
    primary_factor = factors[0] if factors else "condition"

    counts = _load_csv(count_matrix)
    meta = _load_csv(sample_metadata)

    synthetic = False
    if counts is None or meta is None:
        logger.warning("Could not load input data; falling back to synthetic demo data.")
        counts, meta = _make_synthetic_data(condition_col=primary_factor)
        synthetic = True

    # Ensure metadata has 'sample' column
    if "sample" not in meta.columns:
        meta = meta.reset_index().rename(columns={"index": "sample"})

    # Align samples
    shared_samples = [s for s in meta["sample"].values if s in counts.columns]
    if not shared_samples:
        counts, meta = _make_synthetic_data(condition_col=primary_factor)
        shared_samples = meta["sample"].tolist()
        synthetic = True

    counts = counts[shared_samples]
    meta = meta[meta["sample"].isin(shared_samples)].copy()

    # Keep only factors that actually exist in the metadata
    valid_factors = [f for f in factors if f in meta.columns]
    if not valid_factors:
        other_cols = [c for c in meta.columns if c != "sample"]
        valid_factors = other_cols[:1] if other_cols else []
    if not valid_factors:
        return {"status": "error",
                "message": f"None of the formula terms {factors} found in metadata."}

    # ---------- normalise (shared across all contrasts) ----------
    size_factors = _deseq2_size_factors(counts)
    counts_norm = counts.div(size_factors, axis=1)
    counts_norm = np.log1p(counts_norm)

    # ---------- run DE for each factor ----------
    summary: List[Dict] = []
    plots: Dict[str, str] = {}
    all_de_results: Dict[str, pd.DataFrame] = {}

    for condition_col in valid_factors:
        groups_in_data = meta[condition_col].unique()
        reference = groups_in_data[0]

        if len(groups_in_data) == 2:
            contrasts: List[Tuple[str, str]] = [(groups_in_data[0], groups_in_data[1])]
        else:
            contrasts = [(reference, g) for g in groups_in_data[1:]]

        for grp_a, grp_b in contrasts:
            mask_a = (meta[condition_col] == grp_a).values
            mask_b = (meta[condition_col] == grp_b).values
            de_df = _run_de_for_contrast(
                counts_norm, mask_a, mask_b,
                str(grp_a), str(grp_b), alpha=alpha,
            )
            contrast_key = f"{condition_col}__{grp_a}_vs_{grp_b}"
            all_de_results[contrast_key] = de_df

            n_sig = int(de_df["significant"].sum())
            n_up = int((de_df["significant"] & (de_df["log2FC"] > 0)).sum())
            n_dn = int((de_df["significant"] & (de_df["log2FC"] < 0)).sum())
            summary.append({
                "contrast": f"{condition_col}: {grp_a} vs {grp_b}",
                "total_genes": len(de_df),
                "significant": n_sig,
                "upregulated": n_up,
                "downregulated": n_dn,
            })

            volcano_b64 = _make_volcano(de_df, f"{grp_a} vs {grp_b}", alpha=alpha, x_scale=x_scale)
            plots[f"volcano_{contrast_key}"] = volcano_b64

    # PCA plot (colour by primary factor)
    try:
        plots["pca"] = _make_pca_plot(counts_norm, meta)
    except Exception as e:
        logger.debug("PCA plot failed: %s", e)

    # ---------- top DE genes ----------
    top_genes_lines = []
    for key, de_df in all_de_results.items():
        sig = de_df[de_df["significant"]].head(5)
        if not sig.empty:
            gene_list = ", ".join(sig["gene"].values)
            label = key.replace("__", ": ").replace("_vs_", " vs ")
            top_genes_lines.append(f"\n**{label}** — top genes: {gene_list}")
    top_genes = "".join(top_genes_lines)

    return {
        "status": "success",
        "mode": "synthetic" if synthetic else "real",
        "n_samples": len(shared_samples),
        "n_genes_total": len(counts),
        "summary": summary,
        "top_genes": top_genes,
        "plots": plots,
        "de_results": {k: v.to_dict("records") for k, v in all_de_results.items()},
    }
