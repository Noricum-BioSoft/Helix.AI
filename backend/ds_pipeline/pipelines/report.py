"""
Report generation pipeline step.

Produces:
  - artifacts/{run_id}/report.md       — narrative markdown report
  - artifacts/{run_id}/figures/        — matplotlib PNG plots:
      target_distribution.png
      missingness_bar.png
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional


def generate_report(
    run_data: Dict[str, Any],
    run_dir: Path,
    figures_dir: Path,
) -> str:
    """
    Write report.md and generate plots to figures_dir.

    Returns the path to report.md as a string.
    """
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_paths: List[str] = []

    eda = run_data.get("eda_summary") or {}
    metrics = run_data.get("metrics") or {}
    validation = run_data.get("validation_results") or {}
    comparison = run_data.get("comparison") or {}
    next_steps = run_data.get("next_steps") or []
    ingestion = run_data.get("ingestion_metadata") or {}

    # Generate figures
    figure_paths.extend(_plot_target_distribution(eda, figures_dir))
    figure_paths.extend(_plot_missingness(eda, figures_dir))

    run_data["figure_paths"] = figure_paths

    report_md = _render_markdown(run_data, ingestion, eda, metrics, validation, comparison, next_steps)
    report_path = run_dir / "report.md"
    report_path.write_text(report_md)

    return str(report_path)


def _render_markdown(
    run_data: Dict[str, Any],
    ingestion: Dict[str, Any],
    eda: Dict[str, Any],
    metrics: Dict[str, Any],
    validation: Dict[str, Any],
    comparison: Dict[str, Any],
    next_steps: List[Any],
) -> str:
    run_id = run_data.get("run_id", "unknown")
    ts = run_data.get("timestamp", "")
    objective = run_data.get("objective", "")
    hypothesis = run_data.get("hypothesis", "")
    changes = run_data.get("changes", "")
    decision = run_data.get("decision", "")
    target_col = run_data.get("config", {}).get("target_col") or "N/A"
    task_type = run_data.get("train_info", {}).get("task_type") or "N/A"

    lines: List[str] = [
        f"# Analysis Report — {run_id}",
        f"\n**Generated:** {ts}  \n**Objective:** {objective}  \n**Hypothesis:** {hypothesis}  \n**Changes:** {changes}\n",
    ]

    # Dataset summary
    lines.append("## Dataset Summary\n")
    if ingestion:
        lines.append(f"| Property | Value |")
        lines.append(f"|---|---|")
        lines.append(f"| Rows | {ingestion.get('n_rows', 'N/A')} |")
        lines.append(f"| Columns | {ingestion.get('n_cols', 'N/A')} |")
        lines.append(f"| File size | {_fmt_bytes(ingestion.get('size_bytes', 0))} |")
        lines.append(f"| Target column | {target_col} |")
        lines.append(f"| Task type | {task_type} |\n")

    # Data quality
    lines.append("## Data Quality\n")
    val_summary = validation.get("summary", "")
    if val_summary:
        lines.append(val_summary + "\n")
    if validation.get("results"):
        for r in validation["results"]:
            icon = "🚫" if r.get("level") == "error" else "⚠️"
            lines.append(f"- {icon} **[{r.get('check')}]** {r.get('message')}")
        lines.append("")

    # EDA highlights
    lines.append("## EDA Highlights\n")
    shape = eda.get("shape", [])
    if shape:
        lines.append(f"- **Shape after cleaning:** {shape[0]} rows × {shape[1]} columns")
    n_num = eda.get("n_numeric", 0)
    n_cat = eda.get("n_categorical", 0)
    lines.append(f"- **Numeric columns:** {n_num}  **Categorical columns:** {n_cat}")

    miss = eda.get("missingness", {})
    high_miss = {k: v for k, v in miss.items() if v > 0}
    if high_miss:
        lines.append(f"\n### Remaining Missingness After Cleaning\n")
        for col, pct in sorted(high_miss.items(), key=lambda x: -x[1])[:10]:
            lines.append(f"- `{col}`: {pct:.1%}")
    lines.append("")

    top_corr = eda.get("top_correlations", [])
    if top_corr:
        lines.append(f"### Top Correlations with Target (`{target_col}`)\n")
        lines.append("| Feature | Correlation |")
        lines.append("|---|---|")
        for item in top_corr[:5]:
            lines.append(f"| `{item['feature']}` | {item['correlation']:.4f} |")
        lines.append("")

    # Baseline metrics
    if metrics:
        lines.append("## Baseline Model Metrics\n")
        lines.append("| Metric | Value |")
        lines.append("|---|---|")
        for k, v in metrics.items():
            if isinstance(v, (int, float)):
                lines.append(f"| {k} | {v:.4f} |")
        lines.append("")

    # Comparison vs prior best
    if comparison and comparison.get("prior_best"):
        lines.append("## Comparison vs Prior Best\n")
        deltas = comparison.get("deltas", {})
        pm = comparison.get("primary_metric", "")
        is_new = comparison.get("is_new_best", False)
        lines.append(f"**New best:** {'✅ Yes' if is_new else '❌ No'}  \n**Primary metric:** {pm}\n")
        if deltas:
            lines.append("| Metric | Δ vs Prior Best |")
            lines.append("|---|---|")
            for k, v in deltas.items():
                arrow = "↑" if v > 0 else ("↓" if v < 0 else "→")
                lines.append(f"| {k} | {arrow} {v:+.4f} |")
        lines.append("")

    # Next steps
    if next_steps:
        lines.append("## Next Steps\n")
        for i, step in enumerate(next_steps, 1):
            if isinstance(step, dict):
                name = step.get("name", f"Option {i}")
                desc = step.get("description", "")
                score = step.get("score")
                score_str = f" *(score: {score:.2f})*" if score is not None else ""
                lines.append(f"{i}. **{name}**{score_str} — {desc}")
            else:
                lines.append(f"{i}. {step}")
        lines.append("")

    # Decision
    if decision:
        lines.append(f"## Decision\n\n**{decision}**\n")

    # Figures
    figure_paths = run_data.get("figure_paths", [])
    if figure_paths:
        lines.append("## Figures\n")
        for fp in figure_paths:
            name = Path(fp).stem
            lines.append(f"- `{Path(fp).name}`")
        lines.append("")

    return "\n".join(lines)


def _plot_target_distribution(eda: Dict[str, Any], figures_dir: Path) -> List[str]:
    target_dist = eda.get("target_distribution")
    if not target_dist:
        return []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7, 4))
        if target_dist.get("type") == "categorical":
            counts = target_dist.get("value_counts", {})
            labels = list(counts.keys())[:20]
            values = [counts[k] for k in labels]
            ax.bar(range(len(labels)), values)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels([str(l) for l in labels], rotation=45, ha="right")
            ax.set_title("Target Distribution")
            ax.set_ylabel("Count")
        else:
            hist = target_dist.get("histogram", {})
            counts = hist.get("counts", [])
            edges = hist.get("edges", [])
            if counts and edges:
                widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                ax.bar(edges[:-1], counts, width=widths, align="edge")
                ax.set_title("Target Distribution")
                ax.set_xlabel("Value")
                ax.set_ylabel("Count")

        fig.tight_layout()
        out = figures_dir / "target_distribution.png"
        fig.savefig(out, dpi=120)
        plt.close(fig)
        return [str(out)]
    except Exception:
        return []


def _plot_missingness(eda: Dict[str, Any], figures_dir: Path) -> List[str]:
    missingness = eda.get("missingness", {})
    missing_cols = {k: v for k, v in missingness.items() if v > 0}
    if not missing_cols:
        return []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        cols = list(missing_cols.keys())[:20]
        vals = [missing_cols[c] * 100 for c in cols]
        fig, ax = plt.subplots(figsize=(max(6, len(cols) * 0.5 + 2), 4))
        ax.barh(range(len(cols)), vals)
        ax.set_yticks(range(len(cols)))
        ax.set_yticklabels(cols)
        ax.set_xlabel("Missing (%)")
        ax.set_title("Missingness by Column")
        fig.tight_layout()
        out = figures_dir / "missingness_bar.png"
        fig.savefig(out, dpi=120)
        plt.close(fig)
        return [str(out)]
    except Exception:
        return []


def _fmt_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"
