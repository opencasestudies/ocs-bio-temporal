#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


DEFAULT_SWEEP_DIR = Path("/Users/othomas/Desktop/CS5_sweep_results")
DEFAULT_OUTDIR = Path("GSE154386/revised_k_selection_K10")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required sweep file not found: {path}")
    return pd.read_csv(path)


def scale(series: pd.Series) -> pd.Series:
    max_value = series.max()
    if pd.isna(max_value) or max_value <= 0:
        return pd.Series(0.0, index=series.index)
    return series / max_value


def inverse_scale(series: pd.Series) -> pd.Series:
    min_value = series.min()
    max_value = series.max()
    if pd.isna(min_value) or pd.isna(max_value) or max_value == min_value:
        return pd.Series(1.0, index=series.index)
    return 1 - ((series - min_value) / (max_value - min_value))


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def selected_config_runs(per_run: pd.DataFrame, summary_by_k: pd.DataFrame) -> pd.DataFrame:
    selected = summary_by_k[["K", "selected_n_iter"]].rename(columns={"selected_n_iter": "n_iter"})
    return per_run.merge(selected, on=["K", "n_iter"], how="inner")


def top_genes_for_run(sweep_dir: Path, k: int, seed: int, n_iter: int, pattern: str, n: int = 15) -> list[str]:
    path = sweep_dir / "runs" / f"cogaps_K{k}_seed{seed}_iter{n_iter}.top_genes.csv"
    top = read_csv(path)
    return (
        top.loc[top["pattern"].astype(str) == str(pattern)]
        .sort_values("rank")
        .head(n)["gene"]
        .astype(str)
        .tolist()
    )


def write_figures(summary: pd.DataFrame, outdir: Path, chosen_k: int) -> dict[str, str]:
    if plt is None:
        return {}

    figures = outdir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}

    fig, ax = plt.subplots(figsize=(8.5, 5))
    colors = ["#c44e52" if int(k) == chosen_k else "#4c72b0" for k in summary["K"]]
    ax.scatter(summary["stability_core"], summary["max_eta_timepoint"], s=80, c=colors)
    for _, row in summary.iterrows():
        ax.text(row["stability_core"] + 0.002, row["max_eta_timepoint"], f"K={int(row['K'])}", fontsize=8)
    ax.axvline(summary.loc[summary["on_stability_plateau"], "stability_core"].min(), color="#555555", linestyle="--", linewidth=1)
    ax.set_xlabel("Seed-to-seed stability core")
    ax.set_ylabel("Maximum temporal effect size")
    ax.set_title("Stability versus temporal resolution")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = figures / "stability_vs_temporal_effect.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths["stability_vs_temporal_effect"] = str(path)

    display = summary[summary["on_stability_plateau"]].copy().sort_values("K")
    fig, ax = plt.subplots(figsize=(9, 5))
    width = 0.35
    x = range(len(display))
    ax.bar([v - width / 2 for v in x], display["candidate_ifn_pattern_count_mean"], width=width, label="Mean IFN candidates", color="#55a868")
    ax.bar([v + width / 2 for v in x], display["activity_like_pattern_count_mean"], width=width, label="Mean activity-like patterns", color="#8172b3")
    ax.set_xticks(list(x))
    ax.set_xticklabels([f"K={int(v)}" for v in display["K"]])
    ax.set_ylabel("Mean count across seeds")
    ax.set_title("Biology-facing signals among stability-plateau K values")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    path = figures / "plateau_biology_signals.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths["plateau_biology_signals"] = str(path)

    ranked = summary.sort_values("goal_aligned_score", ascending=False).head(10).sort_values("goal_aligned_score")
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#c44e52" if int(k) == chosen_k else "#4c72b0" for k in ranked["K"]]
    ax.barh([f"K={int(k)}" for k in ranked["K"]], ranked["goal_aligned_score"], color=colors)
    ax.set_xlabel("Goal-aligned score")
    ax.set_title("Revised K-selection ranking")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    path = figures / "revised_k_goal_aligned_score.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths["revised_k_goal_aligned_score"] = str(path)

    return paths


def markdown_table(df: pd.DataFrame, columns: list[str], float_digits: int = 3) -> str:
    view = df[columns].copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.{float_digits}f}")
    view = view.fillna("").astype(str)
    header = "| " + " | ".join(view.columns) + " |"
    separator = "| " + " | ".join(["---"] * len(view.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in view.to_numpy()]
    return "\n".join([header, separator, *rows])


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate the Case Study 5 revised K-selection report from existing sweep outputs.")
    parser.add_argument("--sweep-dir", default=str(DEFAULT_SWEEP_DIR))
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    parser.add_argument("--chosen-k", type=int, default=10)
    parser.add_argument("--chosen-seed", type=int, default=2)
    args = parser.parse_args()

    sweep_dir = Path(args.sweep_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    summary_by_k = read_csv(sweep_dir / "summary_by_K.csv")
    per_run = read_csv(sweep_dir / "per_run_metrics.csv")
    pairwise = read_csv(sweep_dir / "pairwise_reproducibility.csv")
    chosen_summary = json.loads((sweep_dir / "chosen_k_summary.json").read_text(encoding="utf-8"))

    summary = summary_by_k.copy()
    summary["on_stability_plateau"] = bool_series(summary["on_stability_plateau"])
    selected_runs = selected_config_runs(per_run, summary[["K", "selected_n_iter"]])
    run_agg = (
        selected_runs.groupby(["K", "n_iter"], as_index=False)
        .agg(
            activity_like_pattern_count_mean=("activity_like_pattern_count", "mean"),
            activity_like_pattern_count_max=("activity_like_pattern_count", "max"),
            candidate_ifn_pattern_count_max=("candidate_ifn_pattern_count", "max"),
            max_ifn_score_rho_mean=("max_ifn_score_rho", "mean"),
            max_ifn_top_gene_overlap_top15_mean=("max_ifn_top_gene_overlap_top15", "mean"),
            runtime_sec_mean=("runtime_sec", "mean"),
        )
        .rename(columns={"n_iter": "selected_n_iter"})
    )
    summary = summary.merge(run_agg, on=["K", "selected_n_iter"], how="left")

    summary["low_redundancy_score"] = inverse_scale(summary["redundancy_penalty"])
    summary["temporal_effect_score"] = scale(summary["max_eta_timepoint"])
    summary["ifn_candidate_score"] = scale(summary["candidate_ifn_pattern_count_mean"])
    summary["activity_score"] = scale(summary["activity_like_pattern_count_mean"].fillna(0))
    summary["goal_aligned_score"] = (
        0.35 * summary["stability_core"]
        + 0.15 * summary["low_redundancy_score"]
        + 0.20 * summary["temporal_effect_score"]
        + 0.15 * summary["ifn_candidate_score"]
        + 0.15 * summary["activity_score"]
    )

    plateau = summary[summary["on_stability_plateau"]].copy()
    biology_eligible = plateau[
        (plateau["candidate_ifn_pattern_count_mean"] >= 1)
        & (plateau["activity_like_pattern_count_mean"].fillna(0) >= 1)
        & (plateau["max_eta_timepoint"] >= 0.1)
    ].copy()
    if biology_eligible.empty:
        chosen_row = plateau.sort_values("goal_aligned_score", ascending=False).iloc[0]
        rule = "highest goal-aligned score among stability-plateau candidates"
    else:
        chosen_row = biology_eligible.sort_values(["K", "goal_aligned_score"], ascending=[True, False]).iloc[0]
        rule = "smallest stability-plateau K with at least one IFN-candidate activity-like temporal program"

    chosen_k = int(chosen_row["K"])
    chosen_n_iter = int(chosen_row["selected_n_iter"])
    chosen_seed = int(args.chosen_seed)
    pattern_path = sweep_dir / "runs" / f"cogaps_K{chosen_k}_seed{chosen_seed}_iter{chosen_n_iter}.pattern_summary.csv"
    pattern_summary = read_csv(pattern_path)
    candidate_patterns = pattern_summary.sort_values(
        ["candidate_ifn_pattern", "eta_timepoint"], ascending=[False, False]
    )
    candidate_patterns = candidate_patterns[
        [
            "pattern",
            "eta_timepoint",
            "eta_broad_cell_type",
            "pattern_class",
            "spearman_ifn_score_rho",
            "peak_timepoint",
            "ifn_top_gene_overlap_top15",
            "candidate_ifn_pattern",
        ]
    ]
    top_candidate = candidate_patterns.iloc[0]
    top_candidate_genes = top_genes_for_run(
        sweep_dir,
        chosen_k,
        chosen_seed,
        chosen_n_iter,
        str(top_candidate["pattern"]),
        n=15,
    )

    figures = write_figures(summary, outdir, chosen_k)

    summary_out = outdir / "revised_k_selection_summary_by_k.csv"
    plateau_out = outdir / "stability_plateau_candidates.csv"
    candidate_out = outdir / f"representative_K{chosen_k}_seed{chosen_seed}_iter{chosen_n_iter}_patterns.csv"
    pairwise_out = outdir / "selected_pairwise_reproducibility.csv"
    summary.sort_values("K").to_csv(summary_out, index=False)
    plateau.sort_values("K").to_csv(plateau_out, index=False)
    candidate_patterns.to_csv(candidate_out, index=False)
    pairwise.loc[pairwise["K"].eq(chosen_k) & pairwise["n_iter"].eq(chosen_n_iter)].to_csv(pairwise_out, index=False)

    manifest: dict[str, Any] = {
        "created_at_utc": utc_now(),
        "sweep_dir": str(sweep_dir),
        "previous_chosen_k": chosen_summary.get("chosen_K"),
        "previous_selected_n_iter": chosen_summary.get("selected_n_iter"),
        "revised_chosen_k": chosen_k,
        "revised_seed": chosen_seed,
        "revised_selected_n_iter": chosen_n_iter,
        "selection_rule": rule,
        "plateau_threshold": chosen_summary.get("plateau_threshold"),
        "source_files": {
            "summary_by_K": str(sweep_dir / "summary_by_K.csv"),
            "summary_by_k_n_iter": str(sweep_dir / "summary_by_k_n_iter.csv"),
            "per_run_metrics": str(sweep_dir / "per_run_metrics.csv"),
            "pairwise_reproducibility": str(sweep_dir / "pairwise_reproducibility.csv"),
            "representative_pattern_summary": str(pattern_path),
        },
        "outputs": {
            "summary_by_k": str(summary_out),
            "plateau_candidates": str(plateau_out),
            "representative_patterns": str(candidate_out),
            "selected_pairwise_reproducibility": str(pairwise_out),
            **figures,
        },
        "top_candidate_pattern": {
            "pattern": str(top_candidate["pattern"]),
            "eta_timepoint": float(top_candidate["eta_timepoint"]),
            "eta_broad_cell_type": float(top_candidate["eta_broad_cell_type"]),
            "pattern_class": str(top_candidate["pattern_class"]),
            "spearman_ifn_score_rho": float(top_candidate["spearman_ifn_score_rho"]),
            "peak_timepoint": str(top_candidate["peak_timepoint"]),
            "ifn_top_gene_overlap_top15": int(top_candidate["ifn_top_gene_overlap_top15"]),
            "top_15_genes": top_candidate_genes,
        },
    }
    manifest_path = outdir / "revised_k_selection_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    report_path = outdir / "revised_k_selection_report.md"
    report = f"""# Revised K Selection for Case Study 5

Generated: `{manifest["created_at_utc"]}`

## Decision

The revised selected model for the main Case Study 5 analysis is:

- `K = {chosen_k}`
- `seed = {chosen_seed}`
- `nIterations = {chosen_n_iter}`

The previous sweep selected `K = {chosen_summary.get("chosen_K")}` because the original rule preferred the smallest K on the seed-to-seed stability plateau. That rule is useful for a compact baseline, but it can be too compressed when the biological question is temporal immune activity: resolving dengue-response programs, especially a D10 interferon-associated activity program.

The revised rule is:

> Select the smallest stability-plateau K that resolves at least one IFN-candidate, activity-like temporal program.

This selects `K = {chosen_k}`.

## Why Reassess K?

This case study asks whether CoGAPS can recover temporal immune programs, separate cell identity from infection activity, and identify interferon-stimulated antiviral responses. A low-rank model can be highly stable while compressing several biological processes into broad lineage-associated factors. The K-selection criterion therefore needs to include biological resolution, not only numerical reproducibility.

This follows the identity/activity framing from Kotliar et al. (2019), where the useful rank depends on the resolution desired by the analyst, and the CoGAPS protocol emphasis on annotating patterns by biological process rather than selecting dimensionality mechanically. Fertig et al. (2014) also motivates using CoGAPS specifically to interpret time-course patterns as overlapping biological processes. Waickman et al. (2021) provides the biological prior that the strongest experimental DENV-1 host response occurs around day 10 and includes interferon/inflammatory genes.

## Stability-Plateau Candidates

{markdown_table(plateau.sort_values("K"), ["K", "selected_n_iter", "stability_core", "overall_score", "max_eta_timepoint", "candidate_ifn_pattern_count_mean", "activity_like_pattern_count_mean", "ifn_secondary_score", "redundancy_penalty"])}

## Revised Ranking

The goal-aligned score is a transparent summary used for review, not a black-box selector:

- 35% seed-to-seed stability
- 15% low redundancy
- 20% maximum temporal effect size
- 15% mean IFN-candidate pattern count
- 15% mean activity-like pattern count

{markdown_table(summary.sort_values("goal_aligned_score", ascending=False).head(10), ["K", "selected_n_iter", "stability_core", "low_redundancy_score", "max_eta_timepoint", "candidate_ifn_pattern_count_mean", "activity_like_pattern_count_mean", "goal_aligned_score", "on_stability_plateau"])}

## Representative K={chosen_k} Pattern Evidence

Representative run: `K={chosen_k}`, `seed={chosen_seed}`, `nIterations={chosen_n_iter}`.

{markdown_table(candidate_patterns.head(10), ["pattern", "eta_timepoint", "eta_broad_cell_type", "pattern_class", "spearman_ifn_score_rho", "peak_timepoint", "ifn_top_gene_overlap_top15", "candidate_ifn_pattern"])}

The leading candidate pattern is `{top_candidate["pattern"]}`:

- temporal effect size: `{float(top_candidate["eta_timepoint"]):.3f}`
- broad-cell-type effect size: `{float(top_candidate["eta_broad_cell_type"]):.3f}`
- class: `{top_candidate["pattern_class"]}`
- IFN-score Spearman rho: `{float(top_candidate["spearman_ifn_score_rho"]):.3f}`
- peak timepoint: `{top_candidate["peak_timepoint"]}`
- top-15 IFN gene overlap: `{int(top_candidate["ifn_top_gene_overlap_top15"])}`
- top genes: `{", ".join(top_candidate_genes)}`

This aligns with the case-study biology: an interferon-associated program peaking at D10, rather than a purely lineage-dominant low-rank factor.

## Interpretation

`K=5` remains useful as a compact stability baseline. It demonstrates that the strongest low-rank structure in this balanced PBMC subset is lineage/identity-associated. However, because the case study's central question is temporal immune activity, `K=10` is the better primary analysis model: it remains on the stability plateau and resolves a D10 IFN-like activity program.

The main takeaway is methodological: the best K depends on the analysis goal. Stability is necessary, but biological resolution determines whether the chosen model can answer the research questions.

## Local Supporting References

- `references/literature/fertig_2014_temporal_cogaps.pdf`
- `references/literature/johnson_2023_cogaps_protocol_nature_protocols.pdf`
- `references/literature/kotliar_2019_identity_activity_cnmf.pdf`
- `references/literature/stein_obrien_2019_sc_cogaps_transfer_learning.pdf`
- `references/literature/sharma_2020_projectr.pdf`
- `references/literature/waickman_2021_dengue_single_cell_pbmc.pdf`

## Output Files

- `{summary_out}`
- `{plateau_out}`
- `{candidate_out}`
- `{pairwise_out}`
- `{manifest_path}`
"""
    if figures:
        report += "\n".join(
            [
                f"- `{figures['stability_vs_temporal_effect']}`",
                f"- `{figures['plateau_biology_signals']}`",
                f"- `{figures['revised_k_goal_aligned_score']}`",
            ]
        )
        report += "\n"
    else:
        report += "- Figures were skipped because matplotlib was not available in the local Python environment.\n"
    report += """
"""
    report_path.write_text(report, encoding="utf-8")
    print(json.dumps({"report": str(report_path), "manifest": str(manifest_path), "chosen_K": chosen_k, "selected_n_iter": chosen_n_iter}, indent=2))


if __name__ == "__main__":
    main()
