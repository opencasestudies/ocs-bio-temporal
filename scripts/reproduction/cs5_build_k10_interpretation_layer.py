#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


DEFAULT_MODEL_OUTDIR = Path("data/processed/selected_model_k10/r")
DEFAULT_COMPARISON_OUTDIR = Path("data/processed/r_python_comparison")
DEFAULT_K_SELECTION_OUTDIR = Path("data/processed/k_selection")
DEFAULT_DIRECTIONALITY_DIR = Path("data/processed/directionality")
DEFAULT_OUTDIR = Path("data/processed/interpretation_rebuilt")
DEFAULT_RUN_STEM = "cogaps_K10_seed2_iter2000"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def require_columns(table_name: str, columns: Iterable[str], required: Iterable[str]) -> None:
    missing = [column for column in required if column not in columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {', '.join(missing)}")


def pattern_sort_key(pattern: str) -> tuple[int, str]:
    suffix = str(pattern).replace("Pattern", "")
    try:
        return (int(suffix), str(pattern))
    except ValueError:
        return (10_000, str(pattern))


def time_sort_key(row_or_label) -> tuple[float, str]:
    if isinstance(row_or_label, pd.Series):
        label = str(row_or_label["timepoint_merged"])
        day = row_or_label.get("day_merged_numeric", np.nan)
    else:
        label = str(row_or_label)
        day = np.nan
    if pd.isna(day):
        if label.startswith("D"):
            try:
                return (float(label.replace("D", "").replace("/15", ".5")), label)
            except ValueError:
                return (math.inf, label)
        return (math.inf, label)
    return (float(day), label)


def read_csv(path: Path, table_name: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required {table_name} file not found: {path}")
    return pd.read_csv(path)


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def pattern_columns(cells: pd.DataFrame) -> list[str]:
    patterns = [column for column in cells.columns if str(column).startswith("Pattern")]
    if not patterns:
        raise ValueError("No Pattern* columns found in discovery cell table.")
    return sorted(patterns, key=pattern_sort_key)


def group_pattern_means(
    cells: pd.DataFrame,
    group_cols: list[str],
    patterns: list[str],
    extra_mean_cols: list[str] | None = None,
) -> pd.DataFrame:
    extra_mean_cols = [column for column in (extra_mean_cols or []) if column in cells.columns]
    grouped = cells.groupby(group_cols, dropna=False, sort=False)
    counts = grouped.size().rename("n_cells").reset_index()
    means = grouped[patterns + extra_mean_cols].mean(numeric_only=True).reset_index()
    out = counts.merge(means, on=group_cols, how="left")
    if "day_merged_numeric" in out.columns:
        out = out.sort_values(group_cols, key=lambda col: col if col.name != "day_merged_numeric" else col.astype(float))
    return out


def melt_patterns(wide: pd.DataFrame, id_cols: list[str], patterns: list[str]) -> pd.DataFrame:
    long = wide.melt(
        id_vars=id_cols,
        value_vars=patterns,
        var_name="pattern",
        value_name="mean_pattern_score",
    )
    long["pattern"] = pd.Categorical(long["pattern"], categories=patterns, ordered=True)
    return long.sort_values(id_cols + ["pattern"]).reset_index(drop=True)


def subject_summary(
    subject_long: pd.DataFrame,
    cell_long: pd.DataFrame,
    group_cols: list[str],
) -> pd.DataFrame:
    summary = (
        subject_long.groupby(group_cols, dropna=False, observed=True)
        .agg(
            n_subjects=("subject", "nunique"),
            n_subject_groups=("subject", "size"),
            total_cells=("n_cells", "sum"),
            mean_subject_score=("mean_pattern_score", "mean"),
            sd_subject_score=("mean_pattern_score", "std"),
            median_subject_score=("mean_pattern_score", "median"),
            min_subject_score=("mean_pattern_score", "min"),
            max_subject_score=("mean_pattern_score", "max"),
        )
        .reset_index()
    )
    summary["se_subject_score"] = summary["sd_subject_score"] / np.sqrt(summary["n_subject_groups"].replace(0, np.nan))
    cell_summary = (
        cell_long.groupby(group_cols, dropna=False, observed=True)
        .agg(
            n_cells_cell_level=("pattern_score", "size"),
            mean_cell_score=("pattern_score", "mean"),
            median_cell_score=("pattern_score", "median"),
        )
        .reset_index()
    )
    out = summary.merge(cell_summary, on=group_cols, how="left")
    sort_cols = [column for column in ["pattern", "day_merged_numeric", "broad_cell_type"] if column in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols).reset_index(drop=True)
    return out


def make_cell_long(cells: pd.DataFrame, id_cols: list[str], patterns: list[str]) -> pd.DataFrame:
    long = cells[id_cols + patterns].melt(
        id_vars=id_cols,
        value_vars=patterns,
        var_name="pattern",
        value_name="pattern_score",
    )
    long["pattern"] = pd.Categorical(long["pattern"], categories=patterns, ordered=True)
    return long


def top_gene_lists(top_genes: pd.DataFrame, n_values: tuple[int, ...] = (10, 15, 25)) -> pd.DataFrame:
    require_columns("top genes", top_genes.columns, ["pattern", "rank", "gene", "weight"])
    rows = []
    for pattern, group in top_genes.sort_values(["pattern", "rank"], key=lambda col: col.map(pattern_sort_key) if col.name == "pattern" else col).groupby("pattern", sort=False):
        row = {"pattern": pattern}
        ordered = group.sort_values("rank")
        for n in n_values:
            row[f"top_{n}_genes"] = ", ".join(ordered.head(n)["gene"].astype(str).tolist())
        rows.append(row)
    return pd.DataFrame(rows)


def add_time_values(table: pd.DataFrame, patterns: list[str], value_column: str = "mean_subject_score") -> pd.DataFrame:
    rows = []
    for pattern in patterns:
        group = table[table["pattern"].astype(str) == pattern]
        lookup = {str(row["timepoint_merged"]): float(row[value_column]) for _, row in group.iterrows()}
        rows.append(
            {
                "pattern": pattern,
                "subject_mean_D0": lookup.get("D0", np.nan),
                "subject_mean_D8": lookup.get("D8", np.nan),
                "subject_mean_D10": lookup.get("D10", np.nan),
                "subject_mean_D14_15": lookup.get("D14/15", np.nan),
                "subject_mean_D28": lookup.get("D28", np.nan),
                "D10_minus_D0": lookup.get("D10", np.nan) - lookup.get("D0", np.nan),
                "D14_15_minus_D0": lookup.get("D14/15", np.nan) - lookup.get("D0", np.nan),
            }
        )
    return pd.DataFrame(rows)


def peak_time_table(pattern_by_time: pd.DataFrame, patterns: list[str]) -> pd.DataFrame:
    rows = []
    for pattern in patterns:
        group = pattern_by_time[pattern_by_time["pattern"].astype(str) == pattern].copy()
        if group.empty:
            continue
        group = group.sort_values("mean_subject_score", ascending=False)
        peak = group.iloc[0]
        rows.append(
            {
                "pattern": pattern,
                "subject_peak_timepoint": peak["timepoint_merged"],
                "subject_peak_day": peak["day_merged_numeric"],
                "subject_peak_score": peak["mean_subject_score"],
                "subject_score_range": group["mean_subject_score"].max() - group["mean_subject_score"].min(),
            }
        )
    return pd.DataFrame(rows)


def dominant_cell_type_table(pattern_by_cell_type: pd.DataFrame, patterns: list[str]) -> pd.DataFrame:
    rows = []
    for pattern in patterns:
        group = pattern_by_cell_type[pattern_by_cell_type["pattern"].astype(str) == pattern].copy()
        if group.empty:
            continue
        group = group.sort_values("mean_subject_score", ascending=False)
        top = group.iloc[0]
        rows.append(
            {
                "pattern": pattern,
                "dominant_broad_cell_type": top["broad_cell_type"],
                "dominant_broad_cell_type_mean": top["mean_subject_score"],
                "dominant_broad_cell_type_n_cells": top["total_cells"],
            }
        )
    return pd.DataFrame(rows)


def observed_label(row: pd.Series) -> str:
    if bool(row.get("candidate_ifn_pattern", False)):
        return "IFN-associated temporal activity program"
    pattern_class = str(row.get("pattern_class", ""))
    dominant = str(row.get("dominant_broad_cell_type", ""))
    if pattern_class == "identity-like" and dominant and dominant != "nan":
        return f"{dominant} identity/composition-associated program"
    if pattern_class == "activity-like":
        return "Non-IFN temporal activity-associated program"
    return "Low-specificity pattern"


def interpretation_note(row: pd.Series) -> str:
    if bool(row.get("candidate_ifn_pattern", False)):
        return (
            "Primary IFN-associated K=10 program; high temporal effect, positive IFN score correlation, "
            "and peak usage at the late acute timepoints."
        )
    if str(row.get("pattern_class", "")) == "identity-like":
        return "Interpreted mainly as identity/composition-associated until stratified analyses show otherwise."
    return "Temporal signal requires follow-up before biological naming."


def load_directionality(directionality_dir: Path | None) -> dict[str, pd.DataFrame]:
    if directionality_dir is None:
        return {}
    tables = directionality_dir / "tables"
    pattern_path = tables / "pattern_directionality_summary.csv"
    gene_path = tables / "gene_directionality_summary.csv"
    if not pattern_path.exists() or not gene_path.exists():
        return {}
    return {
        "pattern": pd.read_csv(pattern_path),
        "gene": pd.read_csv(gene_path),
    }


def directionality_for_rq3(candidates: pd.DataFrame, directionality: dict[str, pd.DataFrame]) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    if not directionality:
        candidates = candidates.copy()
        candidates["directionality_source"] = ""
        candidates["D10_median_top_gene_log2fc"] = np.nan
        candidates["D10_frac_top_genes_up"] = np.nan
        candidates["D10_direction_consensus"] = ""
        candidates["D14_15_median_top_gene_log2fc"] = np.nan
        candidates["D14_15_frac_top_genes_up"] = np.nan
        candidates["D14_15_direction_consensus"] = ""
        return candidates

    pattern_direction = directionality["pattern"].copy()
    pattern_direction = pattern_direction[
        (pattern_direction["scope"] == "all_cells")
        & (pattern_direction["stratum"] == "all_cells")
    ]
    rows = []
    for _, row in candidates.iterrows():
        pattern = str(row["pattern"])
        matches = pattern_direction[pattern_direction["pattern"].astype(str) == pattern]
        lookup = {str(r["timepoint_merged"]): r for _, r in matches.iterrows()}
        out = row.to_dict()
        out["directionality_source"] = "all_cells_subject_pseudobulk_top25"
        for timepoint, prefix in [("D10", "D10"), ("D14/15", "D14_15")]:
            match = lookup.get(timepoint)
            if match is None:
                out[f"{prefix}_median_top_gene_log2fc"] = np.nan
                out[f"{prefix}_frac_top_genes_up"] = np.nan
                out[f"{prefix}_direction_consensus"] = ""
            else:
                out[f"{prefix}_median_top_gene_log2fc"] = match["median_gene_log2fc"]
                out[f"{prefix}_frac_top_genes_up"] = match["frac_genes_up"]
                out[f"{prefix}_direction_consensus"] = match["pattern_direction_consensus"]
        rows.append(out)
    return pd.DataFrame(rows)


def write_readme(
    *,
    path: Path,
    model_outdir: Path,
    comparison_outdir: Path,
    k_selection_outdir: Path,
    directionality_dir: Path | None,
    run_stem: str,
    cells: pd.DataFrame,
    pattern_summary: pd.DataFrame,
    rq3: pd.DataFrame,
) -> None:
    ifn_pattern = rq3.iloc[0]["pattern"] if not rq3.empty else "none"
    lines = [
        "# Case Study 5 K=10 Interpretation Layer",
        "",
        "This directory contains additive interpretation tables for the revised K=10 selected model.",
        "The tables summarize CoGAPS pattern usage at the subject level before aggregating across time or cell type.",
        "",
        "## Assumptions",
        "",
        "- The R K=10 selected-model output is used as the primary table source.",
        "- R/Python agreement is tracked separately from the comparison directory.",
        "- Subjects, not cells, are treated as the replicate unit for the interpretation summaries.",
        "- CoGAPS weights remain nonnegative; expression direction is supplied only when the directionality output is available.",
        "",
        "## Inputs",
        "",
        f"- Model output: `{model_outdir}`",
        f"- Run stem: `{run_stem}`",
        f"- R/Python comparison: `{comparison_outdir}`",
        f"- Revised K-selection: `{k_selection_outdir}`",
        f"- Directionality: `{directionality_dir or ''}`",
        "",
        "## Dataset Snapshot",
        "",
        f"- Discovery cells: `{len(cells)}`",
        f"- Subjects: `{cells['subject'].nunique()}`",
        f"- Timepoints: `{cells['timepoint_merged'].nunique()}`",
        f"- Broad cell types: `{cells['broad_cell_type'].nunique()}`",
        f"- IFN candidate pattern: `{ifn_pattern}`",
        "",
        "## Key Output Tables",
        "",
        "- `subject_level/subject_pattern_means.csv`",
        "- `subject_level/subject_time_pattern_means.csv`",
        "- `subject_level/subject_celltype_pattern_means.csv`",
        "- `subject_level/subject_time_celltype_pattern_means.csv`",
        "- `pattern_summaries/pattern_by_time_subject_summary.csv`",
        "- `pattern_summaries/pattern_by_cell_type_subject_summary.csv`",
        "- `pattern_summaries/pattern_by_time_cell_type_subject_summary.csv`",
        "- `rq_tables/RQ1_temporal_programs.csv`",
        "- `rq_tables/RQ2_identity_vs_activity.csv`",
        "- `rq_tables/RQ3_ifn_associated_program.csv`",
        "- `pattern_annotation_table.csv`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build subject-level K=10 interpretation tables for Case Study 5.")
    parser.add_argument("--model-outdir", type=Path, default=DEFAULT_MODEL_OUTDIR)
    parser.add_argument("--comparison-outdir", type=Path, default=DEFAULT_COMPARISON_OUTDIR)
    parser.add_argument("--k-selection-outdir", type=Path, default=DEFAULT_K_SELECTION_OUTDIR)
    parser.add_argument("--directionality-dir", type=Path, default=DEFAULT_DIRECTIONALITY_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--run-stem", default=DEFAULT_RUN_STEM)
    args = parser.parse_args()

    model_outdir = args.model_outdir.resolve()
    comparison_outdir = args.comparison_outdir.resolve()
    k_selection_outdir = args.k_selection_outdir.resolve()
    directionality_dir = args.directionality_dir.resolve() if args.directionality_dir else None
    outdir = args.outdir.resolve()

    cells_path = model_outdir / f"{args.run_stem}.discovery_cells_with_patterns.csv"
    pattern_summary_path = model_outdir / f"{args.run_stem}.pattern_summary.csv"
    top_genes_path = model_outdir / f"{args.run_stem}.top_genes.csv"
    comparison_metrics_path = comparison_outdir / "comparison_metrics.json"

    cells = read_csv(cells_path, "discovery cells")
    pattern_summary = read_csv(pattern_summary_path, "pattern summary")
    top_genes = read_csv(top_genes_path, "top genes")
    require_columns(
        "discovery cells",
        cells.columns,
        ["subject", "cohort", "timepoint_merged", "day_merged_numeric", "broad_cell_type"],
    )
    require_columns(
        "pattern summary",
        pattern_summary.columns,
        [
            "pattern",
            "eta_timepoint",
            "eta_broad_cell_type",
            "pattern_class",
            "spearman_ifn_score_rho",
            "peak_timepoint",
            "ifn_top_gene_overlap_top15",
            "candidate_ifn_pattern",
        ],
    )

    patterns = pattern_columns(cells)
    pattern_summary = pattern_summary.copy()
    pattern_summary["candidate_ifn_pattern"] = bool_series(pattern_summary["candidate_ifn_pattern"])
    pattern_summary["pattern"] = pattern_summary["pattern"].astype(str)
    top_lists = top_gene_lists(top_genes)

    extra_mean_cols = ["ifn_program_score", "translation_program_score", "mito_program_score", "plasmablast_program_score"]
    subject = group_pattern_means(cells, ["subject"], patterns, extra_mean_cols)
    subject_time = group_pattern_means(cells, ["subject", "timepoint_merged", "day_merged_numeric"], patterns, extra_mean_cols)
    subject_celltype = group_pattern_means(cells, ["subject", "broad_cell_type"], patterns, extra_mean_cols)
    subject_time_celltype = group_pattern_means(
        cells,
        ["subject", "broad_cell_type", "timepoint_merged", "day_merged_numeric"],
        patterns,
        extra_mean_cols,
    )

    subject_time_long = melt_patterns(
        subject_time,
        ["subject", "timepoint_merged", "day_merged_numeric", "n_cells"],
        patterns,
    )
    subject_celltype_long = melt_patterns(
        subject_celltype,
        ["subject", "broad_cell_type", "n_cells"],
        patterns,
    )
    subject_time_celltype_long = melt_patterns(
        subject_time_celltype,
        ["subject", "broad_cell_type", "timepoint_merged", "day_merged_numeric", "n_cells"],
        patterns,
    )
    cell_time_long = make_cell_long(cells, ["timepoint_merged", "day_merged_numeric"], patterns)
    cell_celltype_long = make_cell_long(cells, ["broad_cell_type"], patterns)
    cell_time_celltype_long = make_cell_long(cells, ["broad_cell_type", "timepoint_merged", "day_merged_numeric"], patterns)

    pattern_by_time = subject_summary(
        subject_time_long,
        cell_time_long,
        ["pattern", "timepoint_merged", "day_merged_numeric"],
    )
    pattern_by_cell_type = subject_summary(
        subject_celltype_long,
        cell_celltype_long,
        ["pattern", "broad_cell_type"],
    )
    pattern_by_time_cell_type = subject_summary(
        subject_time_celltype_long,
        cell_time_celltype_long,
        ["pattern", "broad_cell_type", "timepoint_merged", "day_merged_numeric"],
    )

    peaks = peak_time_table(pattern_by_time, patterns)
    time_values = add_time_values(pattern_by_time, patterns)
    dominant = dominant_cell_type_table(pattern_by_cell_type, patterns)

    annotation = (
        pattern_summary.merge(peaks, on="pattern", how="left")
        .merge(time_values, on="pattern", how="left")
        .merge(dominant, on="pattern", how="left")
        .merge(top_lists, on="pattern", how="left")
    )
    annotation["eta_time_to_cell_ratio"] = annotation["eta_timepoint"] / annotation["eta_broad_cell_type"].replace(0, np.nan)
    annotation["observed_label"] = annotation.apply(observed_label, axis=1)
    annotation["interpretation_note"] = annotation.apply(interpretation_note, axis=1)
    annotation = annotation.sort_values("pattern", key=lambda col: col.map(pattern_sort_key)).reset_index(drop=True)

    rq1 = annotation[
        [
            "pattern",
            "observed_label",
            "pattern_class",
            "eta_timepoint",
            "kruskal_time_p_adj",
            "subject_peak_timepoint",
            "subject_peak_score",
            "subject_score_range",
            "subject_mean_D0",
            "subject_mean_D8",
            "subject_mean_D10",
            "subject_mean_D14_15",
            "subject_mean_D28",
            "D10_minus_D0",
            "D14_15_minus_D0",
            "spearman_ifn_score_rho",
            "candidate_ifn_pattern",
            "top_10_genes",
            "interpretation_note",
        ]
    ].sort_values(["eta_timepoint", "subject_score_range"], ascending=False)

    rq2 = annotation[
        [
            "pattern",
            "observed_label",
            "pattern_class",
            "eta_timepoint",
            "eta_broad_cell_type",
            "eta_time_to_cell_ratio",
            "dominant_broad_cell_type",
            "dominant_broad_cell_type_mean",
            "subject_peak_timepoint",
            "candidate_ifn_pattern",
            "top_10_genes",
            "interpretation_note",
        ]
    ].sort_values("pattern", key=lambda col: col.map(pattern_sort_key))

    rq3_base = annotation[annotation["candidate_ifn_pattern"]].copy()
    if rq3_base.empty:
        rq3_base = annotation.sort_values("spearman_ifn_score_rho", ascending=False).head(1).copy()
    rq3 = rq3_base[
        [
            "pattern",
            "observed_label",
            "eta_timepoint",
            "eta_broad_cell_type",
            "spearman_ifn_score_rho",
            "ifn_score_p_adj",
            "ifn_top_gene_overlap_top15",
            "peak_timepoint",
            "subject_peak_timepoint",
            "subject_mean_D0",
            "subject_mean_D8",
            "subject_mean_D10",
            "subject_mean_D14_15",
            "subject_mean_D28",
            "D10_minus_D0",
            "D14_15_minus_D0",
            "dominant_broad_cell_type",
            "top_15_genes",
            "top_25_genes",
            "interpretation_note",
        ]
    ].copy()
    rq3 = directionality_for_rq3(rq3, load_directionality(directionality_dir))

    outdir.mkdir(parents=True, exist_ok=True)
    subject_dir = outdir / "subject_level"
    summary_dir = outdir / "pattern_summaries"
    rq_dir = outdir / "rq_tables"
    subject_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)
    rq_dir.mkdir(parents=True, exist_ok=True)

    subject.to_csv(subject_dir / "subject_pattern_means.csv", index=False)
    subject_time.to_csv(subject_dir / "subject_time_pattern_means.csv", index=False)
    subject_celltype.to_csv(subject_dir / "subject_celltype_pattern_means.csv", index=False)
    subject_time_celltype.to_csv(subject_dir / "subject_time_celltype_pattern_means.csv", index=False)
    subject_time_long.to_csv(subject_dir / "subject_time_pattern_means_long.csv", index=False)
    subject_celltype_long.to_csv(subject_dir / "subject_celltype_pattern_means_long.csv", index=False)
    subject_time_celltype_long.to_csv(subject_dir / "subject_time_celltype_pattern_means_long.csv", index=False)

    pattern_by_time.to_csv(summary_dir / "pattern_by_time_subject_summary.csv", index=False)
    pattern_by_cell_type.to_csv(summary_dir / "pattern_by_cell_type_subject_summary.csv", index=False)
    pattern_by_time_cell_type.to_csv(summary_dir / "pattern_by_time_cell_type_subject_summary.csv", index=False)

    rq1.to_csv(rq_dir / "RQ1_temporal_programs.csv", index=False)
    rq2.to_csv(rq_dir / "RQ2_identity_vs_activity.csv", index=False)
    rq3.to_csv(rq_dir / "RQ3_ifn_associated_program.csv", index=False)
    annotation.to_csv(outdir / "pattern_annotation_table.csv", index=False)

    comparison_metrics = {}
    if comparison_metrics_path.exists():
        comparison_metrics = json.loads(comparison_metrics_path.read_text(encoding="utf-8"))
    manifest = {
        "created_utc": utc_now(),
        "model_outdir": str(model_outdir),
        "comparison_outdir": str(comparison_outdir),
        "k_selection_outdir": str(k_selection_outdir),
        "directionality_dir": str(directionality_dir) if directionality_dir else "",
        "run_stem": args.run_stem,
        "n_cells": int(len(cells)),
        "n_subjects": int(cells["subject"].nunique()),
        "n_timepoints": int(cells["timepoint_merged"].nunique()),
        "n_broad_cell_types": int(cells["broad_cell_type"].nunique()),
        "patterns": patterns,
        "ifn_candidate_patterns": rq3["pattern"].astype(str).tolist(),
        "r_python_agreement": {
            "mean_matched_gene_loading_spearman": comparison_metrics.get("mean_matched_gene_loading_spearman"),
            "mean_matched_cell_score_spearman": comparison_metrics.get("mean_matched_cell_score_spearman"),
            "mean_matched_top_50_gene_jaccard": comparison_metrics.get("mean_matched_top_50_gene_jaccard"),
        },
        "outputs": {
            "subject_level_dir": str(subject_dir),
            "pattern_summaries_dir": str(summary_dir),
            "rq_tables_dir": str(rq_dir),
            "pattern_annotation_table": str(outdir / "pattern_annotation_table.csv"),
        },
    }
    (outdir / "interpretation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    write_readme(
        path=outdir / "README_interpretation_layer.md",
        model_outdir=model_outdir,
        comparison_outdir=comparison_outdir,
        k_selection_outdir=k_selection_outdir,
        directionality_dir=directionality_dir,
        run_stem=args.run_stem,
        cells=cells,
        pattern_summary=pattern_summary,
        rq3=rq3,
    )

    print(f"[done] wrote K=10 interpretation layer to {outdir}")
    print(f"[done] IFN-associated pattern(s): {', '.join(rq3['pattern'].astype(str).tolist())}")
    print(f"[done] subject-level rows: time={len(subject_time)}, celltype={len(subject_celltype)}, time_celltype={len(subject_time_celltype)}")


if __name__ == "__main__":
    main()
