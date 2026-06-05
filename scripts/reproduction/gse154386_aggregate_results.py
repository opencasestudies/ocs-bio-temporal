#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from gse154386_cogaps_sweep_common import (
    DEFAULT_OUTDIR,
    cache_paths,
    parse_int_csv,
    read_json,
    resolve_path,
    utc_now_iso,
    write_json,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate GSE154386 CoGAPS sweep metrics, score reproducibility and redundancy, "
            "and choose K using stability-first logic."
        )
    )
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Sweep output directory.")
    parser.add_argument(
        "--k-grid",
        default=None,
        help="Optional expected K grid used for coverage warnings only, for example 6,8,10,12,14.",
    )
    parser.add_argument(
        "--seeds",
        default=None,
        help="Optional expected seeds used for coverage warnings only, for example 1,2,3,4,5.",
    )
    parser.add_argument(
        "--iters",
        default=None,
        help="Optional expected iteration counts used for coverage warnings only, for example 4000,10000.",
    )
    parser.add_argument(
        "--top-genes",
        type=int,
        default=50,
        help="Top genes per pattern used for matched top-gene overlap calculations.",
    )
    parser.add_argument(
        "--plateau-tol",
        type=float,
        default=0.03,
        help="Absolute tolerance on the stability score when defining the stability plateau.",
    )
    parser.add_argument(
        "--min-successful-runs",
        type=int,
        default=2,
        help="Minimum successful seeds required before a (K, n_iter) config is eligible for selection.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation.",
    )
    return parser


def metrics_records_to_dataframe(records):
    import pandas as pd

    rows = []
    for record in records:
        row = {}
        for key, value in record.items():
            if isinstance(value, (dict, list)):
                row[key] = json.dumps(value, sort_keys=True)
            else:
                row[key] = value
        rows.append(row)
    return pd.DataFrame(rows)


def load_gene_loadings(path, pd):
    frame = pd.read_csv(path, index_col=0)
    pattern_cols = [col for col in frame.columns if str(col).startswith("Pattern")]
    if not pattern_cols:
        raise ValueError(f"No Pattern columns found in gene loading CSV: {path}")
    frame = frame[pattern_cols].copy()
    frame.index = frame.index.astype(str)
    return frame


def load_pattern_summary(path, pd):
    frame = pd.read_csv(path)
    if "pattern" not in frame.columns:
        raise ValueError(f"pattern column missing from pattern summary CSV: {path}")
    frame["pattern"] = frame["pattern"].astype(str)
    return frame.set_index("pattern", drop=False)


def top_gene_sets(gene_loadings_df, top_n):
    out = {}
    for pattern in gene_loadings_df.columns:
        top = gene_loadings_df[pattern].sort_values(ascending=False).head(top_n)
        out[str(pattern)] = set(top.index.astype(str))
    return out


def greedy_assignment(similarity_matrix):
    flat = []
    for i in range(similarity_matrix.shape[0]):
        for j in range(similarity_matrix.shape[1]):
            flat.append((float(similarity_matrix[i, j]), i, j))
    flat.sort(reverse=True, key=lambda row: row[0])
    used_i = set()
    used_j = set()
    matches = []
    for _score, i, j in flat:
        if i in used_i or j in used_j:
            continue
        used_i.add(i)
        used_j.add(j)
        matches.append((i, j))
    return matches


def optimal_assignment(similarity_matrix):
    try:
        from scipy.optimize import linear_sum_assignment
    except Exception:
        return greedy_assignment(similarity_matrix)

    row_ind, col_ind = linear_sum_assignment(-similarity_matrix)
    return list(zip(row_ind.tolist(), col_ind.tolist()))


def safe_spearman(left, right):
    corr = left.corr(right, method="spearman")
    if corr is None:
        return 0.0
    try:
        value = float(corr)
    except Exception:
        return 0.0
    if value != value:
        return 0.0
    if value == float("inf") or value == float("-inf"):
        return 0.0
    return value


def pairwise_reproducibility(row_a, row_b, pd, np, top_n, table_cache):
    key_a = str(row_a["metrics_path"])
    key_b = str(row_b["metrics_path"])
    if key_a not in table_cache:
        table_cache[key_a] = (
            load_gene_loadings(Path(row_a["gene_loadings_csv"]), pd),
            load_pattern_summary(Path(row_a["pattern_summary_csv"]), pd),
        )
    if key_b not in table_cache:
        table_cache[key_b] = (
            load_gene_loadings(Path(row_b["gene_loadings_csv"]), pd),
            load_pattern_summary(Path(row_b["pattern_summary_csv"]), pd),
        )

    gene_a, summary_a = table_cache[key_a]
    gene_b, summary_b = table_cache[key_b]

    common_genes = gene_a.index.intersection(gene_b.index)
    if len(common_genes) == 0:
        raise ValueError(
            f"No common genes between {row_a['gene_loadings_csv']} and {row_b['gene_loadings_csv']}"
        )

    gene_a = gene_a.loc[common_genes]
    gene_b = gene_b.loc[common_genes]
    patterns_a = list(gene_a.columns)
    patterns_b = list(gene_b.columns)
    set_a = top_gene_sets(gene_a, top_n=top_n)
    set_b = top_gene_sets(gene_b, top_n=top_n)

    similarity = np.zeros((len(patterns_a), len(patterns_b)), dtype=float)
    jaccard = np.zeros((len(patterns_a), len(patterns_b)), dtype=float)
    for i, pattern_a in enumerate(patterns_a):
        left = gene_a[pattern_a]
        for j, pattern_b in enumerate(patterns_b):
            right = gene_b[pattern_b]
            similarity[i, j] = safe_spearman(left, right)
            union = set_a[str(pattern_a)].union(set_b[str(pattern_b)])
            jaccard[i, j] = (
                float(len(set_a[str(pattern_a)].intersection(set_b[str(pattern_b)])) / len(union))
                if union else 0.0
            )

    matched_pairs = optimal_assignment(similarity)
    matched_corr_raw = []
    matched_corr_clipped = []
    matched_jaccard = []
    peak_agreements = []
    class_agreements = []
    for i, j in matched_pairs:
        pattern_a = str(patterns_a[i])
        pattern_b = str(patterns_b[j])
        raw_corr = float(similarity[i, j])
        matched_corr_raw.append(raw_corr)
        matched_corr_clipped.append(max(0.0, raw_corr))
        matched_jaccard.append(float(jaccard[i, j]))

        peak_a = summary_a.loc[pattern_a, "peak_timepoint"] if pattern_a in summary_a.index else None
        peak_b = summary_b.loc[pattern_b, "peak_timepoint"] if pattern_b in summary_b.index else None
        if pd.notna(peak_a) and pd.notna(peak_b):
            peak_agreements.append(float(str(peak_a) == str(peak_b)))

        class_a = summary_a.loc[pattern_a, "pattern_class"] if pattern_a in summary_a.index else None
        class_b = summary_b.loc[pattern_b, "pattern_class"] if pattern_b in summary_b.index else None
        if pd.notna(class_a) and pd.notna(class_b):
            class_agreements.append(float(str(class_a) == str(class_b)))

    return {
        "K": int(row_a["K"]),
        "n_iter": int(row_a["n_iter"]),
        "seed_a": int(row_a["seed"]),
        "seed_b": int(row_b["seed"]),
        "matched_pattern_count": len(matched_pairs),
        "matched_loading_corr_mean_raw": float(np.mean(matched_corr_raw)) if matched_corr_raw else 0.0,
        "matched_loading_corr_mean_clipped": float(np.mean(matched_corr_clipped)) if matched_corr_clipped else 0.0,
        "matched_loading_corr_median_raw": float(np.median(matched_corr_raw)) if matched_corr_raw else 0.0,
        "matched_top_gene_jaccard_mean": float(np.mean(matched_jaccard)) if matched_jaccard else 0.0,
        "matched_peak_timepoint_agreement": float(np.mean(peak_agreements)) if peak_agreements else 0.0,
        "matched_pattern_class_agreement": float(np.mean(class_agreements)) if class_agreements else 0.0,
    }


def coverage_warnings(ok_df, expected_ks, expected_seeds, expected_iters):
    warnings = []
    if expected_ks is None or expected_seeds is None or expected_iters is None:
        return warnings

    expected = {
        (int(k_value), int(seed), int(n_iter))
        for k_value in expected_ks
        for seed in expected_seeds
        for n_iter in expected_iters
    }
    observed = {
        (int(row.K), int(row.seed), int(row.n_iter))
        for row in ok_df.itertuples(index=False)
    }
    missing = sorted(expected.difference(observed))
    if missing:
        preview = ", ".join(str(entry) for entry in missing[:10])
        suffix = "" if len(missing) <= 10 else f" ... and {len(missing) - 10} more"
        warnings.append(f"Missing successful runs for {len(missing)} expected configs: {preview}{suffix}")
    return warnings


def summarize_config(group, pair_df, np):
    k_value = int(group["K"].iloc[0])
    n_iter = int(group["n_iter"].iloc[0])
    successful_seed_count = int(group["seed"].nunique())
    pair_subset = pair_df[(pair_df["K"] == k_value) & (pair_df["n_iter"] == n_iter)].copy()

    effective_fraction = (group["n_patterns_effective"].astype(float) / float(k_value)).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    ifn_secondary = (
        (
            group["candidate_ifn_pattern_count"].astype(float).clip(lower=0.0) / 1.0
        ).clip(upper=1.0) +
        (
            group["max_ifn_top_gene_overlap_top15"].astype(float).clip(lower=0.0) / 5.0
        ).clip(upper=1.0)
    ) / 2.0

    mean_corr = float(pair_subset["matched_loading_corr_mean_clipped"].mean()) if len(pair_subset) else 0.0
    mean_jaccard = float(pair_subset["matched_top_gene_jaccard_mean"].mean()) if len(pair_subset) else 0.0
    stability_core = (0.7 * mean_corr) + (0.3 * mean_jaccard)

    redundancy_penalty = (
        0.50 * float(group["within_run_pattern_redundancy_mean"].astype(float).mean()) +
        0.20 * float(group["within_run_pattern_redundancy_max"].astype(float).mean()) +
        0.20 * float(group["within_run_top_gene_jaccard_mean"].astype(float).mean()) +
        0.10 * float((1.0 - effective_fraction).mean())
    )

    secondary_score = (
        0.50 * float(group["pattern_mix_balance"].astype(float).mean()) +
        0.35 * (float(pair_subset["matched_peak_timepoint_agreement"].mean()) if len(pair_subset) else 0.0) +
        0.15 * float(ifn_secondary.mean())
    )

    overall_score = stability_core - (0.35 * redundancy_penalty) + (0.08 * secondary_score)

    return {
        "K": k_value,
        "n_iter": n_iter,
        "successful_seed_count": successful_seed_count,
        "pairwise_seed_comparisons": int(len(pair_subset)),
        "mean_matched_loading_corr": mean_corr,
        "mean_matched_top_gene_jaccard": mean_jaccard,
        "peak_timepoint_agreement_mean": float(pair_subset["matched_peak_timepoint_agreement"].mean()) if len(pair_subset) else 0.0,
        "pattern_class_agreement_mean": float(pair_subset["matched_pattern_class_agreement"].mean()) if len(pair_subset) else 0.0,
        "effective_fraction_mean": float(effective_fraction.mean()),
        "degenerate_fraction_mean": float((1.0 - effective_fraction).mean()),
        "within_run_pattern_redundancy_mean": float(group["within_run_pattern_redundancy_mean"].astype(float).mean()),
        "within_run_pattern_redundancy_max": float(group["within_run_pattern_redundancy_max"].astype(float).mean()),
        "within_run_top_gene_jaccard_mean": float(group["within_run_top_gene_jaccard_mean"].astype(float).mean()),
        "pattern_mix_balance_mean": float(group["pattern_mix_balance"].astype(float).mean()),
        "mean_eta_timepoint": float(group["mean_eta_timepoint"].astype(float).mean()),
        "max_eta_timepoint": float(group["max_eta_timepoint"].astype(float).max()),
        "candidate_ifn_pattern_count_mean": float(group["candidate_ifn_pattern_count"].astype(float).mean()),
        "ifn_secondary_score": float(ifn_secondary.mean()),
        "stability_core": float(stability_core),
        "redundancy_penalty": float(redundancy_penalty),
        "secondary_score": float(secondary_score),
        "overall_score": float(overall_score),
    }


def pick_best_per_k(summary_config_df):
    best_rows = (
        summary_config_df
        .sort_values(
            ["selection_eligible", "overall_score", "successful_seed_count", "stability_core", "redundancy_penalty", "n_iter"],
            ascending=[False, False, False, False, True, True],
        )
        .drop_duplicates(subset=["K"], keep="first")
        .copy()
        .sort_values("K")
    )
    best_rows = best_rows.rename(columns={"n_iter": "selected_n_iter"})
    return best_rows


def generate_plots(summary_config_df, best_by_k_df, figures_dir, plateau_threshold):
    import matplotlib.pyplot as plt

    figures_dir.mkdir(parents=True, exist_ok=True)

    if not best_by_k_df.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(best_by_k_df["K"], best_by_k_df["stability_core"], marker="o", label="stability_core")
        ax.plot(best_by_k_df["K"], best_by_k_df["overall_score"], marker="s", label="overall_score")
        ax.axhline(plateau_threshold, color="tab:gray", linestyle="--", linewidth=1.2, label="plateau threshold")
        ax.set_xlabel("K")
        ax.set_ylabel("score")
        ax.set_title("Best-per-K stability and overall scores")
        ax.legend()
        fig.tight_layout()
        fig.savefig(figures_dir / "best_config_scores_by_k.png", dpi=160)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(best_by_k_df["redundancy_penalty"], best_by_k_df["stability_core"], s=70)
        for row in best_by_k_df.itertuples(index=False):
            ax.annotate(f"K={row.K}\niter={row.selected_n_iter}", (row.redundancy_penalty, row.stability_core))
        ax.set_xlabel("redundancy_penalty")
        ax.set_ylabel("stability_core")
        ax.set_title("Best-per-K stability vs redundancy")
        fig.tight_layout()
        fig.savefig(figures_dir / "best_config_stability_vs_redundancy.png", dpi=160)
        plt.close(fig)

    if not summary_config_df.empty:
        fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(summary_config_df)), 5))
        labels = [f"K={int(row.K)}\niter={int(row.n_iter)}" for row in summary_config_df.itertuples(index=False)]
        ax.bar(range(len(summary_config_df)), summary_config_df["successful_seed_count"])
        ax.set_xticks(range(len(summary_config_df)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylabel("successful seeds")
        ax.set_title("Completed successful runs by (K, n_iter)")
        fig.tight_layout()
        fig.savefig(figures_dir / "successful_runs_by_config.png", dpi=160)
        plt.close(fig)


def write_report(report_path, source_script_path, prep_manifest_path, coverage_messages, per_run_df, summary_config_df, best_by_k_df, chosen_row, plateau_threshold):
    lines = [
        "# GSE154386 CoGAPS sweep report",
        "",
        f"- Generated at: `{utc_now_iso()}`",
        f"- Source Case Study 5 script: `{source_script_path}`",
        f"- Prep manifest: `{prep_manifest_path}`",
        "",
        "## Coverage",
        f"- Metrics JSON files discovered: {len(per_run_df)}",
        f"- Successful runs: {int((per_run_df['status'].astype(str) == 'ok').sum())}",
    ]
    if coverage_messages:
        lines.append("- Coverage warnings:")
        for message in coverage_messages:
            lines.append(f"  - {message}")
    else:
        lines.append("- Coverage warnings: none")

    lines.extend(
        [
            "",
            "## Selection logic",
            "- Primary score is stability across seeds: matched gene-loading correlations plus matched top-gene overlap.",
            "- Redundancy penalties come from within-run pattern similarity and effectively empty/degenerate patterns.",
            "- Secondary signals reward a useful mix of identity-like and activity-like patterns and stable time-varying behavior across seeds.",
            "- IFN-like signal is recorded, but only as a low-weight secondary metric.",
            f"- Stability plateau threshold score: `{plateau_threshold:.3f}`. K selection then prefers the smallest K on or above that plateau.",
            "",
            "## Chosen K",
            f"- Chosen K: `{int(chosen_row['K'])}`",
            f"- Representative n_iter for that K: `{int(chosen_row['selected_n_iter'])}`",
            f"- stability_core: `{chosen_row['stability_core']:.4f}`",
            f"- redundancy_penalty: `{chosen_row['redundancy_penalty']:.4f}`",
            f"- overall_score: `{chosen_row['overall_score']:.4f}`",
            "",
            "## Top candidate K values",
            "```text",
            best_by_k_df[
                [
                    "K",
                    "selected_n_iter",
                    "successful_seed_count",
                    "stability_core",
                    "redundancy_penalty",
                    "overall_score",
                    "on_stability_plateau",
                    "chosen_K",
                ]
            ].to_string(index=False),
            "```",
            "",
            "## Summary by (K, n_iter)",
            "```text",
            summary_config_df[
                [
                    "K",
                    "n_iter",
                    "successful_seed_count",
                    "pairwise_seed_comparisons",
                    "stability_core",
                    "redundancy_penalty",
                    "overall_score",
                    "selection_eligible",
                ]
            ].to_string(index=False),
            "```",
            "",
            "## Interpretation",
            "- Prefer the smallest K whose seed-to-seed solutions are already stable.",
            "- Reject larger K values when they mostly add redundant or weak patterns instead of new reproducible structure.",
            "- Use the chosen K as the default, but review nearby plateau candidates if biological interpretability or runtime tradeoffs matter.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = build_parser().parse_args()
    if args.top_genes <= 0:
        raise ValueError("--top-genes must be positive")
    if args.plateau_tol < 0:
        raise ValueError("--plateau-tol must be non-negative")
    if args.min_successful_runs < 2:
        raise ValueError("--min-successful-runs must be at least 2")

    outdir = resolve_path(args.outdir)
    cache = cache_paths(outdir)
    runs_dir = outdir / "runs"
    if not runs_dir.exists():
        raise FileNotFoundError(f"Runs directory does not exist: {runs_dir}")

    import numpy as np  # lazy imports so --help works without cluster deps
    import pandas as pd

    metric_paths = sorted(runs_dir.glob("*.metrics.json"))
    if not metric_paths:
        raise FileNotFoundError(f"No metrics JSON files were found under {runs_dir}")

    raw_records = []
    for path in metric_paths:
        record = read_json(path)
        record.setdefault("metrics_path", str(path.resolve()))
        raw_records.append(record)

    per_run_df = metrics_records_to_dataframe(raw_records)
    if per_run_df.empty:
        raise ValueError("No metrics records could be parsed.")

    per_run_csv = outdir / "per_run_metrics.csv"
    per_run_df.to_csv(per_run_csv, index=False)

    ok_df = per_run_df.loc[per_run_df["status"].astype(str) == "ok"].copy()
    if ok_df.empty:
        raise ValueError(
            f"Found {len(per_run_df)} metrics JSON files but none with status=ok. "
            "Check the per-run metrics for actionable errors."
        )

    for column in (
        "K",
        "seed",
        "n_iter",
        "n_patterns_effective",
        "within_run_pattern_redundancy_mean",
        "within_run_pattern_redundancy_max",
        "within_run_top_gene_jaccard_mean",
        "pattern_mix_balance",
        "candidate_ifn_pattern_count",
        "max_ifn_top_gene_overlap_top15",
        "mean_eta_timepoint",
        "max_eta_timepoint",
    ):
        ok_df[column] = pd.to_numeric(ok_df[column], errors="coerce")

    expected_ks = parse_int_csv(args.k_grid) if args.k_grid else None
    expected_seeds = parse_int_csv(args.seeds) if args.seeds else None
    expected_iters = parse_int_csv(args.iters) if args.iters else None
    warnings = coverage_warnings(ok_df=ok_df, expected_ks=expected_ks, expected_seeds=expected_seeds, expected_iters=expected_iters)

    pair_rows = []
    table_cache = {}
    for (k_value, n_iter), group in ok_df.groupby(["K", "n_iter"], observed=True):
        group = group.sort_values("seed")
        rows = list(group.to_dict(orient="records"))
        for row_a, row_b in itertools.combinations(rows, 2):
            pair_rows.append(
                pairwise_reproducibility(
                    row_a=row_a,
                    row_b=row_b,
                    pd=pd,
                    np=np,
                    top_n=args.top_genes,
                    table_cache=table_cache,
                )
            )
    pair_df = pd.DataFrame(
        pair_rows,
        columns=[
            "K",
            "n_iter",
            "seed_a",
            "seed_b",
            "matched_pattern_count",
            "matched_loading_corr_mean_raw",
            "matched_loading_corr_mean_clipped",
            "matched_loading_corr_median_raw",
            "matched_top_gene_jaccard_mean",
            "matched_peak_timepoint_agreement",
            "matched_pattern_class_agreement",
        ],
    )
    pair_csv = outdir / "pairwise_reproducibility.csv"
    pair_df.to_csv(pair_csv, index=False)

    summary_rows = []
    for (_k_value, _n_iter), group in ok_df.groupby(["K", "n_iter"], observed=True):
        summary = summarize_config(group=group, pair_df=pair_df, np=np)
        summary["selection_eligible"] = bool(
            summary["successful_seed_count"] >= args.min_successful_runs and summary["pairwise_seed_comparisons"] > 0
        )
        summary_rows.append(summary)

    summary_config_df = pd.DataFrame(summary_rows).sort_values(["K", "n_iter"]).reset_index(drop=True)
    if summary_config_df.empty:
        raise ValueError("No successful (K, n_iter) configurations were available to summarize.")

    summary_config_csv = outdir / "summary_by_k_n_iter.csv"
    summary_config_df.to_csv(summary_config_csv, index=False)

    best_by_k_df = pick_best_per_k(summary_config_df)
    eligible_best = best_by_k_df.loc[best_by_k_df["selection_eligible"]].copy()
    if eligible_best.empty:
        raise ValueError(
            "No K values are eligible for selection. "
            "At least two successful seeds are required for one (K, n_iter) configuration."
        )

    max_stability = float(eligible_best["stability_core"].max())
    plateau_threshold = max_stability - float(args.plateau_tol)
    best_by_k_df["on_stability_plateau"] = best_by_k_df["stability_core"] >= plateau_threshold

    plateau_candidates = best_by_k_df.loc[
        best_by_k_df["selection_eligible"] & best_by_k_df["on_stability_plateau"]
    ].copy()
    if plateau_candidates.empty:
        plateau_candidates = eligible_best.sort_values("stability_core", ascending=False).head(1).copy()

    chosen_row = (
        plateau_candidates
        .sort_values(
            ["K", "redundancy_penalty", "degenerate_fraction_mean", "overall_score", "selected_n_iter"],
            ascending=[True, True, True, False, True],
        )
        .iloc[0]
    )
    chosen_k = int(chosen_row["K"])
    best_by_k_df["chosen_K"] = best_by_k_df["K"].astype(int) == chosen_k
    summary_by_k_csv = outdir / "summary_by_K.csv"
    best_by_k_df.to_csv(summary_by_k_csv, index=False)

    selection_summary_json = outdir / "chosen_k_summary.json"
    write_json(
        selection_summary_json,
        {
            "generated_at_utc": utc_now_iso(),
            "chosen_K": chosen_k,
            "selected_n_iter": int(chosen_row["selected_n_iter"]),
            "plateau_threshold": plateau_threshold,
            "source_metrics_csv": str(per_run_csv),
            "source_summary_by_k_csv": str(summary_by_k_csv),
        },
    )

    if not args.no_plots:
        generate_plots(
            summary_config_df=summary_config_df,
            best_by_k_df=best_by_k_df,
            figures_dir=outdir / "figures",
            plateau_threshold=plateau_threshold,
        )

    prep_manifest = cache["prep_manifest_json"]
    source_script_path = read_json(prep_manifest).get("source_script_path", "unknown") if prep_manifest.exists() else "unknown"
    report_path = outdir / "report.md"
    write_report(
        report_path=report_path,
        source_script_path=source_script_path,
        prep_manifest_path=str(prep_manifest),
        coverage_messages=warnings,
        per_run_df=per_run_df,
        summary_config_df=summary_config_df,
        best_by_k_df=best_by_k_df.sort_values(["on_stability_plateau", "overall_score"], ascending=[False, False]).reset_index(drop=True),
        chosen_row=chosen_row,
        plateau_threshold=plateau_threshold,
    )

    print(f"[aggregate] wrote {per_run_csv}")
    print(f"[aggregate] wrote {pair_csv}")
    print(f"[aggregate] wrote {summary_config_csv}")
    print(f"[aggregate] wrote {summary_by_k_csv}")
    print(f"[aggregate] wrote {selection_summary_json}")
    print(f"[aggregate] wrote {report_path}")
    if not args.no_plots:
        print(f"[aggregate] wrote figures under {outdir / 'figures'}")
    print(f"[aggregate] chosen K={chosen_k} using selected_n_iter={int(chosen_row['selected_n_iter'])}")


if __name__ == "__main__":
    main()
