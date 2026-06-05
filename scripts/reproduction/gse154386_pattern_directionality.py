#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_SWEEP_OUTDIR = Path("data/cs5_sweep_results")
DEFAULT_RUN_STEM = "cogaps_K10_seed2_iter2000"
DEFAULT_OUTDIR = Path("GSE154386/pattern_directionality_revised_model_K10_seed2_iter2000")


@dataclass(frozen=True)
class GeneSpec:
    pattern: str
    rank: int
    gene: str
    weight: float


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate whether selected CoGAPS pattern genes are up- or down-regulated "
            "relative to D0 using pseudobulk expression contrasts."
        )
    )
    parser.add_argument(
        "--sweep-outdir",
        type=Path,
        default=DEFAULT_SWEEP_OUTDIR,
        help="Sweep output directory containing cache/ and runs/. Defaults to the synced local results.",
    )
    parser.add_argument(
        "--run-stem",
        default=DEFAULT_RUN_STEM,
        help="Selected run stem used to locate top-gene and pattern summary CSVs.",
    )
    parser.add_argument(
        "--expression-h5ad",
        type=Path,
        default=None,
        help="Expression AnnData to use. Defaults to cache/gse154386_preprocessed_hvg.h5ad.",
    )
    parser.add_argument(
        "--top-genes-csv",
        type=Path,
        default=None,
        help="Selected run top genes CSV. Defaults to runs/<run-stem>.top_genes.csv.",
    )
    parser.add_argument(
        "--pattern-summary-csv",
        type=Path,
        default=None,
        help="Selected run pattern summary CSV. Defaults to runs/<run-stem>.pattern_summary.csv.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory for directionality tables and figures.",
    )
    parser.add_argument("--top-n", type=int, default=25, help="Top genes per pattern to test.")
    parser.add_argument(
        "--baseline-timepoint",
        default="D0",
        help="Experimental baseline timepoint for directionality contrasts.",
    )
    parser.add_argument(
        "--cohort",
        default="experimental",
        help="Cohort to analyze from the expression object.",
    )
    parser.add_argument(
        "--layer",
        default="counts",
        help="AnnData layer used for pseudobulk expression. Use 'X' for adata.X.",
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        default=10,
        help="Minimum cells required in both baseline and contrast pseudobulk groups.",
    )
    parser.add_argument(
        "--log2fc-threshold",
        type=float,
        default=0.25,
        help="Absolute log2FC threshold used for up/down/near-zero labels.",
    )
    parser.add_argument(
        "--cpm-pseudocount",
        type=float,
        default=1.0,
        help="Pseudocount added to CPM before log2 transformation.",
    )
    parser.add_argument(
        "--heatmap-abs-cap",
        type=float,
        default=1.5,
        help="Absolute log2FC cap for diverging heatmap colors.",
    )
    return parser


def require_columns(table_name: str, columns: Iterable[str], required: Iterable[str]) -> None:
    missing = [column for column in required if column not in columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {', '.join(missing)}")


def read_top_genes(path: Path, top_n: int) -> list[GeneSpec]:
    rows: list[GeneSpec] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns("top genes CSV", reader.fieldnames or [], ["pattern", "rank", "gene", "weight"])
        for row in reader:
            rank = int(row["rank"])
            if rank <= top_n:
                rows.append(
                    GeneSpec(
                        pattern=str(row["pattern"]),
                        rank=rank,
                        gene=str(row["gene"]),
                        weight=float(row["weight"]),
                    )
                )
    if not rows:
        raise ValueError(f"No top genes found in {path} with top_n={top_n}")
    return rows


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def clean_string_series(series):
    return series.astype("string").fillna("NA").astype(str)


def direction_label(value: float, threshold: float) -> str:
    if math.isnan(value):
        return "missing"
    if value >= threshold:
        return "up"
    if value <= -threshold:
        return "down"
    return "near_zero"


def consensus_label(frac_up: float, frac_down: float) -> str:
    if frac_up >= 0.60 and frac_down <= 0.20:
        return "mostly_up"
    if frac_down >= 0.60 and frac_up <= 0.20:
        return "mostly_down"
    if frac_up >= 0.30 and frac_down >= 0.30:
        return "mixed"
    return "mostly_near_zero"


def mean(values: Iterable[float]) -> float:
    vals = [float(v) for v in values if not math.isnan(float(v))]
    return sum(vals) / len(vals) if vals else float("nan")


def median(values: Iterable[float]) -> float:
    vals = sorted(float(v) for v in values if not math.isnan(float(v)))
    if not vals:
        return float("nan")
    mid = len(vals) // 2
    if len(vals) % 2:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2.0


def safe_float(value: float) -> str:
    if value is None or math.isnan(float(value)):
        return ""
    return f"{float(value):.6g}"


def numeric_sort_key(label: str, day_lookup: dict[str, float]) -> tuple[float, str]:
    return (day_lookup.get(label, float("inf")), label)


def load_expression(path: Path):
    try:
        import anndata as ad
        import numpy as np
        from scipy import sparse
    except ImportError as exc:
        raise SystemExit(
            "This script requires anndata, numpy, and scipy. Run it in the Rhino oshane-jlab "
            "environment or install those packages in a local analysis environment."
        ) from exc
    return ad, np, sparse, ad.read_h5ad(path)


def get_matrix(adata, layer: str):
    if layer == "X":
        return adata.X
    if layer not in adata.layers:
        raise ValueError(f"Requested layer {layer!r} not found. Available layers: {list(adata.layers.keys())}")
    return adata.layers[layer]


def make_group_summaries(
    *,
    adata,
    matrix,
    obs,
    gene_indices: list[int],
    gene_specs_by_gene: dict[str, list[GeneSpec]],
    np,
    sparse,
    cpm_pseudocount: float,
) -> tuple[list[dict], dict[tuple[str, str, str, str, str], dict]]:
    """Return per-group expression rows and a lookup keyed by scope/stratum/subject/timepoint/gene."""
    genes_by_index = {idx: gene for gene, idx in zip(gene_specs_by_gene, gene_indices)}
    # Preserve the requested gene order and allow duplicated genes across patterns in downstream joins.
    selected_genes = list(gene_specs_by_gene.keys())
    group_rows: list[dict] = []
    lookup: dict[tuple[str, str, str, str, str], dict] = {}

    grouping_specs = [
        ("all_cells", ["subject", "timepoint_merged"]),
        ("cell_type", ["subject", "broad_cell_type", "timepoint_merged"]),
    ]

    for scope, group_cols in grouping_specs:
        grouped = obs.groupby(group_cols, observed=True, sort=False).indices
        for group_key, row_indices in grouped.items():
            if not isinstance(group_key, tuple):
                group_key = (group_key,)
            group_values = dict(zip(group_cols, [str(v) for v in group_key]))
            subject = group_values["subject"]
            timepoint = group_values["timepoint_merged"]
            stratum = group_values.get("broad_cell_type", "all_cells")

            idx = np.fromiter(row_indices, dtype=np.int64)
            submatrix = matrix[idx, :]
            n_cells = int(len(idx))
            total_counts = float(submatrix.sum())
            gene_counts = submatrix[:, gene_indices].sum(axis=0)
            if sparse.issparse(gene_counts):
                gene_counts = gene_counts.A1
            else:
                gene_counts = np.asarray(gene_counts).ravel()
            denominator = total_counts if total_counts > 0 else 1.0
            cpm = (gene_counts / denominator) * 1_000_000.0
            log2_cpm = np.log2(cpm + cpm_pseudocount)

            for pos, gene in enumerate(selected_genes):
                row = {
                    "scope": scope,
                    "stratum": stratum,
                    "subject": subject,
                    "timepoint_merged": timepoint,
                    "day_merged_numeric": float(obs.iloc[int(idx[0])]["day_merged_numeric"]),
                    "gene": gene,
                    "n_cells": n_cells,
                    "total_counts": total_counts,
                    "gene_counts": float(gene_counts[pos]),
                    "cpm": float(cpm[pos]),
                    "log2_cpm": float(log2_cpm[pos]),
                }
                group_rows.append(row)
                lookup[(scope, stratum, subject, timepoint, gene)] = row

    return group_rows, lookup


def make_contrast_rows(
    *,
    gene_specs: list[GeneSpec],
    group_lookup: dict[tuple[str, str, str, str, str], dict],
    baseline_timepoint: str,
    min_cells: int,
    threshold: float,
) -> list[dict]:
    rows: list[dict] = []
    gene_spec_lookup = defaultdict(list)
    for spec in gene_specs:
        gene_spec_lookup[spec.gene].append(spec)

    for key, contrast in group_lookup.items():
        scope, stratum, subject, timepoint, gene = key
        baseline = group_lookup.get((scope, stratum, subject, baseline_timepoint, gene))
        if baseline is None:
            continue
        if int(baseline["n_cells"]) < min_cells or int(contrast["n_cells"]) < min_cells:
            continue
        log2fc = float(contrast["log2_cpm"]) - float(baseline["log2_cpm"])
        for spec in gene_spec_lookup[gene]:
            rows.append(
                {
                    "pattern": spec.pattern,
                    "rank": spec.rank,
                    "gene": spec.gene,
                    "weight": spec.weight,
                    "scope": scope,
                    "stratum": stratum,
                    "subject": subject,
                    "timepoint_merged": timepoint,
                    "day_merged_numeric": contrast["day_merged_numeric"],
                    "baseline_timepoint": baseline_timepoint,
                    "baseline_n_cells": baseline["n_cells"],
                    "contrast_n_cells": contrast["n_cells"],
                    "baseline_log2_cpm": baseline["log2_cpm"],
                    "contrast_log2_cpm": contrast["log2_cpm"],
                    "log2fc_vs_baseline": log2fc,
                    "direction": "baseline" if timepoint == baseline_timepoint else direction_label(log2fc, threshold),
                }
            )
    return rows


def summarize_gene_directionality(rows: list[dict], threshold: float) -> list[dict]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        if row["timepoint_merged"] == row["baseline_timepoint"]:
            continue
        key = (
            row["pattern"],
            row["rank"],
            row["gene"],
            row["weight"],
            row["scope"],
            row["stratum"],
            row["timepoint_merged"],
            row["day_merged_numeric"],
        )
        grouped[key].append(row)

    summary: list[dict] = []
    for key, group in grouped.items():
        pattern, rank, gene, weight, scope, stratum, timepoint, day = key
        values = [float(row["log2fc_vs_baseline"]) for row in group]
        n_up = sum(1 for v in values if v >= threshold)
        n_down = sum(1 for v in values if v <= -threshold)
        n_near = len(values) - n_up - n_down
        summary.append(
            {
                "pattern": pattern,
                "rank": int(rank),
                "gene": gene,
                "weight": float(weight),
                "scope": scope,
                "stratum": stratum,
                "timepoint_merged": timepoint,
                "day_merged_numeric": float(day),
                "n_subject_contrasts": len(values),
                "mean_log2fc": mean(values),
                "median_log2fc": median(values),
                "n_up": n_up,
                "n_down": n_down,
                "n_near_zero": n_near,
                "frac_up": n_up / len(values) if values else float("nan"),
                "frac_down": n_down / len(values) if values else float("nan"),
                "gene_direction_consensus": consensus_label(
                    n_up / len(values) if values else 0.0,
                    n_down / len(values) if values else 0.0,
                ),
            }
        )
    return sorted(summary, key=lambda r: (r["scope"], r["stratum"], r["pattern"], float(r["day_merged_numeric"]), int(r["rank"])))


def summarize_pattern_directionality(gene_summary_rows: list[dict], threshold: float) -> list[dict]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in gene_summary_rows:
        key = (
            row["pattern"],
            row["scope"],
            row["stratum"],
            row["timepoint_merged"],
            row["day_merged_numeric"],
        )
        grouped[key].append(row)

    summary: list[dict] = []
    for key, group in grouped.items():
        pattern, scope, stratum, timepoint, day = key
        values = [float(row["mean_log2fc"]) for row in group]
        n_up = sum(1 for v in values if v >= threshold)
        n_down = sum(1 for v in values if v <= -threshold)
        n_near = len(values) - n_up - n_down
        frac_up = n_up / len(values) if values else float("nan")
        frac_down = n_down / len(values) if values else float("nan")
        summary.append(
            {
                "pattern": pattern,
                "scope": scope,
                "stratum": stratum,
                "timepoint_merged": timepoint,
                "day_merged_numeric": float(day),
                "n_genes_tested": len(values),
                "mean_gene_log2fc": mean(values),
                "median_gene_log2fc": median(values),
                "n_genes_up": n_up,
                "n_genes_down": n_down,
                "n_genes_near_zero": n_near,
                "frac_genes_up": frac_up,
                "frac_genes_down": frac_down,
                "pattern_direction_consensus": consensus_label(frac_up if not math.isnan(frac_up) else 0.0, frac_down if not math.isnan(frac_down) else 0.0),
            }
        )
    return sorted(summary, key=lambda r: (r["scope"], r["stratum"], r["pattern"], float(r["day_merged_numeric"])))


def attach_pattern_metadata(pattern_summary_path: Path, pattern_rows: list[dict]) -> list[dict]:
    if not pattern_summary_path.exists():
        return pattern_rows
    metadata = {row["pattern"]: row for row in read_csv_dicts(pattern_summary_path)}
    for row in pattern_rows:
        meta = metadata.get(row["pattern"], {})
        for column in ["pattern_class", "peak_timepoint", "eta_timepoint", "eta_broad_cell_type", "spearman_ifn_score_rho"]:
            row[column] = meta.get(column, "")
    return pattern_rows


def color_for_value(value: float, cap: float) -> str:
    if value is None or math.isnan(float(value)):
        return "#e5e7eb"
    v = max(-cap, min(cap, float(value))) / cap
    if v >= 0:
        # white to red
        r, g, b = 185, 28, 28
        t = v
    else:
        # white to blue
        r, g, b = 37, 99, 235
        t = -v
    base = 255
    rr = round(base * (1 - t) + r * t)
    gg = round(base * (1 - t) + g * t)
    bb = round(base * (1 - t) + b * t)
    return f"#{rr:02x}{gg:02x}{bb:02x}"


def write_svg_heatmap(
    *,
    path: Path,
    title: str,
    row_labels: list[str],
    col_labels: list[str],
    values: dict[tuple[str, str], float],
    cap: float,
    row_height: int = 18,
    col_width: int = 70,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    left = 230
    top = 70
    width = left + col_width * len(col_labels) + 40
    height = top + row_height * len(row_labels) + 80
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Helvetica,Arial,sans-serif;font-size:11px;fill:#111827}.title{font-size:16px;font-weight:700}.small{font-size:10px;fill:#4b5563}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text class="title" x="20" y="28">{html.escape(title)}</text>',
        f'<text class="small" x="20" y="48">Color scale: blue=down, white=near 0, red=up; capped at +/-{cap:g} log2FC</text>',
    ]
    for j, col in enumerate(col_labels):
        x = left + j * col_width + col_width / 2
        parts.append(f'<text text-anchor="middle" x="{x:.1f}" y="{top - 12}">{html.escape(col)}</text>')
    for i, row in enumerate(row_labels):
        y = top + i * row_height
        parts.append(f'<text text-anchor="end" x="{left - 8}" y="{y + row_height - 5}">{html.escape(row)}</text>')
        for j, col in enumerate(col_labels):
            x = left + j * col_width
            value = values.get((row, col), float("nan"))
            parts.append(
                f'<rect x="{x}" y="{y}" width="{col_width}" height="{row_height}" '
                f'fill="{color_for_value(value, cap)}" stroke="#ffffff" stroke-width="1"/>'
            )
            if not math.isnan(float(value)):
                parts.append(
                    f'<text class="small" text-anchor="middle" x="{x + col_width / 2:.1f}" '
                    f'y="{y + row_height - 5}">{value:.2f}</text>'
                )
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def create_figures(
    *,
    outdir: Path,
    gene_summary_rows: list[dict],
    pattern_summary_rows: list[dict],
    gene_specs: list[GeneSpec],
    cap: float,
) -> None:
    figures_dir = outdir / "figures"
    all_cell_gene_rows = [
        row for row in gene_summary_rows if row["scope"] == "all_cells" and row["stratum"] == "all_cells"
    ]
    day_lookup = {row["timepoint_merged"]: float(row["day_merged_numeric"]) for row in all_cell_gene_rows}
    timepoints = sorted({row["timepoint_merged"] for row in all_cell_gene_rows}, key=lambda x: numeric_sort_key(x, day_lookup))

    for pattern in sorted({spec.pattern for spec in gene_specs}):
        specs = sorted([spec for spec in gene_specs if spec.pattern == pattern], key=lambda spec: spec.rank)
        rows = [f"{pattern} {spec.rank}. {spec.gene}" for spec in specs]
        values = {}
        gene_value_lookup = {
            (row["gene"], row["timepoint_merged"]): float(row["mean_log2fc"])
            for row in all_cell_gene_rows
            if row["pattern"] == pattern
        }
        for spec in specs:
            row_label = f"{pattern} {spec.rank}. {spec.gene}"
            for tp in timepoints:
                values[(row_label, tp)] = gene_value_lookup.get((spec.gene, tp), float("nan"))
        write_svg_heatmap(
            path=figures_dir / f"{pattern}_top_gene_log2fc_all_cells.svg",
            title=f"{pattern} top gene directionality, all experimental cells",
            row_labels=rows,
            col_labels=timepoints,
            values=values,
            cap=cap,
        )

    combined_specs = sorted(gene_specs, key=lambda spec: (spec.pattern, spec.rank, spec.gene))
    combined_rows = [f"{spec.pattern} {spec.rank:02d}. {spec.gene}" for spec in combined_specs]
    combined_values = {}
    combined_lookup = {
        (row["pattern"], row["gene"], row["timepoint_merged"]): float(row["mean_log2fc"])
        for row in all_cell_gene_rows
    }
    for spec in combined_specs:
        row_label = f"{spec.pattern} {spec.rank:02d}. {spec.gene}"
        for tp in timepoints:
            combined_values[(row_label, tp)] = combined_lookup.get((spec.pattern, spec.gene, tp), float("nan"))
    write_svg_heatmap(
        path=figures_dir / "all_patterns_top_gene_log2fc_all_cells.svg",
        title="All selected CoGAPS pattern genes, log2FC vs D0",
        row_labels=combined_rows,
        col_labels=timepoints,
        values=combined_values,
        cap=cap,
        row_height=16,
        col_width=72,
    )

    pattern_rows = [
        row for row in pattern_summary_rows if row["scope"] == "all_cells" and row["stratum"] == "all_cells"
    ]
    pattern_day_lookup = {row["timepoint_merged"]: float(row["day_merged_numeric"]) for row in pattern_rows}
    pattern_timepoints = sorted(
        {row["timepoint_merged"] for row in pattern_rows}, key=lambda x: numeric_sort_key(x, pattern_day_lookup)
    )
    patterns = sorted({row["pattern"] for row in pattern_rows})
    values = {
        (row["pattern"], row["timepoint_merged"]): float(row["median_gene_log2fc"])
        for row in pattern_rows
    }
    write_svg_heatmap(
        path=figures_dir / "pattern_level_median_log2fc_all_cells.svg",
        title="Pattern-level median top-gene log2FC, all experimental cells",
        row_labels=patterns,
        col_labels=pattern_timepoints,
        values=values,
        cap=cap,
    )

    cell_type_rows = [row for row in pattern_summary_rows if row["scope"] == "cell_type"]
    for pattern in sorted({row["pattern"] for row in cell_type_rows}):
        rows = [row for row in cell_type_rows if row["pattern"] == pattern]
        row_labels = sorted({row["stratum"] for row in rows})
        day_lookup = {row["timepoint_merged"]: float(row["day_merged_numeric"]) for row in rows}
        col_labels = sorted({row["timepoint_merged"] for row in rows}, key=lambda x: numeric_sort_key(x, day_lookup))
        values = {
            (row["stratum"], row["timepoint_merged"]): float(row["median_gene_log2fc"])
            for row in rows
        }
        write_svg_heatmap(
            path=figures_dir / f"{pattern}_pattern_level_log2fc_by_cell_type.svg",
            title=f"{pattern} median top-gene log2FC by cell type",
            row_labels=row_labels,
            col_labels=col_labels,
            values=values,
            cap=cap,
        )


def write_report(
    *,
    path: Path,
    args,
    expression_path: Path,
    top_genes_path: Path,
    pattern_summary_path: Path,
    n_cells: int,
    n_cohort_cells: int,
    n_genes: int,
    gene_specs: list[GeneSpec],
    pattern_summary_rows: list[dict],
) -> None:
    all_cells = [
        row for row in pattern_summary_rows if row["scope"] == "all_cells" and row["stratum"] == "all_cells"
    ]
    strongest = sorted(all_cells, key=lambda row: abs(float(row["median_gene_log2fc"])), reverse=True)[:10]
    lines = [
        "# GSE154386 Selected CoGAPS Pattern Directionality",
        "",
        "This post-sweep analysis asks whether the genes defining the selected CoGAPS patterns are up- or down-regulated relative to D0.",
        "",
        "## Inputs",
        "",
        f"- Expression object: `{expression_path}`",
        f"- Selected top genes: `{top_genes_path}`",
        f"- Pattern summary: `{pattern_summary_path}`",
        f"- Run stem: `{args.run_stem}`",
        f"- Top genes per pattern tested: `{args.top_n}`",
        f"- Baseline timepoint: `{args.baseline_timepoint}`",
        f"- Direction threshold: `abs(log2FC) >= {args.log2fc_threshold}`",
        "",
        "## Expression Matrix",
        "",
        f"- Cells analyzed before cohort filtering: `{n_cells}`",
        f"- Cells analyzed after `{args.cohort}` cohort filtering: `{n_cohort_cells}`",
        f"- HVGs available: `{n_genes}`",
        "",
        "## Interpretation Notes",
        "",
        "- CoGAPS gene weights are nonnegative, so they identify pattern-defining genes but do not encode up/down regulation by themselves.",
        "- Directionality here is estimated from pseudobulk expression contrasts using raw counts from `layers['counts']`.",
        "- The primary contrast is each post-baseline experimental timepoint versus the same subject's D0 baseline.",
        "- Cell-type-stratified outputs should be used to avoid confusing cell-type composition with within-cell-type regulation.",
        "",
        "## Strongest All-Cell Pattern-Timepoint Summaries",
        "",
        "| pattern | timepoint | median gene log2FC | frac genes up | frac genes down | consensus |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in strongest:
        lines.append(
            "| {pattern} | {timepoint_merged} | {median_gene_log2fc:.3f} | {frac_genes_up:.2f} | {frac_genes_down:.2f} | {pattern_direction_consensus} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Output Files",
            "",
            "- `gene_directionality_by_subject.csv`: subject-level log2FC contrasts for each selected gene.",
            "- `gene_directionality_summary.csv`: gene-level summaries across subjects.",
            "- `pattern_directionality_summary.csv`: pattern-level summaries of top-gene direction.",
            "- `figures/all_patterns_top_gene_log2fc_all_cells.svg`: combined heatmap of every selected top-gene row across all patterns.",
            "- `figures/`: additional per-pattern and cell-type-stratified SVG heatmaps.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = build_parser().parse_args()
    if args.top_n <= 0:
        raise ValueError("--top-n must be positive")
    if args.min_cells <= 0:
        raise ValueError("--min-cells must be positive")
    if args.cpm_pseudocount <= 0:
        raise ValueError("--cpm-pseudocount must be positive")

    sweep_outdir = args.sweep_outdir.resolve()
    expression_h5ad = (
        args.expression_h5ad.resolve()
        if args.expression_h5ad
        else sweep_outdir / "cache" / "gse154386_preprocessed_hvg.h5ad"
    )
    top_genes_csv = (
        args.top_genes_csv.resolve()
        if args.top_genes_csv
        else sweep_outdir / "runs" / f"{args.run_stem}.top_genes.csv"
    )
    pattern_summary_csv = (
        args.pattern_summary_csv.resolve()
        if args.pattern_summary_csv
        else sweep_outdir / "runs" / f"{args.run_stem}.pattern_summary.csv"
    )
    outdir = args.outdir.resolve()

    for path in [expression_h5ad, top_genes_csv, pattern_summary_csv]:
        if not path.exists():
            raise FileNotFoundError(path)

    gene_specs = read_top_genes(top_genes_csv, args.top_n)
    unique_genes = []
    for spec in gene_specs:
        if spec.gene not in unique_genes:
            unique_genes.append(spec.gene)

    ad, np, sparse, adata = load_expression(expression_h5ad)
    n_cells, n_genes = adata.shape

    required_obs = ["cohort", "subject", "timepoint_merged", "day_merged_numeric", "broad_cell_type"]
    require_columns("expression obs", adata.obs.columns, required_obs)

    obs = adata.obs.copy()
    for column in ["cohort", "subject", "timepoint_merged", "broad_cell_type"]:
        obs[column] = clean_string_series(obs[column])

    cohort_mask = obs["cohort"] == args.cohort
    if not bool(cohort_mask.any()):
        raise ValueError(f"No cells found for cohort={args.cohort!r}")
    obs = obs.loc[cohort_mask].copy()

    matrix = get_matrix(adata, args.layer)
    matrix = matrix[cohort_mask.to_numpy(), :]

    var_names = [str(v) for v in adata.var_names]
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    missing_genes = [gene for gene in unique_genes if gene not in gene_to_idx]
    present_genes = [gene for gene in unique_genes if gene in gene_to_idx]
    if not present_genes:
        raise ValueError("None of the selected top genes were present in the expression object.")

    filtered_specs = [spec for spec in gene_specs if spec.gene in gene_to_idx]
    gene_specs_by_gene: dict[str, list[GeneSpec]] = defaultdict(list)
    for spec in filtered_specs:
        gene_specs_by_gene[spec.gene].append(spec)
    gene_indices = [gene_to_idx[gene] for gene in gene_specs_by_gene.keys()]

    group_rows, group_lookup = make_group_summaries(
        adata=adata,
        matrix=matrix,
        obs=obs,
        gene_indices=gene_indices,
        gene_specs_by_gene=gene_specs_by_gene,
        np=np,
        sparse=sparse,
        cpm_pseudocount=args.cpm_pseudocount,
    )
    contrast_rows = make_contrast_rows(
        gene_specs=filtered_specs,
        group_lookup=group_lookup,
        baseline_timepoint=args.baseline_timepoint,
        min_cells=args.min_cells,
        threshold=args.log2fc_threshold,
    )
    gene_summary_rows = summarize_gene_directionality(contrast_rows, args.log2fc_threshold)
    pattern_summary_rows = summarize_pattern_directionality(gene_summary_rows, args.log2fc_threshold)
    pattern_summary_rows = attach_pattern_metadata(pattern_summary_csv, pattern_summary_rows)

    outdir.mkdir(parents=True, exist_ok=True)
    tables_dir = outdir / "tables"

    write_csv(
        tables_dir / "group_expression_log2cpm.csv",
        group_rows,
        [
            "scope",
            "stratum",
            "subject",
            "timepoint_merged",
            "day_merged_numeric",
            "gene",
            "n_cells",
            "total_counts",
            "gene_counts",
            "cpm",
            "log2_cpm",
        ],
    )
    write_csv(
        tables_dir / "gene_directionality_by_subject.csv",
        contrast_rows,
        [
            "pattern",
            "rank",
            "gene",
            "weight",
            "scope",
            "stratum",
            "subject",
            "timepoint_merged",
            "day_merged_numeric",
            "baseline_timepoint",
            "baseline_n_cells",
            "contrast_n_cells",
            "baseline_log2_cpm",
            "contrast_log2_cpm",
            "log2fc_vs_baseline",
            "direction",
        ],
    )
    write_csv(
        tables_dir / "gene_directionality_summary.csv",
        gene_summary_rows,
        [
            "pattern",
            "rank",
            "gene",
            "weight",
            "scope",
            "stratum",
            "timepoint_merged",
            "day_merged_numeric",
            "n_subject_contrasts",
            "mean_log2fc",
            "median_log2fc",
            "n_up",
            "n_down",
            "n_near_zero",
            "frac_up",
            "frac_down",
            "gene_direction_consensus",
        ],
    )
    write_csv(
        tables_dir / "pattern_directionality_summary.csv",
        pattern_summary_rows,
        [
            "pattern",
            "scope",
            "stratum",
            "timepoint_merged",
            "day_merged_numeric",
            "n_genes_tested",
            "mean_gene_log2fc",
            "median_gene_log2fc",
            "n_genes_up",
            "n_genes_down",
            "n_genes_near_zero",
            "frac_genes_up",
            "frac_genes_down",
            "pattern_direction_consensus",
            "pattern_class",
            "peak_timepoint",
            "eta_timepoint",
            "eta_broad_cell_type",
            "spearman_ifn_score_rho",
        ],
    )

    create_figures(
        outdir=outdir,
        gene_summary_rows=gene_summary_rows,
        pattern_summary_rows=pattern_summary_rows,
        gene_specs=filtered_specs,
        cap=args.heatmap_abs_cap,
    )

    manifest = {
        "sweep_outdir": str(sweep_outdir),
        "run_stem": args.run_stem,
        "expression_h5ad": str(expression_h5ad),
        "top_genes_csv": str(top_genes_csv),
        "pattern_summary_csv": str(pattern_summary_csv),
        "outdir": str(outdir),
        "top_n": args.top_n,
        "baseline_timepoint": args.baseline_timepoint,
        "cohort": args.cohort,
        "layer": args.layer,
        "min_cells": args.min_cells,
        "log2fc_threshold": args.log2fc_threshold,
        "cpm_pseudocount": args.cpm_pseudocount,
        "expression_shape": [int(n_cells), int(n_genes)],
        "n_selected_gene_pattern_rows": len(gene_specs),
        "n_unique_selected_genes": len(unique_genes),
        "n_present_selected_genes": len(present_genes),
        "missing_selected_genes": missing_genes,
        "output_tables": {
            "group_expression_log2cpm": str(tables_dir / "group_expression_log2cpm.csv"),
            "gene_directionality_by_subject": str(tables_dir / "gene_directionality_by_subject.csv"),
            "gene_directionality_summary": str(tables_dir / "gene_directionality_summary.csv"),
            "pattern_directionality_summary": str(tables_dir / "pattern_directionality_summary.csv"),
        },
    }
    (outdir / "directionality_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    write_report(
        path=outdir / "README_directionality.md",
        args=args,
        expression_path=expression_h5ad,
        top_genes_path=top_genes_csv,
        pattern_summary_path=pattern_summary_csv,
        n_cells=int(n_cells),
        n_cohort_cells=int(obs.shape[0]),
        n_genes=int(n_genes),
        gene_specs=filtered_specs,
        pattern_summary_rows=pattern_summary_rows,
    )

    print(f"[done] wrote directionality results to {outdir}")
    print(f"[done] selected top-gene rows: {len(gene_specs)}; present: {len(filtered_specs)}; missing unique genes: {len(missing_genes)}")


if __name__ == "__main__":
    main()
