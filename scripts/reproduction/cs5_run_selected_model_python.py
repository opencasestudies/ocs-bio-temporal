#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Dict

import anndata as ad
import numpy as np
import pandas as pd

from PyCoGAPS.parameters import CoParams, setParams
from PyCoGAPS.helper_functions import isCompiledWithOpenMPSupport
from PyCoGAPS.pycogaps_main import CoGAPS

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
WORKSPACE_ROOT = REPO_ROOT
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

from cs5_cogaps_python_utils import (  # noqa: E402
    atomic_write_h5ad,
    atomic_write_json,
    effective_pattern_summary,
    redundancy_summary,
    setup_logger,
    set_thread_env,
    summarize_result_uns,
    top_genes_by_pattern,
    trace_dataframe,
    utc_now_iso,
)


DEFAULT_SWEEP_DIR = Path("/home/rstudio/project/data/cs5_sweep_results")
DEFAULT_WORKSPACE_DIR = Path("/home/rstudio/project")
DEFAULT_OUTDIR = DEFAULT_WORKSPACE_DIR / "GSE154386" / "cogaps_python_revised_model_K10_seed2_iter2000"
DEFAULT_SOURCE_SCRIPT = SCRIPT_DIR / "gse154386_sparse_distributed_cogaps.py"


def load_source_module(source_script: Path):
    if not source_script.exists():
        raise FileNotFoundError(f"Case Study 5 source script does not exist: {source_script}")
    module_name = f"cs5_source_{abs(hash(str(source_script)))}"
    spec = importlib.util.spec_from_file_location(module_name, source_script)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load import spec for {source_script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the Case Study 5 selected PyCoGAPS model with full diagnostics.")
    parser.add_argument(
        "--cogaps-input-h5ad",
        default=str(DEFAULT_SWEEP_DIR / "cache" / "gse154386_experimental_discovery_genes_x_cells.h5ad"),
        help="Genes x cells AnnData for PyCoGAPS.",
    )
    parser.add_argument(
        "--discovery-h5ad",
        default=str(DEFAULT_SWEEP_DIR / "cache" / "gse154386_experimental_discovery_cells_x_genes.h5ad"),
        help="Cells x genes discovery AnnData used for metadata and orientation checks.",
    )
    parser.add_argument("--source-script", default=str(DEFAULT_SOURCE_SCRIPT), help="Case Study 5 analysis source module.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory for Python selected-model artifacts.")
    parser.add_argument("--k", type=int, default=10, help="Selected CoGAPS K.")
    parser.add_argument("--seed", type=int, default=2, help="Selected CoGAPS seed.")
    parser.add_argument("--n-iter", type=int, default=2000, help="Selected CoGAPS iteration count.")
    parser.add_argument("--top-genes", type=int, default=50, help="Top genes exported per pattern.")
    parser.add_argument("--cogaps-threads", type=int, default=4, help="Threads passed to CoGAPS.")
    parser.add_argument("--blas-threads", type=int, default=1, help="Threads exported to BLAS-like libraries.")
    parser.add_argument("--output-frequency", type=int, default=100, help="CoGAPS trace output frequency.")
    parser.add_argument("--checkpoint-interval", type=int, default=500, help="CoGAPS checkpoint interval.")
    parser.add_argument("--checkpoint-out-file", default="", help="Optional checkpoint path. Defaults to outdir/tag.checkpoint.out.")
    parser.add_argument("--checkpoint-in-file", default="", help="Optional checkpoint input path for resuming.")
    parser.add_argument("--n-snapshots", type=int, default=10, help="Snapshots per enabled phase.")
    parser.add_argument(
        "--snapshot-phase",
        default="all",
        choices=["sampling", "equilibration", "all"],
        help="CoGAPS snapshot phase.",
    )
    parser.add_argument("--take-pump-samples", action="store_true", default=True, help="Collect pump diagnostics.")
    parser.add_argument("--no-pump-samples", action="store_false", dest="take_pump_samples", help="Disable pump diagnostics.")
    parser.add_argument("--asynchronous-updates", action="store_true", default=True, help="Enable asynchronous updates.")
    parser.add_argument("--sync-updates", action="store_false", dest="asynchronous_updates", help="Disable asynchronous updates.")
    parser.add_argument("--use-sparse-opt", action="store_true", default=True, help="Enable sparse optimization.")
    parser.add_argument("--no-sparse-opt", action="store_false", dest="use_sparse_opt", help="Disable sparse optimization.")
    parser.add_argument("--force-rerun", action="store_true", help="Rerun even if an ok metrics file exists.")
    return parser


def write_case_study_tables(
    *,
    source_module,
    outdir: Path,
    stem: str,
    result: ad.AnnData,
    discovery: ad.AnnData,
    top_genes: int,
) -> Dict[str, Any]:
    gene_patterns_df, cell_patterns_df, orientation_msg = source_module.extract_cogaps_matrices(result, discovery)
    gene_patterns_df = source_module.clean_pattern_df(gene_patterns_df, clip_negative=True)
    cell_patterns_df = source_module.clean_pattern_df(cell_patterns_df, clip_negative=True)
    top_gene_df = source_module.top_genes_per_pattern(gene_patterns_df, top_n=top_genes)
    cell_meta_patterns = discovery.obs.join(cell_patterns_df, how="inner")
    pattern_summary_df = source_module.summarize_cogaps_patterns(cell_meta_patterns, top_gene_df)

    pattern_names = list(gene_patterns_df.columns)
    gene_loadings_csv = outdir / f"{stem}.gene_loadings.csv"
    cell_scores_csv = outdir / f"{stem}.cell_scores.csv"
    top_genes_csv = outdir / f"{stem}.top_genes.csv"
    cell_metadata_patterns_csv = outdir / f"{stem}.discovery_cells_with_patterns.csv"
    pattern_summary_csv = outdir / f"{stem}.pattern_summary.csv"
    rq1_csv = outdir / "RQ1_time_varying_patterns.csv"
    rq2_csv = outdir / "RQ2_identity_vs_activity_patterns.csv"
    rq3_csv = outdir / "RQ3_interferon_candidate_patterns.csv"
    means_by_time_csv = outdir / "experimental_pattern_means_by_time.csv"
    means_by_time_celltype_csv = outdir / "experimental_pattern_means_by_time_and_celltype.csv"

    gene_patterns_df.to_csv(gene_loadings_csv)
    cell_patterns_df.to_csv(cell_scores_csv)
    top_gene_df.to_csv(top_genes_csv, index=False)
    cell_meta_patterns.to_csv(cell_metadata_patterns_csv)
    pattern_summary_df.to_csv(pattern_summary_csv, index=False)
    pattern_summary_df.sort_values("eta_timepoint", ascending=False).to_csv(rq1_csv, index=False)
    pattern_summary_df[
        ["pattern", "eta_timepoint", "eta_broad_cell_type", "pattern_class", "peak_timepoint"]
    ].sort_values(["pattern_class", "eta_timepoint"], ascending=[True, False]).to_csv(rq2_csv, index=False)
    pattern_summary_df.loc[pattern_summary_df["candidate_ifn_pattern"]].to_csv(rq3_csv, index=False)

    (
        cell_meta_patterns.groupby(["timepoint_merged", "day_merged_numeric"], observed=True)[pattern_names]
        .mean()
        .reset_index()
        .sort_values("day_merged_numeric")
        .to_csv(means_by_time_csv, index=False)
    )
    (
        cell_meta_patterns.groupby(["broad_cell_type", "timepoint_merged", "day_merged_numeric"], observed=True)[pattern_names]
        .mean()
        .reset_index()
        .sort_values(["broad_cell_type", "day_merged_numeric"])
        .to_csv(means_by_time_celltype_csv, index=False)
    )

    effective_rows, n_patterns_effective, degenerate_patterns = effective_pattern_summary(gene_patterns_df, cell_patterns_df)
    redundancy = redundancy_summary(gene_patterns_df, top_gene_df)
    pattern_counts = pattern_summary_df["pattern_class"].astype(str).value_counts(dropna=False)
    eta_series = pattern_summary_df["eta_timepoint"].replace([np.inf, -np.inf], np.nan).dropna()
    ifn_overlap_series = pattern_summary_df["ifn_top_gene_overlap_top15"].replace([np.inf, -np.inf], np.nan).dropna()
    ifn_rho_series = pattern_summary_df["spearman_ifn_score_rho"].replace([np.inf, -np.inf], np.nan).dropna()

    return {
        "orientation_msg": orientation_msg,
        "pattern_names": pattern_names,
        "gene_loadings_csv": str(gene_loadings_csv),
        "cell_scores_csv": str(cell_scores_csv),
        "top_genes_csv": str(top_genes_csv),
        "cell_metadata_patterns_csv": str(cell_metadata_patterns_csv),
        "pattern_summary_csv": str(pattern_summary_csv),
        "rq1_time_varying_patterns_csv": str(rq1_csv),
        "rq2_identity_vs_activity_patterns_csv": str(rq2_csv),
        "rq3_interferon_candidate_patterns_csv": str(rq3_csv),
        "experimental_pattern_means_by_time_csv": str(means_by_time_csv),
        "experimental_pattern_means_by_time_and_celltype_csv": str(means_by_time_celltype_csv),
        "n_patterns_effective": int(n_patterns_effective),
        "degenerate_patterns": degenerate_patterns,
        "effective_pattern_details": effective_rows,
        "top_genes_by_pattern": top_genes_by_pattern(top_gene_df),
        "activity_like_pattern_count": int(pattern_counts.get("activity-like", 0)),
        "identity_like_pattern_count": int(pattern_counts.get("identity-like", 0)),
        "mixed_pattern_count": int(pattern_counts.get("mixed", 0)),
        "unclassified_pattern_count": int(pattern_counts.get("unclassified", 0)),
        "candidate_ifn_pattern_count": int(pattern_summary_df["candidate_ifn_pattern"].fillna(False).astype(bool).sum()),
        "max_ifn_top_gene_overlap_top15": int(ifn_overlap_series.max()) if len(ifn_overlap_series) else 0,
        "max_ifn_score_rho": float(ifn_rho_series.max()) if len(ifn_rho_series) else 0.0,
        "mean_eta_timepoint": float(eta_series.mean()) if len(eta_series) else 0.0,
        "max_eta_timepoint": float(eta_series.max()) if len(eta_series) else 0.0,
        **redundancy,
    }


def main() -> None:
    args = build_parser().parse_args()
    if args.k <= 0 or args.n_iter <= 0 or args.top_genes <= 0:
        raise ValueError("--k, --n-iter, and --top-genes must be positive")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    tag = f"K{args.k}_seed{args.seed}_iter{args.n_iter}"
    stem = f"cogaps_{tag}"
    log_path = outdir / f"{stem}.log"
    logger = setup_logger(log_path)

    result_path = outdir / f"{stem}.h5ad"
    metrics_path = outdir / f"{stem}.metrics.json"
    diagnostics_path = outdir / f"{stem}.diagnostics.json"
    trace_path = outdir / f"{stem}.trace.csv"
    checkpoint_path = Path(args.checkpoint_out_file) if args.checkpoint_out_file else outdir / f"{stem}.checkpoint.out"
    snapshot_frequency = 0 if args.n_snapshots <= 0 else max(1, int(args.n_iter) // int(args.n_snapshots))

    if result_path.exists() and metrics_path.exists() and not args.force_rerun:
        try:
            prior = json.loads(metrics_path.read_text(encoding="utf-8"))
            if prior.get("status") == "ok":
                logger.info("[CACHE] status=ok; skipping %s", tag)
                return
        except Exception:
            logger.info("[CACHE] Existing metrics are unreadable; rerunning.")

    set_thread_env(cogaps_threads=args.cogaps_threads, blas_threads=args.blas_threads)
    source_script = Path(args.source_script).resolve()
    cogaps_input_h5ad = Path(args.cogaps_input_h5ad).resolve()
    discovery_h5ad = Path(args.discovery_h5ad).resolve()
    source_module = load_source_module(source_script)

    payload: Dict[str, Any] = {
        "created_at_utc": utc_now_iso(),
        "status": "ok",
        "language": "Python",
        "package": "PyCoGAPS",
        "K": int(args.k),
        "seed": int(args.seed),
        "n_iter": int(args.n_iter),
        "transposeData": False,
        "distributed_mode": False,
        "use_sparse_optimization": bool(args.use_sparse_opt),
        "cogaps_threads": int(args.cogaps_threads),
        "blas_threads": int(args.blas_threads),
        "outputFrequency": int(args.output_frequency),
        "checkpointInterval": int(args.checkpoint_interval),
        "checkpointOutFile": str(checkpoint_path),
        "checkpointInFile": args.checkpoint_in_file,
        "nSnapshots": int(args.n_snapshots),
        "snapshotFrequency": int(snapshot_frequency),
        "snapshotPhase": args.snapshot_phase,
        "takePumpSamples": bool(args.take_pump_samples),
        "asynchronousUpdates": bool(args.asynchronous_updates),
        "source_script_path": str(source_script),
        "cogaps_input_h5ad": str(cogaps_input_h5ad),
        "discovery_h5ad": str(discovery_h5ad),
        "result_path": str(result_path),
        "metrics_path": str(metrics_path),
        "diagnostics_path": str(diagnostics_path),
        "trace_path": str(trace_path),
        "log_path": str(log_path),
    }

    t0 = time.time()
    try:
        pycogaps_openmp_support = bool(isCompiledWithOpenMPSupport())
        payload["pycogapsOpenMPSupport"] = pycogaps_openmp_support
        if args.cogaps_threads > 1 and not pycogaps_openmp_support:
            raise RuntimeError(
                "PyCoGAPS was built without OpenMP support, so this run cannot honor "
                f"--cogaps-threads={args.cogaps_threads}. Rebuild the runtime with "
                "OpenMP-enabled PyCoGAPS before running the selected model."
            )

        logger.info("[RUN] Starting selected PyCoGAPS model %s", tag)
        logger.info("[RUN] cogaps_input_h5ad=%s", cogaps_input_h5ad)
        logger.info("[RUN] discovery_h5ad=%s", discovery_h5ad)
        cogaps_input = ad.read_h5ad(cogaps_input_h5ad)
        discovery = ad.read_h5ad(discovery_h5ad)

        params = CoParams(adata=cogaps_input)
        setParams(
            params,
            {
                "nPatterns": int(args.k),
                "nIterations": int(args.n_iter),
                "seed": int(args.seed),
                "useSparseOptimization": bool(args.use_sparse_opt),
                "transposeData": False,
                "nThreads": int(args.cogaps_threads),
                "takePumpSamples": bool(args.take_pump_samples),
            },
        )
        params.printParams()

        result = CoGAPS(
            cogaps_input,
            params,
            nThreads=int(args.cogaps_threads),
            outputFrequency=int(args.output_frequency),
            checkpointOutFile=str(checkpoint_path),
            checkpointInterval=int(args.checkpoint_interval),
            checkpointInFile=args.checkpoint_in_file,
            asynchronousUpdates=(True if args.asynchronous_updates else None),
            nSnapshots=int(snapshot_frequency),
            snapshotPhase=args.snapshot_phase,
        )
        if not isinstance(result, ad.AnnData):
            raise TypeError(f"Expected AnnData from PyCoGAPS, got {type(result)}")

        atomic_write_h5ad(result, result_path)
        diagnostics = summarize_result_uns(result)
        trace_dataframe(result).to_csv(trace_path, index=False)
        atomic_write_json(diagnostics, diagnostics_path)
        table_metrics = write_case_study_tables(
            source_module=source_module,
            outdir=outdir,
            stem=stem,
            result=result,
            discovery=discovery,
            top_genes=int(args.top_genes),
        )
        payload.update(diagnostics)
        payload.update(table_metrics)
        payload["runtime_sec"] = float(time.time() - t0)
        payload["completed_at_utc"] = utc_now_iso()
        atomic_write_json(payload, metrics_path)
        logger.info("[DONE] %s finished in %.2f min", tag, payload["runtime_sec"] / 60)
        logger.info("[DONE] metrics=%s", metrics_path)
    except Exception as exc:
        payload["status"] = "error"
        payload["runtime_sec"] = float(time.time() - t0)
        payload["completed_at_utc"] = utc_now_iso()
        payload["error_type"] = exc.__class__.__name__
        payload["error_message"] = str(exc)
        payload["traceback"] = traceback.format_exc()
        atomic_write_json(payload, metrics_path)
        logger.error("[ERROR] selected PyCoGAPS run failed:\n%s", payload["traceback"])
        raise


if __name__ == "__main__":
    main()
