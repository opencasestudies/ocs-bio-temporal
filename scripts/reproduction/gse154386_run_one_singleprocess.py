#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import os
import sys
import time
import traceback
from pathlib import Path

from gse154386_cogaps_sweep_common import (
    DEFAULT_OUTDIR,
    REPO_ROOT,
    cache_paths,
    find_case_study5_script,
    load_reference_module,
    metrics_status_ok,
    resolve_path,
    run_artifact_paths,
    set_thread_env,
    utc_now_iso,
    write_json,
)


class TeeStream:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text: str) -> int:
        for stream in self.streams:
            stream.write(text)
            stream.flush()
        return len(text)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


@contextlib.contextmanager
def tee_output(log_path: Path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as handle:
        stdout_orig = sys.stdout
        stderr_orig = sys.stderr
        sys.stdout = TeeStream(stdout_orig, handle)
        sys.stderr = TeeStream(stderr_orig, handle)
        try:
            yield
        finally:
            sys.stdout = stdout_orig
            sys.stderr = stderr_orig


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run exactly one non-distributed GSE154386 CoGAPS configuration from the cached "
            "Case Study 5 prep outputs."
        )
    )
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Sweep output directory.")
    parser.add_argument(
        "--case-study-script",
        default=None,
        help="Optional override for the local Case Study 5 source script path.",
    )
    parser.add_argument(
        "--cogaps-input-h5ad",
        default=None,
        help="Prep-generated genes x cells AnnData for CoGAPS. Defaults to the sweep cache file.",
    )
    parser.add_argument(
        "--discovery-h5ad",
        default=None,
        help="Prep-generated experimental discovery AnnData in cells x genes orientation.",
    )
    parser.add_argument("--k", type=int, required=True, help="CoGAPS nPatterns value.")
    parser.add_argument("--seed", type=int, required=True, help="Random seed for CoGAPS.")
    parser.add_argument("--n-iter", type=int, required=True, help="CoGAPS nIterations value.")
    parser.add_argument(
        "--cogaps-threads",
        type=int,
        default=None,
        help="Thread count passed to CoGAPS nThreads. Defaults to available CPU count.",
    )
    parser.add_argument(
        "--blas-threads",
        type=int,
        default=None,
        help="Thread count exported to BLAS/OpenMP environment variables.",
    )
    parser.add_argument(
        "--top-genes",
        type=int,
        default=50,
        help="Number of top genes to record per pattern in CSV and metrics outputs.",
    )
    parser.add_argument(
        "--use-sparse-optimization",
        dest="use_sparse_optimization",
        action="store_true",
        default=True,
        help="Pass useSparseOptimization=True to CoGAPS. Enabled by default.",
    )
    parser.add_argument(
        "--no-sparse-optimization",
        dest="use_sparse_optimization",
        action="store_false",
        help="Disable CoGAPS sparse optimization.",
    )
    parser.add_argument(
        "--force-rerun",
        action="store_true",
        help="Rerun even if a previous metrics JSON already reports status=ok.",
    )
    return parser


def pattern_effectiveness_summary(gene_patterns_df, cell_patterns_df, np):
    rows = []
    for pattern in gene_patterns_df.columns:
        gene_vec = gene_patterns_df[pattern].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)
        cell_vec = cell_patterns_df[pattern].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)
        gene_max = float(np.max(gene_vec)) if gene_vec.size else 0.0
        cell_max = float(np.max(cell_vec)) if cell_vec.size else 0.0
        nonzero_gene_count = int(np.count_nonzero(gene_vec > 0))
        effective = bool(gene_max > 0.0 and cell_max > 0.0 and float(np.std(cell_vec)) > 1e-8 and nonzero_gene_count >= 10)
        rows.append(
            {
                "pattern": str(pattern),
                "effective": effective,
                "gene_loading_max": gene_max,
                "cell_score_max": cell_max,
                "nonzero_gene_count": nonzero_gene_count,
            }
        )
    effective_count = sum(1 for row in rows if row["effective"])
    degenerate = [row["pattern"] for row in rows if not row["effective"]]
    return rows, effective_count, degenerate


def within_run_redundancy(gene_patterns_df, top_gene_df, np):
    corr = gene_patterns_df.corr(method="spearman").to_numpy(dtype=float)
    if corr.shape[0] <= 1:
        return {
            "within_run_pattern_redundancy_mean": 0.0,
            "within_run_pattern_redundancy_max": 0.0,
            "within_run_top_gene_jaccard_mean": 0.0,
            "within_run_top_gene_jaccard_max": 0.0,
        }

    tri = corr[np.triu_indices_from(corr, k=1)]
    tri = np.nan_to_num(tri, nan=0.0, posinf=0.0, neginf=0.0)
    tri = np.clip(tri, 0.0, None)

    top_gene_sets = {
        str(pattern): set(group["gene"].astype(str).tolist())
        for pattern, group in top_gene_df.groupby("pattern", observed=True)
    }
    pattern_names = list(gene_patterns_df.columns)
    jaccards = []
    for i, left in enumerate(pattern_names):
        for right in pattern_names[i + 1:]:
            left_set = top_gene_sets.get(str(left), set())
            right_set = top_gene_sets.get(str(right), set())
            union = left_set.union(right_set)
            jaccards.append(float(len(left_set.intersection(right_set)) / len(union)) if union else 0.0)

    return {
        "within_run_pattern_redundancy_mean": float(np.mean(tri)) if tri.size else 0.0,
        "within_run_pattern_redundancy_max": float(np.max(tri)) if tri.size else 0.0,
        "within_run_top_gene_jaccard_mean": float(np.mean(jaccards)) if jaccards else 0.0,
        "within_run_top_gene_jaccard_max": float(np.max(jaccards)) if jaccards else 0.0,
    }


def main() -> None:
    args = build_parser().parse_args()
    if args.k <= 0:
        raise ValueError("--k must be positive")
    if args.seed < 0:
        raise ValueError("--seed must be non-negative")
    if args.n_iter <= 0:
        raise ValueError("--n-iter must be positive")
    if args.top_genes <= 0:
        raise ValueError("--top-genes must be positive")

    outdir = resolve_path(args.outdir)
    cache = cache_paths(outdir)
    source_script = resolve_path(args.case_study_script) if args.case_study_script else find_case_study5_script()
    cogaps_input_h5ad = resolve_path(args.cogaps_input_h5ad) if args.cogaps_input_h5ad else cache["cogaps_input_h5ad"]
    discovery_h5ad = resolve_path(args.discovery_h5ad) if args.discovery_h5ad else cache["experimental_discovery_h5ad"]

    artifacts = run_artifact_paths(outdir=outdir, k=args.k, seed=args.seed, n_iter=args.n_iter)
    if metrics_status_ok(artifacts["metrics_path"]) and not args.force_rerun:
        print(f"[skip] existing successful metrics found: {artifacts['metrics_path']}")
        return

    with tee_output(artifacts["log_path"]):
        print(f"[run] started at {utc_now_iso()}")
        print(f"[run] outdir={outdir}")
        print(f"[run] source_script={source_script}")
        print(f"[run] cogaps_input_h5ad={cogaps_input_h5ad}")
        print(f"[run] discovery_h5ad={discovery_h5ad}")
        print(f"[run] parameters: K={args.k} seed={args.seed} n_iter={args.n_iter}")

        run_start = time.time()
        try:
            os.chdir(REPO_ROOT)
            requested_threads = args.cogaps_threads or max(1, os.cpu_count() or 1)
            set_thread_env(args.blas_threads or requested_threads)

            import numpy as np  # lazy import so --help works without cluster deps
            import scanpy as sc

            module = load_reference_module(source_script)
            adata_discovery_cogaps = sc.read_h5ad(cogaps_input_h5ad)
            adata_discovery = sc.read_h5ad(discovery_h5ad)

            module.COGAPS_N_PATTERNS = int(args.k)
            module.COGAPS_N_ITER = int(args.n_iter)
            module.RANDOM_SEED = int(args.seed)
            module.COGAPS_N_THREADS = int(requested_threads)
            module.COGAPS_USE_SPARSE_OPT = bool(args.use_sparse_optimization)

            result = module.run_cogaps_nondistributed(adata_discovery_cogaps, artifacts["result_path"])
            gene_patterns_df, cell_patterns_df, orientation_msg = module.extract_cogaps_matrices(result, adata_discovery)
            print(f"[run] orientation: {orientation_msg}")

            gene_patterns_df = module.clean_pattern_df(gene_patterns_df, clip_negative=True)
            cell_patterns_df = module.clean_pattern_df(cell_patterns_df, clip_negative=True)
            top_gene_df = module.top_genes_per_pattern(gene_patterns_df, top_n=args.top_genes)
            cell_meta_patterns = adata_discovery.obs.join(cell_patterns_df, how="inner")
            pattern_summary_df = module.summarize_cogaps_patterns(cell_meta_patterns, top_gene_df)

            gene_patterns_df.to_csv(artifacts["gene_loadings_csv"])
            cell_patterns_df.to_csv(artifacts["cell_scores_csv"])
            top_gene_df.to_csv(artifacts["top_genes_csv"], index=False)
            pattern_summary_df.to_csv(artifacts["pattern_summary_csv"], index=False)
            cell_meta_patterns.to_csv(artifacts["cell_metadata_patterns_csv"])

            effective_rows, n_patterns_effective, degenerate_patterns = pattern_effectiveness_summary(
                gene_patterns_df=gene_patterns_df,
                cell_patterns_df=cell_patterns_df,
                np=np,
            )
            redundancy = within_run_redundancy(
                gene_patterns_df=gene_patterns_df,
                top_gene_df=top_gene_df,
                np=np,
            )

            pattern_counts = pattern_summary_df["pattern_class"].astype(str).value_counts(dropna=False)
            activity_like_count = int(pattern_counts.get("activity-like", 0))
            identity_like_count = int(pattern_counts.get("identity-like", 0))
            mixed_count = int(pattern_counts.get("mixed", 0))
            unclassified_count = int(pattern_counts.get("unclassified", 0))
            pattern_mix_balance = float(min(activity_like_count, identity_like_count) / max(1, args.k))

            top_genes_by_pattern = {
                str(pattern): group["gene"].astype(str).tolist()
                for pattern, group in top_gene_df.groupby("pattern", observed=True)
            }
            eta_series = pattern_summary_df["eta_timepoint"].replace([np.inf, -np.inf], np.nan).dropna()
            ifn_overlap_series = pattern_summary_df["ifn_top_gene_overlap_top15"].replace([np.inf, -np.inf], np.nan).dropna()
            ifn_rho_series = pattern_summary_df["spearman_ifn_score_rho"].replace([np.inf, -np.inf], np.nan).dropna()

            runtime_sec = float(time.time() - run_start)
            metrics = {
                "created_at_utc": utc_now_iso(),
                "status": "ok",
                "source_script_path": str(source_script),
                "K": int(args.k),
                "seed": int(args.seed),
                "n_iter": int(args.n_iter),
                "runtime_sec": runtime_sec,
                "transposeData": False,
                "distributed_mode": False,
                "use_sparse_optimization": bool(args.use_sparse_optimization),
                "cogaps_threads": int(requested_threads),
                "blas_threads": int(args.blas_threads or requested_threads),
                "orientation_msg": orientation_msg,
                "result_path": str(artifacts["result_path"]),
                "metrics_path": str(artifacts["metrics_path"]),
                "log_path": str(artifacts["log_path"]),
                "gene_loadings_csv": str(artifacts["gene_loadings_csv"]),
                "cell_scores_csv": str(artifacts["cell_scores_csv"]),
                "top_genes_csv": str(artifacts["top_genes_csv"]),
                "pattern_summary_csv": str(artifacts["pattern_summary_csv"]),
                "cell_metadata_patterns_csv": str(artifacts["cell_metadata_patterns_csv"]),
                "n_patterns_effective": int(n_patterns_effective),
                "degenerate_patterns": degenerate_patterns,
                "top_genes_by_pattern": top_genes_by_pattern,
                "effective_pattern_details": effective_rows,
                "activity_like_pattern_count": activity_like_count,
                "identity_like_pattern_count": identity_like_count,
                "mixed_pattern_count": mixed_count,
                "unclassified_pattern_count": unclassified_count,
                "pattern_mix_balance": pattern_mix_balance,
                "candidate_ifn_pattern_count": int(pattern_summary_df["candidate_ifn_pattern"].fillna(False).astype(bool).sum()),
                "max_ifn_top_gene_overlap_top15": int(ifn_overlap_series.max()) if len(ifn_overlap_series) else 0,
                "max_ifn_score_rho": float(ifn_rho_series.max()) if len(ifn_rho_series) else 0.0,
                "mean_eta_timepoint": float(eta_series.mean()) if len(eta_series) else 0.0,
                "max_eta_timepoint": float(eta_series.max()) if len(eta_series) else 0.0,
                **redundancy,
            }
            write_json(artifacts["metrics_path"], metrics)
            print(f"[run] wrote result: {artifacts['result_path']}")
            print(f"[run] wrote metrics: {artifacts['metrics_path']}")
        except Exception as exc:
            failure_metrics = {
                "created_at_utc": utc_now_iso(),
                "status": "error",
                "source_script_path": str(source_script),
                "K": int(args.k),
                "seed": int(args.seed),
                "n_iter": int(args.n_iter),
                "runtime_sec": float(time.time() - run_start),
                "result_path": str(artifacts["result_path"]),
                "metrics_path": str(artifacts["metrics_path"]),
                "log_path": str(artifacts["log_path"]),
                "error_type": exc.__class__.__name__,
                "error_message": str(exc),
                "traceback": traceback.format_exc(),
            }
            write_json(artifacts["metrics_path"], failure_metrics)
            print(f"[run] wrote failure metrics: {artifacts['metrics_path']}", file=sys.stderr)
            raise


if __name__ == "__main__":
    main()
