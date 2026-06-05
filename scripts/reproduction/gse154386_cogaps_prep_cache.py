#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict

from gse154386_cogaps_sweep_common import (
    DEFAULT_OUTDIR,
    REPO_ROOT,
    cache_paths,
    ensure_dir,
    find_case_study5_script,
    load_reference_module,
    read_json,
    resolve_path,
    summarize_missing_paths,
    utc_now_iso,
    write_json,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare reusable Case Study 5 CoGAPS sweep caches by reusing the local "
            "preprocessing and discovery-set logic from the existing monolithic script."
        )
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTDIR),
        help="Output directory for sweep cache files and manifest.",
    )
    parser.add_argument(
        "--case-study-script",
        default=None,
        help=(
            "Optional override for the local Case Study 5 source script. "
            "If omitted, the script auto-detects the on-disk file."
        ),
    )
    parser.add_argument(
        "--discovery-group",
        default="sample_id",
        help="obs column used for balanced discovery sampling. Default preserves the local script.",
    )
    parser.add_argument(
        "--discovery-max-cells-per-sample",
        type=int,
        default=200,
        help="Maximum cells per discovery sampling group. Use the local script default unless intentionally changing it.",
    )
    parser.add_argument(
        "--discovery-seed",
        type=int,
        default=42,
        help="Random seed for balanced discovery sampling.",
    )
    parser.add_argument(
        "--force-refresh",
        action="store_true",
        help="Rebuild the sweep cache outputs even if an existing manifest matches.",
    )
    parser.add_argument(
        "--force-rebuild-hvg",
        action="store_true",
        help="Ignore the existing local HVG cache and rebuild preprocessing from the source script.",
    )
    parser.add_argument(
        "--reuse-local-hvg-cache",
        dest="reuse_local_hvg_cache",
        action="store_true",
        default=True,
        help="Reuse GSE154386/GSE154386_all_preprocessed_HVG.h5ad when available. Enabled by default.",
    )
    parser.add_argument(
        "--no-reuse-local-hvg-cache",
        dest="reuse_local_hvg_cache",
        action="store_false",
        help="Do not reuse the existing local HVG cache; rebuild or use the sweep-local cache instead.",
    )
    return parser


def manifest_signature(args: argparse.Namespace, source_script: Path) -> Dict[str, object]:
    return {
        "source_script_path": str(source_script),
        "discovery_group": str(args.discovery_group),
        "discovery_max_cells_per_sample": int(args.discovery_max_cells_per_sample),
        "discovery_seed": int(args.discovery_seed),
        "reuse_local_hvg_cache": bool(args.reuse_local_hvg_cache),
    }


def maybe_skip_existing_cache(paths: Dict[str, Path], expected_signature: Dict[str, object], force_refresh: bool) -> bool:
    if force_refresh or not paths["prep_manifest_json"].exists():
        return False

    missing = summarize_missing_paths(
        [
            paths["preprocessed_hvg_h5ad"],
            paths["experimental_discovery_h5ad"],
            paths["cogaps_input_h5ad"],
            paths["natural_target_h5ad"],
            paths["discovery_counts_csv"],
        ]
    )
    if missing:
        print("[cache] manifest exists but cache files are missing; rebuilding:")
        for path in missing:
            print(f"  - {path}")
        return False

    manifest = read_json(paths["prep_manifest_json"])
    if manifest.get("prep_parameters") != expected_signature:
        print("[cache] existing prep manifest parameters do not match the requested configuration; rebuilding")
        return False

    print(f"[cache] existing prep cache is valid: {paths['prep_manifest_json']}")
    return True


def main() -> None:
    args = build_parser().parse_args()
    if args.discovery_max_cells_per_sample <= 0:
        raise ValueError("--discovery-max-cells-per-sample must be positive")

    outdir = resolve_path(args.outdir)
    ensure_dir(outdir)
    paths = cache_paths(outdir)

    source_script = resolve_path(args.case_study_script) if args.case_study_script else find_case_study5_script()
    expected_signature = manifest_signature(args, source_script)

    if maybe_skip_existing_cache(paths=paths, expected_signature=expected_signature, force_refresh=args.force_refresh):
        return

    os.chdir(REPO_ROOT)

    module = load_reference_module(source_script)
    import scanpy as sc  # lazy import so --help works without cluster deps

    source_hvg_cache = resolve_path(getattr(module, "OUT_ALL_HVG"))
    print(f"[source] using Case Study 5 script: {source_script}")
    if source_script.name != "gse154386_sparse_nondistributed_cogaps.py":
        print(
            "[source] note: the task docs mention gse154386_sparse_nondistributed_cogaps.py, "
            f"but the actual on-disk source script is {source_script.name}"
        )

    if paths["preprocessed_hvg_h5ad"].exists() and not args.force_refresh and not args.force_rebuild_hvg:
        print(f"[cache] loading existing sweep-local HVG cache: {paths['preprocessed_hvg_h5ad']}")
        adata_hvg = sc.read_h5ad(paths["preprocessed_hvg_h5ad"])
    elif args.reuse_local_hvg_cache and source_hvg_cache.exists() and not args.force_rebuild_hvg:
        print(f"[cache] reusing existing local HVG cache: {source_hvg_cache}")
        adata_hvg = sc.read_h5ad(source_hvg_cache)
    else:
        print("[prep] building HVG object from the source Case Study 5 script")
        adata_hvg = module.build_preprocessed_hvg()

    module.add_merged_timepoints(adata_hvg)

    if args.discovery_group not in adata_hvg.obs.columns:
        raise KeyError(
            f"Requested discovery group column {args.discovery_group!r} is not present in adata_hvg.obs. "
            f"Available columns include: {list(adata_hvg.obs.columns)[:20]}"
        )

    adata_exp_all = adata_hvg[adata_hvg.obs["cohort"] == "experimental"].copy()
    adata_nat = adata_hvg[adata_hvg.obs["cohort"] == "natural"].copy()
    if adata_exp_all.n_obs == 0:
        raise ValueError("No experimental cells were found after preprocessing.")
    if adata_nat.n_obs == 0:
        raise ValueError("No natural-cohort cells were found after preprocessing.")

    adata_discovery = module.balanced_sample_by_group(
        adata_exp_all,
        groupby=args.discovery_group,
        max_cells_per_group=args.discovery_max_cells_per_sample,
        seed=args.discovery_seed,
    )
    module.add_merged_timepoints(adata_discovery)

    discovery_counts = (
        adata_discovery.obs
        .groupby(
            ["sample_id", "subject", "timepoint_merged", "day_merged_numeric", "broad_cell_type"],
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
        .sort_values(["subject", "day_merged_numeric", "broad_cell_type"])
    )

    adata_discovery_cogaps = module.make_cogaps_ready_pretransposed_dense(adata_discovery)
    adata_nat_target = module.make_cogaps_ready_sparse(adata_nat)
    module.check_matrix("prep_discovery_cogaps_X", adata_discovery_cogaps.X)
    module.check_matrix("prep_natural_target_X", adata_nat_target.X)

    print(f"[write] {paths['preprocessed_hvg_h5ad']}")
    adata_hvg.write_h5ad(paths["preprocessed_hvg_h5ad"])
    print(f"[write] {paths['experimental_discovery_h5ad']}")
    adata_discovery.write_h5ad(paths["experimental_discovery_h5ad"])
    print(f"[write] {paths['cogaps_input_h5ad']}")
    adata_discovery_cogaps.write_h5ad(paths["cogaps_input_h5ad"])
    print(f"[write] {paths['natural_target_h5ad']}")
    adata_nat_target.write_h5ad(paths["natural_target_h5ad"])
    discovery_counts.to_csv(paths["discovery_counts_csv"], index=False)
    print(f"[write] {paths['discovery_counts_csv']}")

    manifest = {
        "created_at_utc": utc_now_iso(),
        "source_script_path": str(source_script),
        "source_script_name": source_script.name,
        "source_note": (
            "Task docs reference gse154386_sparse_nondistributed_cogaps.py. "
            "The actual on-disk source script was auto-detected and recorded here."
        ),
        "prep_parameters": expected_signature,
        "local_case_study_cache_paths": {
            "source_hvg_cache": str(source_hvg_cache),
            "source_hvg_cache_exists": bool(source_hvg_cache.exists()),
        },
        "cached_outputs": {key: str(value) for key, value in paths.items() if key not in {"outdir", "cache_dir"}},
        "shapes": {
            "preprocessed_hvg": list(adata_hvg.shape),
            "experimental_discovery_cells_x_genes": list(adata_discovery.shape),
            "experimental_discovery_genes_x_cells": list(adata_discovery_cogaps.shape),
            "natural_target_cells_x_genes": list(adata_nat_target.shape),
        },
        "columns": {
            "preprocessed_hvg_obs": [str(col) for col in adata_hvg.obs.columns],
            "preprocessed_hvg_var": [str(col) for col in adata_hvg.var.columns],
            "experimental_discovery_obs": [str(col) for col in adata_discovery.obs.columns],
            "experimental_discovery_var": [str(col) for col in adata_discovery.var.columns],
            "cogaps_input_obs": [str(col) for col in adata_discovery_cogaps.obs.columns],
            "cogaps_input_var": [str(col) for col in adata_discovery_cogaps.var.columns],
        },
    }
    write_json(paths["prep_manifest_json"], manifest)
    print(f"[write] {paths['prep_manifest_json']}")
    print("[prep] cache preparation complete")


if __name__ == "__main__":
    main()
