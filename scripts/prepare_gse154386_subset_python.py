#!/usr/bin/env python3
"""Download, preprocess, and subset GSE154386 for Case Study 5.

This is an optional full-reproduction script. It prepares the large local
objects used before CoGAPS is run, but it does not run CoGAPS.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import tarfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import requests
import scanpy as sc
from anndata import AnnData
from scipy import sparse
from scipy.io import mmread


RAW_TAR_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE154386&format=file"
RANDOM_SEED = 42

MIN_GENES = 200
MIN_COUNTS = 500
MAX_PCT_MT = 20.0
MIN_CELLS_PER_GENE = 10
N_HVGS = 5000
DISCOVERY_MAX_CELLS_PER_SAMPLE = 600

BROAD_MARKERS = {
    "Monocyte": ["LST1", "FCN1", "CTSS", "SAT1", "TYMP", "S100A8", "S100A9"],
    "T_cell": ["CD3D", "CD3E", "TRAC", "LTB", "IL7R", "MALAT1"],
    "NK_cell": ["NKG7", "GNLY", "PRF1", "CTSW", "KLRD1", "TYROBP"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "BANK1"],
    "Plasmablast": ["MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"],
    "Dendritic": ["FCER1A", "CST3", "CLEC10A", "CD1C", "HLA-DRA"],
    "Neutrophil": ["FCGR3B", "CXCR2", "CSF3R", "MNDA", "S100A8", "S100A9"],
}

PROGRAMS = {
    "ifn_program": [
        "IFITM1", "IFI6", "ISG15", "MX1", "IFIT1", "IFIT3",
        "IFI44L", "ISG20", "LY6E", "TRIM22", "OAS1", "OASL",
    ],
    "plasmablast_program": ["MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"],
    "translation_program": ["RPL4", "RPL5", "RPL6", "RPS3", "RPS8", "EEF2"],
    "mito_program": ["MT-CYB", "MT-ND4", "MT-ND1", "MT-CO1", "MT-CO2"],
}


def download(url: str, dest: Path, chunk_mb: int = 16) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest} already exists")
        return
    print(f"[download] {url} -> {dest}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=180) as response:
        response.raise_for_status()
        with open(tmp, "wb") as handle:
            for chunk in response.iter_content(chunk_size=chunk_mb * 1024 * 1024):
                if chunk:
                    handle.write(chunk)
    tmp.rename(dest)


def safe_extract_tar(tar_path: Path, out_dir: Path) -> None:
    if any(out_dir.iterdir()):
        print(f"[skip] {out_dir} already has extracted files")
        return
    print(f"[extract] {tar_path} -> {out_dir}")
    out_dir_resolved = out_dir.resolve()
    with tarfile.open(tar_path, "r") as tf:
        for member in tf.getmembers():
            target = (out_dir / member.name).resolve()
            if not str(target).startswith(str(out_dir_resolved)):
                raise RuntimeError(f"Unsafe tar member path: {member.name}")
        tf.extractall(out_dir)


def parse_sample_key(sample_key: str) -> Dict[str, object]:
    match = re.match(r"^(GSM\d+)_([^_]+)_(D(?:neg)?\d+)$", sample_key)
    if not match:
        raise ValueError(f"Could not parse sample key: {sample_key}")

    gsm_id, subject, tp = match.groups()
    cohort = "experimental" if subject.startswith("Subject") else "natural"
    day = -int(tp.replace("Dneg", "")) if tp.startswith("Dneg") else int(tp.replace("D", ""))

    return {
        "gsm_id": gsm_id,
        "sample_key": sample_key,
        "sample_id": gsm_id,
        "subject": subject,
        "cohort": cohort,
        "timepoint_raw": tp,
        "timepoint_display": f"D{day}",
        "day_numeric": day,
        "is_reference": (cohort == "experimental" and day == 0)
        or (cohort == "natural" and day == 180),
    }


def add_merged_timepoints(adata: AnnData) -> None:
    labels: List[str] = []
    days: List[float] = []
    for cohort, day, display in zip(
        adata.obs["cohort"].astype(str),
        adata.obs["day_numeric"].astype(float),
        adata.obs["timepoint_display"].astype(str),
    ):
        if cohort == "experimental" and day in (14.0, 15.0):
            labels.append("D14/15")
            days.append(14.5)
        else:
            labels.append(display)
            days.append(float(day))
    adata.obs["timepoint_merged"] = labels
    adata.obs["day_merged_numeric"] = days


def strip_known_suffix(fname: str) -> Tuple[str, Optional[str]]:
    suffix_map = {
        "_matrix.mtx.gz": "matrix",
        "_matrix.mtx": "matrix",
        "_barcodes.tsv.gz": "barcodes",
        "_barcodes.tsv": "barcodes",
        "_features.tsv.gz": "features",
        "_features.tsv": "features",
        "_genes.tsv.gz": "genes",
        "_genes.tsv": "genes",
    }
    for suffix, role in suffix_map.items():
        if fname.endswith(suffix):
            return fname[: -len(suffix)], role
    return fname, None


def discover_10x_triplets(root: Path) -> List[Dict[str, Path]]:
    grouped: Dict[str, Dict[str, Path]] = {}
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        key, role = strip_known_suffix(path.name)
        if role is None:
            continue
        grouped.setdefault(key, {})[role] = path

    triplets: List[Dict[str, Path]] = []
    for key, files in sorted(grouped.items()):
        if "matrix" in files and "barcodes" in files and ("features" in files or "genes" in files):
            entry = {"sample_key": key, **files}
            triplets.append(entry)
    return triplets


def read_text_table(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t", header=None)


def read_mtx(path: Path) -> sparse.csr_matrix:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        return sparse.csr_matrix(mmread(handle))


def build_anndata_from_triplet(triplet: Dict[str, Path]) -> AnnData:
    sample_key = str(triplet["sample_key"])
    meta = parse_sample_key(sample_key)

    matrix = read_mtx(Path(triplet["matrix"]))  # features x cells
    barcodes = read_text_table(Path(triplet["barcodes"])).iloc[:, 0].astype(str).tolist()
    feature_path = Path(triplet.get("features", triplet.get("genes")))
    features = read_text_table(feature_path)

    gene_ids = features.iloc[:, 0].astype(str).tolist()
    gene_symbols = features.iloc[:, 1].astype(str).tolist() if features.shape[1] > 1 else gene_ids
    feature_types = (
        features.iloc[:, 2].astype(str).tolist()
        if features.shape[1] > 2
        else ["Gene Expression"] * len(gene_symbols)
    )
    keep = np.array([x == "Gene Expression" for x in feature_types])
    if keep.sum() == 0:
        keep = np.ones(len(gene_symbols), dtype=bool)

    matrix = matrix[keep, :]
    gene_ids = [g for g, k in zip(gene_ids, keep) if k]
    gene_symbols = [g for g, k in zip(gene_symbols, keep) if k]

    adata = AnnData(X=matrix.T.tocsr())  # cells x genes
    adata.obs_names = [f"{meta['sample_id']}:{barcode}" for barcode in barcodes]
    adata.var_names = pd.Index(gene_symbols)
    adata.var_names_make_unique()
    adata.var["gene_id"] = gene_ids
    adata.var["gene_symbol"] = gene_symbols
    adata.obs["barcode"] = barcodes
    for key, value in meta.items():
        adata.obs[key] = value
    adata.layers["counts"] = adata.X.copy()
    return adata


def flag_upper_outliers_by_sample(adata: AnnData, column: str, q: float = 99.5) -> pd.Series:
    flags = pd.Series(False, index=adata.obs_names)
    for _, idx in adata.obs.groupby("sample_id", observed=True).groups.items():
        vals = adata.obs.loc[idx, column].astype(float).to_numpy()
        if len(vals) >= 10:
            threshold = float(np.nanpercentile(vals, q))
            flags.loc[idx] = adata.obs.loc[idx, column] > threshold
    return flags


def score_gene_sets(adata: AnnData, gene_sets: Dict[str, List[str]]) -> None:
    for label, genes in gene_sets.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=f"{label}_score", use_raw=False)


def annotate_broad_cell_types(adata: AnnData, cluster_key: str = "leiden") -> pd.DataFrame:
    score_cols = [f"{label}_score" for label in BROAD_MARKERS if f"{label}_score" in adata.obs]
    cluster_means = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
    best = cluster_means.idxmax(axis=1).str.replace("_score", "", regex=False)
    adata.obs["broad_cell_type"] = adata.obs[cluster_key].map(best.to_dict()).astype("category")
    out = cluster_means.copy()
    out["assigned_label"] = best
    return out.reset_index()


def balanced_sample_by_group(
    adata: AnnData,
    groupby: str = "sample_id",
    max_cells_per_group: Optional[int] = DISCOVERY_MAX_CELLS_PER_SAMPLE,
    seed: int = RANDOM_SEED,
) -> AnnData:
    if max_cells_per_group is None:
        return adata.copy()
    rng = np.random.default_rng(seed)
    keep: List[str] = []
    for _, idx in adata.obs.groupby(groupby, observed=True).groups.items():
        idx_list = list(idx)
        if len(idx_list) <= max_cells_per_group:
            keep.extend(idx_list)
        else:
            keep.extend(rng.choice(idx_list, max_cells_per_group, replace=False).tolist())
    return adata[keep].copy()


def make_cogaps_ready(adata: AnnData) -> AnnData:
    out = adata.copy()
    x = out.X.toarray() if sparse.issparse(out.X) else np.asarray(out.X)
    x = np.nan_to_num(x.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    x[x < 0] = 0.0
    out.X = x
    return out


def write_manifest(out_path: Path, payload: Dict[str, object]) -> None:
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=Path("data/external/GSE154386"))
    parser.add_argument("--max-cells-per-sample", type=int, default=DISCOVERY_MAX_CELLS_PER_SAMPLE)
    parser.add_argument("--n-hvgs", type=int, default=N_HVGS)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    args = parser.parse_args()

    workdir = args.workdir
    raw_tar = workdir / "GSE154386_RAW.tar"
    extract_dir = workdir / "extracted"
    results_dir = workdir / "results"
    for path in [workdir, extract_dir, results_dir]:
        path.mkdir(parents=True, exist_ok=True)

    download(RAW_TAR_URL, raw_tar)
    safe_extract_tar(raw_tar, extract_dir)

    triplets = discover_10x_triplets(extract_dir)
    if not triplets:
        raise RuntimeError(f"No 10x matrix triplets found under {extract_dir}")

    adatas = [build_anndata_from_triplet(triplet) for triplet in triplets]
    adata_all = sc.concat(
        adatas,
        join="outer",
        label="batch",
        keys=[a.obs["sample_id"].iloc[0] for a in adatas],
        fill_value=0,
        index_unique=None,
    )
    adata_all.layers["counts"] = adata_all.X.copy()
    add_merged_timepoints(adata_all)

    adata_all.var["mt"] = adata_all.var_names.str.upper().str.startswith("MT-")
    adata_all.var["ribo"] = adata_all.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata_all.var["hb"] = adata_all.var_names.str.upper().isin(
        ["HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"]
    )
    sc.pp.calculate_qc_metrics(adata_all, qc_vars=["mt", "ribo", "hb"], log1p=False, inplace=True)
    adata_all.obs["high_genes_outlier"] = flag_upper_outliers_by_sample(adata_all, "n_genes_by_counts")
    adata_all.obs["high_counts_outlier"] = flag_upper_outliers_by_sample(adata_all, "total_counts")

    keep_cells = (
        (adata_all.obs["n_genes_by_counts"] >= MIN_GENES)
        & (adata_all.obs["total_counts"] >= MIN_COUNTS)
        & (adata_all.obs["pct_counts_mt"] <= MAX_PCT_MT)
        & (~adata_all.obs["high_genes_outlier"])
        & (~adata_all.obs["high_counts_outlier"])
    )
    adata_qc = adata_all[keep_cells].copy()
    sc.pp.filter_genes(adata_qc, min_cells=MIN_CELLS_PER_GENE)
    adata_qc.layers["counts"] = adata_qc.X.copy()

    sc.pp.normalize_total(adata_qc, target_sum=1e4)
    sc.pp.log1p(adata_qc)
    add_merged_timepoints(adata_qc)
    score_gene_sets(adata_qc, BROAD_MARKERS)
    score_gene_sets(adata_qc, PROGRAMS)

    try:
        sc.pp.highly_variable_genes(
            adata_qc,
            n_top_genes=args.n_hvgs,
            flavor="seurat_v3",
            layer="counts",
            batch_key="sample_id",
        )
    except Exception:
        sc.pp.highly_variable_genes(adata_qc, n_top_genes=args.n_hvgs, flavor="cell_ranger", layer="counts")

    adata_hvg = adata_qc[:, adata_qc.var["highly_variable"]].copy()
    add_merged_timepoints(adata_hvg)
    sc.pp.pca(adata_hvg, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_hvg, min_dist=0.3)
    sc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden")
    cluster_annotation = annotate_broad_cell_types(adata_hvg)
    cluster_annotation.to_csv(results_dir / "cluster_to_broad_cell_type.csv", index=False)

    comp = (
        adata_hvg.obs.groupby(["cohort", "timepoint_merged", "day_merged_numeric", "broad_cell_type"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    comp["fraction"] = comp.groupby(["cohort", "timepoint_merged"], observed=True)["n_cells"].transform(lambda x: x / x.sum())
    comp.to_csv(results_dir / "cell_composition_by_timepoint.csv", index=False)

    adata_exp = adata_hvg[adata_hvg.obs["cohort"] == "experimental"].copy()
    adata_nat = adata_hvg[adata_hvg.obs["cohort"] == "natural"].copy()
    adata_discovery = balanced_sample_by_group(
        adata_exp,
        groupby="sample_id",
        max_cells_per_group=args.max_cells_per_sample,
        seed=args.seed,
    )
    add_merged_timepoints(adata_discovery)

    discovery_counts = (
        adata_discovery.obs.groupby(
            ["sample_id", "subject", "timepoint_merged", "day_merged_numeric", "broad_cell_type"],
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
    )
    discovery_counts.to_csv(results_dir / "experimental_discovery_cell_counts.csv", index=False)

    adata_discovery_cogaps = make_cogaps_ready(adata_discovery)
    adata_nat_target = make_cogaps_ready(adata_nat)

    adata_hvg.write_h5ad(workdir / "gse154386_preprocessed_hvg.h5ad")
    adata_discovery_cogaps.write_h5ad(workdir / "gse154386_experimental_discovery_cells_x_genes.h5ad")
    adata_discovery_cogaps.T.write_h5ad(workdir / "gse154386_experimental_discovery_genes_x_cells.h5ad")
    adata_nat_target.write_h5ad(workdir / "gse154386_natural_projection_target.h5ad")

    write_manifest(
        results_dir / "preprocessing_manifest.json",
        {
            "source": "GSE154386",
            "raw_tar_url": RAW_TAR_URL,
            "n_samples": len(triplets),
            "all_cells_after_qc": int(adata_qc.n_obs),
            "hvg_shape": [int(adata_hvg.n_obs), int(adata_hvg.n_vars)],
            "experimental_discovery_shape": [int(adata_discovery_cogaps.n_obs), int(adata_discovery_cogaps.n_vars)],
            "natural_projection_shape": [int(adata_nat_target.n_obs), int(adata_nat_target.n_vars)],
            "max_cells_per_sample": args.max_cells_per_sample,
            "seed": args.seed,
            "outputs": {
                "preprocessed_hvg": str(workdir / "gse154386_preprocessed_hvg.h5ad"),
                "discovery_cells_x_genes": str(workdir / "gse154386_experimental_discovery_cells_x_genes.h5ad"),
                "discovery_genes_x_cells": str(workdir / "gse154386_experimental_discovery_genes_x_cells.h5ad"),
                "natural_projection_target": str(workdir / "gse154386_natural_projection_target.h5ad"),
            },
        },
    )


if __name__ == "__main__":
    main()
