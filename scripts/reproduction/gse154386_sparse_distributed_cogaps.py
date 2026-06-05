# =========================
# GSE154386: NON-DISTRIBUTED PyCoGAPS workflow for temporal dengue case study
#
# Main workflow
#   1) Build or load preprocessed HVG AnnData
#   2) Split experimental / natural cohorts
#   3) Build a balanced experimental discovery set
#   4) Pretranspose discovery input to genes x cells for CoGAPS
#   5) Run NON-DISTRIBUTED PyCoGAPS
#   6) Analyze CoGAPS patterns for Q1-Q3
#   7) Project learned gene weights into natural cohort for Q4
#
# Key fixes
#   - NO distributed mode
#   - discovery input is PRETRANSPOSED to genes x cells so transposeData=False works
#   - discovery input is dense float32 because your installed PyCoGAPS constructor
#     accepts numpy.ndarray[numpy.float32] and rejected SciPy sparse matrices
#   - CoGAPS gene/cell pattern matrices are cleaned for NaN/Inf before downstream use
#   - NNLS projection uses cleaned nonnegative finite gene loadings
#
# Outputs
#   - GSE154386/GSE154386_all_preprocessed_HVG.h5ad
#   - GSE154386/GSE154386_experimental_discovery_cogaps_input_nondist.h5ad
#   - GSE154386/GSE154386_experimental_discovery_cogaps_result_nondist.h5ad
#   - GSE154386/GSE154386_natural_projection_target_nondist.h5ad
#   - GSE154386/results/*
#
# Recommended run command:
#   python -u gse154386_nondistributed_cogaps.py | tee GSE154386/results/cogaps_nondist_run.log
# =========================

from __future__ import annotations

import gzip
import os
import re
import sys
import tarfile
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import scanpy as sc
from anndata import AnnData
from scipy import sparse
from scipy.io import mmread
from scipy.optimize import nnls
from scipy.stats import kruskal, spearmanr

# -------------------------
# Make stdout less annoying in terminals
# -------------------------
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass

# -------------------------
# PyCoGAPS is REQUIRED
# -------------------------
try:
    from PyCoGAPS.parameters import CoParams, setParams
    from PyCoGAPS.pycogaps_main import CoGAPS
except ImportError as e:
    raise ImportError(
        "PyCoGAPS is required for this script. Activate the correct environment first, "
        "for example: source ~/PycharmProjects/pycogaps/.venv/bin/activate"
    ) from e

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'.*")

# -------------------------
# 0) Global config
# -------------------------
WORKDIR = Path("GSE154386")
WORKDIR.mkdir(exist_ok=True, parents=True)

RAW_TAR = WORKDIR / "GSE154386_RAW.tar"
EXTRACT_DIR = WORKDIR / "extracted"
EXTRACT_DIR.mkdir(exist_ok=True)

RAW_TAR_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE154386&format=file"

RESULTS_DIR = WORKDIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

COGAPS_DIR = RESULTS_DIR / "cogaps_nondist"
COGAPS_DIR.mkdir(exist_ok=True)

OUT_ALL_HVG = WORKDIR / "GSE154386_all_preprocessed_HVG.h5ad"
OUT_DISCOVERY_COGAPS = WORKDIR / "GSE154386_experimental_discovery_cogaps_input_nondist.h5ad"
OUT_COGAPS_RESULT = WORKDIR / "GSE154386_experimental_discovery_cogaps_result_nondist.h5ad"
OUT_NAT_TARGET = WORKDIR / "GSE154386_natural_projection_target_nondist.h5ad"

# Toggle: reuse the HVG object you already built
USE_CACHED_HVG = True

# QC / preprocessing
MIN_GENES = 200
MIN_COUNTS = 500
MAX_PCT_MT = 20.0
MIN_CELLS_PER_GENE = 10
N_HVGS = 5000
RANDOM_SEED = 42

# Discovery-set construction for the actual CoGAPS run
# Lower than your distributed run because non-distributed is slower.
DISCOVERY_MAX_CELLS_PER_SAMPLE = 200
DISCOVERY_SAMPLE_GROUP = "sample_id"

# CoGAPS parameters
COGAPS_N_PATTERNS = 12
COGAPS_N_ITER = 4000
COGAPS_N_THREADS = max(1, min(8, os.cpu_count() or 4))
COGAPS_USE_SPARSE_OPT = True

# Projection into natural cohort
PROJECTION_TOP_GENES_PER_PATTERN = 250
PROJECTION_PROGRESS_EVERY = 500
MIN_COMMON_PROJECTION_GENES = 100

# Marker sets for annotation and interpretation
BROAD_MARKERS = {
    "Monocyte":    ["LST1", "FCN1", "CTSS", "SAT1", "TYMP", "S100A8", "S100A9"],
    "T_cell":      ["CD3D", "CD3E", "TRAC", "LTB", "IL7R", "MALAT1"],
    "NK_cell":     ["NKG7", "GNLY", "PRF1", "CTSW", "KLRD1", "TYROBP"],
    "B_cell":      ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "BANK1"],
    "Plasmablast": ["MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"],
    "Dendritic":   ["FCER1A", "CST3", "CLEC10A", "CD1C", "HLA-DRA"],
    "Neutrophil":  ["FCGR3B", "CXCR2", "CSF3R", "MNDA", "S100A8", "S100A9"],
}

PROGRAMS = {
    "ifn_program": [
        "IFITM1", "IFI6", "ISG15", "MX1", "IFIT1", "IFIT3",
        "IFI44L", "ISG20", "LY6E", "TRIM22", "OAS1", "OASL"
    ],
    "translation_program": [
        "RPL4", "RPL5", "RPL6", "RPS3", "RPS8", "EEF2", "EIF3L", "EIF4B"
    ],
    "mito_program": [
        "MT-CYB", "MT-ND4", "MT-ND1", "MT-CO1", "MT-CO2"
    ],
    "plasmablast_program": [
        "MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"
    ],
}

# -------------------------
# 1) Download + extract
# -------------------------
def download(url: str, dest: Path, chunk_mb: int = 16) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest} exists ({dest.stat().st_size/1e6:.1f} MB)", flush=True)
        return
    print(f"[download] {url}\n          -> {dest}", flush=True)
    with requests.get(url, stream=True, timeout=180) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_mb * 1024 * 1024):
                if chunk:
                    f.write(chunk)
        tmp.rename(dest)
    print("[done]", flush=True)

def extract_tar(tar_path: Path, out_dir: Path) -> None:
    if any(out_dir.iterdir()):
        print(f"[skip] {out_dir} already has files", flush=True)
        return
    print(f"[extract] {tar_path} -> {out_dir}", flush=True)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(out_dir)
    print("[done]", flush=True)

# -------------------------
# 2) GEO sample parsing
# -------------------------
def parse_sample_key(sample_key: str) -> Dict[str, object]:
    meta: Dict[str, object] = {
        "gsm_id": None,
        "sample_key": sample_key,
        "sample_id": None,
        "subject": None,
        "cohort": None,
        "timepoint_raw": None,
        "timepoint_display": None,
        "day_numeric": None,
        "is_reference": False,
        "reference_label": None,
    }

    m = re.match(r"^(GSM\d+)_([^_]+)_(D(?:neg)?\d+)$", sample_key)
    if not m:
        raise ValueError(f"Could not parse sample key: {sample_key}")

    gsm_id, subject, tp = m.groups()
    meta["gsm_id"] = gsm_id
    meta["sample_id"] = gsm_id
    meta["subject"] = subject
    meta["timepoint_raw"] = tp

    if subject.startswith("Subject"):
        meta["cohort"] = "experimental"
    elif subject.startswith("Natural"):
        meta["cohort"] = "natural"
    else:
        meta["cohort"] = "unknown"

    if tp.startswith("Dneg"):
        day = -int(tp.replace("Dneg", ""))
        display = f"D{day}"
    else:
        day = int(tp.replace("D", ""))
        display = f"D{day}"

    meta["day_numeric"] = day
    meta["timepoint_display"] = display

    if meta["cohort"] == "experimental" and day == 0:
        meta["is_reference"] = True
        meta["reference_label"] = "experimental_baseline"
    elif meta["cohort"] == "natural" and day == 180:
        meta["is_reference"] = True
        meta["reference_label"] = "natural_convalescent"
    else:
        meta["is_reference"] = False
        meta["reference_label"] = "not_reference"

    return meta

def add_merged_timepoints(adata: AnnData) -> None:
    merged_labels = []
    merged_days = []

    for cohort, raw_day, disp in zip(
        adata.obs["cohort"].astype(str),
        adata.obs["day_numeric"].astype(float),
        adata.obs["timepoint_display"].astype(str),
    ):
        if cohort == "experimental" and raw_day in (14.0, 15.0):
            merged_labels.append("D14/15")
            merged_days.append(14.5)
        else:
            merged_labels.append(disp)
            merged_days.append(float(raw_day))

    adata.obs["timepoint_merged"] = merged_labels
    adata.obs["day_merged_numeric"] = merged_days

    display_order = (
        adata.obs[["timepoint_display", "day_numeric"]]
        .drop_duplicates()
        .sort_values("day_numeric")["timepoint_display"]
        .tolist()
    )
    merged_order = (
        adata.obs[["timepoint_merged", "day_merged_numeric"]]
        .drop_duplicates()
        .sort_values("day_merged_numeric")["timepoint_merged"]
        .tolist()
    )

    adata.obs["timepoint_display"] = pd.Categorical(
        adata.obs["timepoint_display"],
        categories=display_order,
        ordered=True,
    )
    adata.obs["timepoint_merged"] = pd.Categorical(
        adata.obs["timepoint_merged"],
        categories=merged_order,
        ordered=True,
    )

# -------------------------
# 3) Discover and read 10x triplets
# -------------------------
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
    for suf, role in suffix_map.items():
        if fname.endswith(suf):
            return fname[:-len(suf)], role
    return fname, None

def discover_10x_triplets(root: Path) -> List[Dict[str, object]]:
    groups: Dict[str, Dict[str, Path]] = {}

    for p in root.rglob("*"):
        if not p.is_file():
            continue
        key, role = strip_known_suffix(p.name)
        if role is None:
            continue
        groups.setdefault(key, {})
        groups[key][role] = p

    triplets: List[Dict[str, object]] = []
    for key, d in sorted(groups.items()):
        if "matrix" not in d or "barcodes" not in d or ("features" not in d and "genes" not in d):
            continue
        entry: Dict[str, object] = {"sample_key": key}
        entry.update(d)
        triplets.append(entry)

    return triplets

def read_text_table(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return pd.read_csv(f, sep="\t", header=None)

def read_mtx(path: Path) -> sparse.csr_matrix:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as f:
        mat = mmread(f)
    return sparse.csr_matrix(mat)

def build_anndata_from_triplet(triplet: Dict[str, object]) -> AnnData:
    sample_key = str(triplet["sample_key"])
    meta = parse_sample_key(sample_key)

    matrix = read_mtx(Path(triplet["matrix"]))  # features x cells
    barcodes = read_text_table(Path(triplet["barcodes"])).iloc[:, 0].astype(str).tolist()

    feature_path = Path(triplet["features"]) if "features" in triplet else Path(triplet["genes"])
    feat = read_text_table(feature_path)

    if feat.shape[1] == 1:
        gene_ids = feat.iloc[:, 0].astype(str).tolist()
        gene_symbols = feat.iloc[:, 0].astype(str).tolist()
        feature_types = ["Gene Expression"] * len(gene_symbols)
    elif feat.shape[1] == 2:
        gene_ids = feat.iloc[:, 0].astype(str).tolist()
        gene_symbols = feat.iloc[:, 1].astype(str).tolist()
        feature_types = ["Gene Expression"] * len(gene_symbols)
    else:
        gene_ids = feat.iloc[:, 0].astype(str).tolist()
        gene_symbols = feat.iloc[:, 1].astype(str).tolist()
        feature_types = feat.iloc[:, 2].astype(str).tolist()

    if matrix.shape[1] != len(barcodes):
        raise ValueError(
            f"Barcode length mismatch for {sample_key}: matrix has {matrix.shape[1]} columns, "
            f"barcodes file has {len(barcodes)} rows"
        )
    if matrix.shape[0] != len(gene_symbols):
        raise ValueError(
            f"Feature length mismatch for {sample_key}: matrix has {matrix.shape[0]} rows, "
            f"features file has {len(gene_symbols)} rows"
        )

    keep = np.array([ft == "Gene Expression" for ft in feature_types], dtype=bool)
    if keep.sum() == 0:
        keep = np.ones(len(gene_symbols), dtype=bool)

    matrix = matrix[keep, :]
    gene_ids = [g for g, k in zip(gene_ids, keep) if k]
    gene_symbols = [g for g, k in zip(gene_symbols, keep) if k]
    feature_types = [g for g, k in zip(feature_types, keep) if k]

    X = matrix.T.tocsr()  # cells x genes

    ad = AnnData(X=X)
    ad.obs_names = [f"{meta['sample_id']}:{bc}" for bc in barcodes]
    ad.var_names = pd.Index([str(x) for x in gene_symbols])
    ad.var_names_make_unique()

    ad.var["gene_id"] = [str(x) for x in gene_ids]
    ad.var["gene_symbol"] = [str(x) for x in gene_symbols]
    ad.var["feature_type"] = [str(x) for x in feature_types]

    ad.obs["barcode"] = barcodes
    for k, v in meta.items():
        ad.obs[k] = v

    ad.layers["counts"] = ad.X.copy()
    return ad

# -------------------------
# 4) Validation + utilities
# -------------------------
def _min_max_sparse(X: sparse.spmatrix) -> Tuple[float, float]:
    if X.nnz == 0:
        return 0.0, 0.0
    return float(X.min()), float(X.max())

def check_matrix(name: str, X) -> None:
    if sparse.issparse(X):
        mn, mx = _min_max_sparse(X)
        data = X.data
        n_nan = int(np.isnan(data).sum())
        n_inf = int(np.isinf(data).sum())
        n_neg = int((data < 0).sum())
        density = float(X.nnz) / float(X.shape[0] * X.shape[1])
        print(
            f"[check:{name}] sparse shape={X.shape} density={density:.5f} "
            f"min={mn} max={mx} NaN={n_nan} Inf={n_inf} Neg={n_neg}",
            flush=True
        )
    else:
        arr = np.asarray(X)
        n_nan = int(np.isnan(arr).sum())
        n_inf = int(np.isinf(arr).sum())
        n_neg = int((arr < 0).sum())
        density = float(np.count_nonzero(arr)) / float(arr.size)
        print(
            f"[check:{name}] dense shape={arr.shape} density={density:.5f} "
            f"min={float(arr.min())} max={float(arr.max())} "
            f"NaN={n_nan} Inf={n_inf} Neg={n_neg}",
            flush=True
        )

def drop_all_zero_genes(adata: AnnData) -> AnnData:
    Xc = adata.layers["counts"]
    gene_sums = np.asarray(Xc.sum(axis=0)).ravel()
    keep = gene_sums > 0
    dropped = int((~keep).sum())
    if dropped > 0:
        print(f"[filter] dropping {dropped} all-zero genes", flush=True)
        adata = adata[:, keep].copy()
    return adata

def flag_upper_outliers_by_sample(adata: AnnData, column: str, q: float = 99.5) -> pd.Series:
    flags = pd.Series(False, index=adata.obs_names, dtype=bool)
    for sample_id, idx in adata.obs.groupby("sample_id", observed=True).groups.items():
        vals = adata.obs.loc[idx, column].astype(float).to_numpy()
        if len(vals) < 10:
            continue
        threshold = float(np.nanpercentile(vals, q))
        flags.loc[idx] = adata.obs.loc[idx, column] > threshold
        print(f"[qc] sample={sample_id} {column} {q}th percentile threshold={threshold:.2f}", flush=True)
    return flags

def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)
    mask = np.isfinite(pvals)
    if mask.sum() == 0:
        return out

    p = pvals[mask]
    order = np.argsort(p)
    ranked = p[order]
    m = len(ranked)

    adj_ranked = ranked * m / np.arange(1, m + 1)
    adj_ranked = np.minimum.accumulate(adj_ranked[::-1])[::-1]
    adj_ranked = np.minimum(adj_ranked, 1.0)

    adj = np.empty_like(adj_ranked)
    adj[order] = adj_ranked
    out[mask] = adj
    return out

def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)

    mask = np.isfinite(values)
    values = values[mask]
    groups = groups[mask]

    if len(values) == 0:
        return np.nan

    grand_mean = values.mean()
    ss_total = np.sum((values - grand_mean) ** 2)
    if ss_total <= 0:
        return 0.0

    ss_between = 0.0
    for g in pd.unique(groups):
        vals = values[groups == g]
        if len(vals) == 0:
            continue
        ss_between += len(vals) * (vals.mean() - grand_mean) ** 2

    return float(ss_between / ss_total)

def classify_pattern(eta_time: float, eta_cell: float) -> str:
    if not np.isfinite(eta_time) or not np.isfinite(eta_cell):
        return "unclassified"
    if eta_cell > 1.5 * eta_time:
        return "identity-like"
    if eta_time > 1.5 * eta_cell:
        return "activity-like"
    return "mixed"

def write_text(path: Path, text: str) -> None:
    with open(path, "w") as f:
        f.write(text)

def get_timepoint_order(df: pd.DataFrame, label_col: str, day_col: str) -> List[str]:
    return (
        df[[label_col, day_col]]
        .drop_duplicates()
        .sort_values(day_col)[label_col]
        .tolist()
    )

def numeric_suffix(s: str) -> int:
    m = re.findall(r"(\d+)", str(s))
    return int(m[0]) if m else 999999

def get_pattern_columns(cols: List[str]) -> List[str]:
    pats = [str(c) for c in cols if re.fullmatch(r"Pattern_?\d+", str(c))]
    return sorted(pats, key=numeric_suffix)

def clean_sparse_nonnegative(X) -> sparse.csr_matrix:
    if not sparse.issparse(X):
        X = sparse.csr_matrix(np.asarray(X))
    else:
        X = X.tocsr().copy()
    X.data = np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    X.data[X.data < 0] = 0.0
    X.eliminate_zeros()
    return X

# FIX: pretranspose discovery data into dense genes x cells for PyCoGAPS
def make_cogaps_ready_pretransposed_dense(adata_subset: AnnData) -> AnnData:
    """
    Build a dense float32 AnnData in genes x cells orientation so PyCoGAPS can be run
    with transposeData=False in installations that do not accept sparse adata.X in CoParams.
    """
    lognorm = adata_subset.X.copy()

    if sparse.issparse(lognorm):
        X = lognorm.toarray()
    else:
        X = np.asarray(lognorm)

    X = np.asarray(X, dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    X[X < 0] = 0.0

    # pretranspose: cells x genes -> genes x cells
    X_gxc = np.ascontiguousarray(X.T, dtype=np.float32)

    out = AnnData(X=X_gxc)

    # genes now live in obs
    out.obs_names = adata_subset.var_names.copy()
    out.obs = adata_subset.var.copy()

    # cells now live in var
    out.var_names = adata_subset.obs_names.copy()
    out.var = adata_subset.obs.copy()

    # keep transposed layers for traceability
    counts = adata_subset.layers["counts"]
    if sparse.issparse(counts):
        out.layers["counts"] = counts.T.tocsr()
    else:
        out.layers["counts"] = sparse.csr_matrix(np.asarray(counts).T)

    if sparse.issparse(lognorm):
        out.layers["lognorm"] = lognorm.T.tocsr()
    else:
        out.layers["lognorm"] = sparse.csr_matrix(np.asarray(lognorm).T)

    return out

def clean_pattern_df(df: pd.DataFrame, clip_negative: bool = True) -> pd.DataFrame:
    arr = df.to_numpy(dtype=np.float64)
    n_nan = int(np.isnan(arr).sum())
    n_inf = int(np.isinf(arr).sum())
    n_neg = int((arr < 0).sum())

    if n_nan or n_inf or n_neg:
        print(
            f"[clean] pattern matrix had NaN={n_nan}, Inf={n_inf}, Neg={n_neg}; replacing with safe values",
            flush=True,
        )

    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    if clip_negative:
        arr[arr < 0] = 0.0

    return pd.DataFrame(arr, index=df.index.copy(), columns=df.columns.copy())

# -------------------------
# 5) Scoring + annotation
# -------------------------
def score_marker_sets(adata: AnnData, gene_sets: Dict[str, List[str]]) -> List[str]:
    score_cols: List[str] = []
    for label, genes in gene_sets.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            print(f"[warn] skipping score {label}; only {len(present)} genes present", flush=True)
            continue
        score_name = f"{label}_score"
        sc.tl.score_genes(adata, gene_list=present, score_name=score_name, use_raw=False)
        score_cols.append(score_name)
        print(f"[score] {label}: {len(present)} genes present", flush=True)
    return score_cols

def annotate_broad_cell_types(adata: AnnData, cluster_key: str = "leiden") -> pd.DataFrame:
    score_cols = [f"{k}_score" for k in BROAD_MARKERS.keys() if f"{k}_score" in adata.obs.columns]
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"{cluster_key} not found in adata.obs")
    if not score_cols:
        raise ValueError("No broad marker score columns found.")

    cluster_means = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
    best_scores = cluster_means.max(axis=1)
    best_labels = cluster_means.idxmax(axis=1).str.replace("_score", "", regex=False)
    best_labels = best_labels.where(best_scores > 0, other="Unassigned")

    mapping = best_labels.to_dict()
    adata.obs["broad_cell_type"] = adata.obs[cluster_key].map(mapping).astype("category")

    out = cluster_means.copy()
    out["assigned_label"] = best_labels
    out["best_score"] = best_scores
    out.index.name = cluster_key
    return out.reset_index()

# -------------------------
# 6) Plot helpers
# -------------------------
def save_umap_panels(adata: AnnData, colors: List[str], out_prefix: str) -> None:
    for color in colors:
        if color not in adata.obs.columns and color not in adata.var.columns:
            continue
        try:
            fig = sc.pl.umap(adata, color=color, show=False, return_fig=True)
            if fig is not None:
                fig.savefig(FIG_DIR / f"{out_prefix}_{color}.png", dpi=200, bbox_inches="tight")
                plt.close(fig)
        except Exception as e:
            print(f"[warn] could not save UMAP for {color}: {e}", flush=True)

def save_composition_plot(comp_df: pd.DataFrame, out_png: Path, label_col: str, day_col: str) -> None:
    if comp_df.empty:
        return

    tmp = (
        comp_df.groupby([label_col, day_col, "broad_cell_type"], as_index=False, observed=True)["fraction"]
        .mean()
    )
    order = get_timepoint_order(tmp, label_col, day_col)
    pivot = tmp.pivot_table(
        index=label_col,
        columns="broad_cell_type",
        values="fraction",
        aggfunc="mean",
        fill_value=0,
    ).reindex(order)

    plt.figure(figsize=(10, 5))
    pivot.plot(kind="bar", stacked=True, ax=plt.gca())
    plt.ylabel("Fraction of cells")
    plt.xlabel("Timepoint")
    plt.title("Cell composition over time")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def save_pattern_lineplot(summary_df: pd.DataFrame, pattern_col: str, out_png: Path, title: str) -> None:
    if summary_df.empty:
        return
    plt.figure(figsize=(10, 5))
    for cell_type, sub in summary_df.groupby("broad_cell_type", observed=True):
        sub = sub.sort_values("day_merged_numeric")
        plt.plot(sub["day_merged_numeric"], sub[pattern_col], marker="o", label=str(cell_type))
    plt.xlabel("Day")
    plt.ylabel(pattern_col)
    plt.title(title)
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def save_heatmap(df: pd.DataFrame, index_col: str, column_col: str, value_col: str, day_col: Optional[str], out_png: Path, title: str) -> None:
    if df.empty:
        return
    pivot = df.pivot_table(
        index=index_col,
        columns=column_col,
        values=value_col,
        aggfunc="mean",
        fill_value=0,
    )
    if day_col is not None and day_col in df.columns:
        order = get_timepoint_order(df[[index_col, day_col]].drop_duplicates(), index_col, day_col)
        pivot = pivot.reindex(order)
    plt.figure(figsize=(10, 6))
    plt.imshow(pivot.values, aspect="auto")
    plt.xticks(range(pivot.shape[1]), pivot.columns, rotation=90)
    plt.yticks(range(pivot.shape[0]), pivot.index)
    plt.colorbar(label=value_col)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

# -------------------------
# 7) Build preprocessed HVG object
# -------------------------
def build_preprocessed_hvg() -> AnnData:
    download(RAW_TAR_URL, RAW_TAR)
    extract_tar(RAW_TAR, EXTRACT_DIR)

    triplets = discover_10x_triplets(EXTRACT_DIR)
    print(f"[files] found {len(triplets)} 10x sample groups", flush=True)
    assert len(triplets) > 0, "No 10x sample groups found."

    adatas: List[AnnData] = []
    for i, triplet in enumerate(triplets, start=1):
        ad = build_anndata_from_triplet(triplet)
        adatas.append(ad)
        if i <= 3:
            print(f"[loaded] {ad.obs['sample_key'][0]} -> {ad.shape}", flush=True)
            print(ad.obs.iloc[0][["gsm_id", "subject", "cohort", "timepoint_display", "day_numeric"]].to_dict(), flush=True)

    adata_all = sc.concat(
        adatas,
        join="outer",
        label="batch",
        keys=[a.obs["sample_id"][0] for a in adatas],
        fill_value=0,
        index_unique=None,
    )
    adata_all.layers["counts"] = adata_all.X.copy()
    add_merged_timepoints(adata_all)

    print("\n[concat]", adata_all, flush=True)
    print("\ncohort counts:\n", adata_all.obs["cohort"].value_counts(dropna=False), flush=True)
    print("\nsubject counts:\n", adata_all.obs["subject"].value_counts(dropna=False), flush=True)
    print("\ntimepoint counts:\n", adata_all.obs["timepoint_display"].value_counts(dropna=False).sort_index(), flush=True)

    check_matrix("counts_pre", adata_all.layers["counts"])
    adata_all = drop_all_zero_genes(adata_all)
    check_matrix("counts_post_zero_gene_drop", adata_all.layers["counts"])

    adata_all.X = adata_all.layers["counts"].copy()
    adata_all.var["mt"] = adata_all.var_names.str.upper().str.startswith("MT-")
    adata_all.var["ribo"] = adata_all.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata_all.var["hb"] = adata_all.var_names.str.upper().isin(
        ["HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"]
    )

    sc.pp.calculate_qc_metrics(
        adata_all,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[20],
        log1p=False,
        inplace=True,
    )

    adata_all.obs["high_genes_outlier"] = flag_upper_outliers_by_sample(adata_all, "n_genes_by_counts", q=99.5)
    adata_all.obs["high_counts_outlier"] = flag_upper_outliers_by_sample(adata_all, "total_counts", q=99.5)

    keep_cells = (
        (adata_all.obs["n_genes_by_counts"] >= MIN_GENES) &
        (adata_all.obs["total_counts"] >= MIN_COUNTS) &
        (adata_all.obs["pct_counts_mt"] <= MAX_PCT_MT) &
        (~adata_all.obs["high_genes_outlier"]) &
        (~adata_all.obs["high_counts_outlier"])
    )

    print(f"\n[qc] keeping {int(keep_cells.sum())} / {adata_all.n_obs} cells", flush=True)
    adata_qc = adata_all[keep_cells].copy()

    sc.pp.filter_genes(adata_qc, min_cells=MIN_CELLS_PER_GENE)
    adata_qc.layers["counts"] = adata_qc.X.copy()
    check_matrix("counts_post_qc", adata_qc.layers["counts"])

    sc.pp.normalize_total(adata_qc, target_sum=1e4)
    sc.pp.log1p(adata_qc)
    add_merged_timepoints(adata_qc)
    check_matrix("X_log1p", adata_qc.X)

    broad_score_cols = score_marker_sets(adata_qc, BROAD_MARKERS)
    program_score_cols = score_marker_sets(adata_qc, PROGRAMS)

    try:
        sc.pp.highly_variable_genes(
            adata_qc,
            n_top_genes=N_HVGS,
            flavor="seurat_v3",
            layer="counts",
            batch_key="sample_id",
        )
        print("[HVG] used flavor='seurat_v3' batch_key='sample_id'", flush=True)
    except Exception as e:
        print(f"[warn] seurat_v3 HVG failed ({e}); falling back to flavor='cell_ranger'", flush=True)
        sc.pp.highly_variable_genes(
            adata_qc,
            n_top_genes=N_HVGS,
            flavor="cell_ranger",
            layer="counts",
        )

    adata_hvg = adata_qc[:, adata_qc.var["highly_variable"]].copy()
    add_merged_timepoints(adata_hvg)

    print("\n[HVG]", adata_hvg, flush=True)
    check_matrix("counts_HVG", adata_hvg.layers["counts"])
    check_matrix("X_HVG", adata_hvg.X)

    sc.pp.pca(adata_hvg, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_hvg, min_dist=0.3)
    sc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden")

    cluster_annotation = annotate_broad_cell_types(adata_hvg, cluster_key="leiden")
    cluster_annotation.to_csv(RESULTS_DIR / "cluster_to_broad_cell_type.csv", index=False)

    save_umap_panels(
        adata_hvg,
        colors=[
            "cohort", "subject", "timepoint_display", "timepoint_merged",
            "leiden", "broad_cell_type"
        ] + [c for c in ["ifn_program_score", "plasmablast_program_score", "translation_program_score"] if c in adata_hvg.obs.columns],
        out_prefix="full_umap",
    )

    comp = (
        adata_hvg.obs
        .groupby(["cohort", "timepoint_merged", "day_merged_numeric", "broad_cell_type"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    comp["fraction"] = (
        comp.groupby(["cohort", "timepoint_merged"], observed=True)["n_cells"]
        .transform(lambda x: x / x.sum())
    )
    comp.to_csv(RESULTS_DIR / "cell_composition_by_timepoint.csv", index=False)

    exp_comp = comp.loc[comp["cohort"] == "experimental"].copy()
    nat_comp = comp.loc[comp["cohort"] == "natural"].copy()

    save_composition_plot(
        exp_comp,
        FIG_DIR / "experimental_cell_composition_over_time.png",
        label_col="timepoint_merged",
        day_col="day_merged_numeric",
    )
    if len(nat_comp) > 0:
        save_composition_plot(
            nat_comp,
            FIG_DIR / "natural_cell_composition_over_time.png",
            label_col="timepoint_merged",
            day_col="day_merged_numeric",
        )

    adata_hvg.write_h5ad(OUT_ALL_HVG)
    print(f"[saved HVG] {OUT_ALL_HVG}", flush=True)
    return adata_hvg

# -------------------------
# 8) Discovery-set construction
# -------------------------
def balanced_sample_by_group(
    adata: AnnData,
    groupby: str,
    max_cells_per_group: Optional[int],
    seed: int = 42,
) -> AnnData:
    if max_cells_per_group is None:
        print("[discovery] using all cells in the experimental cohort", flush=True)
        return adata.copy()

    rng = np.random.default_rng(seed)
    keep_idx: List[str] = []
    for group, idx in adata.obs.groupby(groupby, observed=True).groups.items():
        idx_list = list(idx)
        if len(idx_list) <= max_cells_per_group:
            keep_idx.extend(idx_list)
            kept = len(idx_list)
        else:
            chosen = rng.choice(idx_list, size=max_cells_per_group, replace=False)
            keep_idx.extend(chosen.tolist())
            kept = max_cells_per_group
        print(f"[discovery] {group}: kept {kept} of {len(idx_list)}", flush=True)
    return adata[keep_idx].copy()

def make_cogaps_ready_sparse(adata_subset: AnnData) -> AnnData:
    lognorm = adata_subset.X.copy()
    X = clean_sparse_nonnegative(lognorm)

    out = adata_subset.copy()
    out.X = X
    out.layers["counts"] = adata_subset.layers["counts"].copy()
    out.layers["lognorm"] = lognorm.copy() if sparse.issparse(lognorm) else sparse.csr_matrix(np.asarray(lognorm))
    return out

# -------------------------
# 9) Run NON-DISTRIBUTED CoGAPS
# -------------------------
def run_cogaps_nondistributed(adata_discovery_cogaps: AnnData, out_path: Path) -> AnnData:
    print("\n[CoGAPS] starting NON-DISTRIBUTED run", flush=True)
    print(
        f"[CoGAPS] discovery input shape={adata_discovery_cogaps.shape}, "
        f"nPatterns={COGAPS_N_PATTERNS}, nIterations={COGAPS_N_ITER}, "
        f"useSparseOptimization={COGAPS_USE_SPARSE_OPT}, nThreads={COGAPS_N_THREADS}",
        flush=True,
    )

    params = CoParams(adata=adata_discovery_cogaps)
    setParams(params, {
        "nPatterns": COGAPS_N_PATTERNS,
        "nIterations": COGAPS_N_ITER,
        "seed": RANDOM_SEED,
        "useSparseOptimization": COGAPS_USE_SPARSE_OPT,
        "transposeData": False,   # discovery input has already been pretransposed to genes x cells
        "nThreads": COGAPS_N_THREADS,
        # intentionally NO distributed key
    })

    params.printParams()

    start = time.time()
    result = CoGAPS(adata_discovery_cogaps, params)
    end = time.time()
    print(f"[CoGAPS] runtime: {(end - start)/60:.2f} min", flush=True)

    if not isinstance(result, AnnData):
        raise TypeError(f"Expected PyCoGAPS to return an AnnData result, got {type(result)} instead.")

    result.write_h5ad(out_path)
    print(f"[saved CoGAPS result] {out_path}", flush=True)
    return result

# -------------------------
# 10) CoGAPS result extraction
# -------------------------
def extract_cogaps_matrices(result: AnnData, discovery_input: AnnData) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    result.obs_names = result.obs_names.astype(str)
    result.var_names = result.var_names.astype(str)

    obs_pats = get_pattern_columns(result.obs.columns.tolist())
    var_pats = get_pattern_columns(result.var.columns.tolist())
    if len(obs_pats) == 0 or len(var_pats) == 0:
        raise ValueError("Could not find Pattern columns in the CoGAPS result object.")

    discovery_cells = set(discovery_input.obs_names.astype(str))
    discovery_genes = set(discovery_input.var_names.astype(str))

    obs_cell_overlap = len(set(result.obs_names).intersection(discovery_cells))
    var_cell_overlap = len(set(result.var_names).intersection(discovery_cells))
    obs_gene_overlap = len(set(result.obs_names).intersection(discovery_genes))
    var_gene_overlap = len(set(result.var_names).intersection(discovery_genes))

    if var_cell_overlap >= obs_cell_overlap and obs_gene_overlap >= var_gene_overlap:
        gene_patterns_df = result.obs[obs_pats].copy()
        cell_patterns_df = result.var[var_pats].copy()
        orientation_msg = "Detected genes in result.obs and cells in result.var."
    elif obs_cell_overlap > var_cell_overlap and var_gene_overlap > obs_gene_overlap:
        gene_patterns_df = result.var[var_pats].copy()
        cell_patterns_df = result.obs[obs_pats].copy()
        orientation_msg = "Detected genes in result.var and cells in result.obs."
    else:
        if result.n_obs == discovery_input.n_vars and result.n_vars == discovery_input.n_obs:
            gene_patterns_df = result.obs[obs_pats].copy()
            cell_patterns_df = result.var[var_pats].copy()
            orientation_msg = "Fallback by dimension: genes in obs, cells in var."
        elif result.n_var == discovery_input.n_vars and result.n_obs == discovery_input.n_obs:
            gene_patterns_df = result.var[var_pats].copy()
            cell_patterns_df = result.obs[obs_pats].copy()
            orientation_msg = "Fallback by dimension: genes in var, cells in obs."
        else:
            raise ValueError(
                "Could not confidently determine CoGAPS result orientation. "
                f"obs_cell_overlap={obs_cell_overlap}, var_cell_overlap={var_cell_overlap}, "
                f"obs_gene_overlap={obs_gene_overlap}, var_gene_overlap={var_gene_overlap}"
            )

    pattern_names = sorted(
        list(set(get_pattern_columns(gene_patterns_df.columns.tolist())).intersection(
            set(get_pattern_columns(cell_patterns_df.columns.tolist()))
        )),
        key=numeric_suffix,
    )
    gene_patterns_df = gene_patterns_df[pattern_names].copy()
    cell_patterns_df = cell_patterns_df[pattern_names].copy()

    return gene_patterns_df, cell_patterns_df, orientation_msg

def top_genes_per_pattern(gene_patterns_df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    rows = []
    for pattern in gene_patterns_df.columns:
        top = gene_patterns_df[pattern].sort_values(ascending=False).head(top_n)
        for rank, (gene, weight) in enumerate(top.items(), start=1):
            rows.append({
                "pattern": pattern,
                "rank": rank,
                "gene": gene,
                "weight": float(weight),
            })
    return pd.DataFrame(rows)

# -------------------------
# 11) Q1-Q3 analysis on direct CoGAPS cell scores
# -------------------------
def summarize_cogaps_patterns(cell_meta_patterns: pd.DataFrame, top_gene_df: pd.DataFrame) -> pd.DataFrame:
    pattern_names = get_pattern_columns(cell_meta_patterns.columns.tolist())
    ifn_set = set(PROGRAMS["ifn_program"])

    tp_order = get_timepoint_order(
        cell_meta_patterns[["timepoint_merged", "day_merged_numeric"]].drop_duplicates(),
        "timepoint_merged",
        "day_merged_numeric",
    )

    rows = []
    for p in pattern_names:
        groups = []
        for tp in tp_order:
            vals = cell_meta_patterns.loc[
                cell_meta_patterns["timepoint_merged"].astype(str) == str(tp), p
            ].to_numpy()
            vals = vals[np.isfinite(vals)]
            if len(vals) > 1:
                groups.append(vals)

        if len(groups) >= 2:
            kw_stat, kw_p = kruskal(*groups)
        else:
            kw_stat, kw_p = np.nan, np.nan

        eta_time = eta_squared(
            cell_meta_patterns[p].to_numpy(),
            cell_meta_patterns["timepoint_merged"].astype(str).to_numpy(),
        )
        eta_cell = eta_squared(
            cell_meta_patterns[p].to_numpy(),
            cell_meta_patterns["broad_cell_type"].astype(str).to_numpy(),
        )

        if "ifn_program_score" in cell_meta_patterns.columns:
            rho_ifn, p_ifn = spearmanr(
                cell_meta_patterns[p].to_numpy(),
                cell_meta_patterns["ifn_program_score"].to_numpy(),
                nan_policy="omit",
            )
        else:
            rho_ifn, p_ifn = np.nan, np.nan

        means_by_time = (
            cell_meta_patterns
            .groupby(["timepoint_merged", "day_merged_numeric"], observed=True)[p]
            .mean()
            .reset_index()
            .sort_values("day_merged_numeric")
        )
        peak_tp = means_by_time.loc[means_by_time[p].idxmax(), "timepoint_merged"] if len(means_by_time) else None

        top_genes = top_gene_df.loc[top_gene_df["pattern"] == p, "gene"].head(15).tolist()
        ifn_overlap = len(ifn_set.intersection(top_genes))

        rows.append({
            "pattern": p,
            "kruskal_time_stat": kw_stat,
            "kruskal_time_p": kw_p,
            "eta_timepoint": eta_time,
            "eta_broad_cell_type": eta_cell,
            "pattern_class": classify_pattern(eta_time, eta_cell),
            "spearman_ifn_score_rho": rho_ifn,
            "spearman_ifn_score_p": p_ifn,
            "peak_timepoint": peak_tp,
            "ifn_top_gene_overlap_top15": ifn_overlap,
            "candidate_ifn_pattern": bool(
                (ifn_overlap >= 3) or
                (pd.notnull(rho_ifn) and rho_ifn > 0.35)
            ),
        })

    summary = pd.DataFrame(rows)
    summary["kruskal_time_p_adj"] = bh_adjust(summary["kruskal_time_p"].to_numpy())
    summary["ifn_score_p_adj"] = bh_adjust(summary["spearman_ifn_score_p"].to_numpy())
    summary = summary.sort_values(["candidate_ifn_pattern", "eta_timepoint"], ascending=[False, False])
    return summary

def write_pattern_report(summary_df: pd.DataFrame, top_gene_df: pd.DataFrame, out_path: Path) -> None:
    lines = []
    lines.append("CoGAPS preliminary interpretation for GSE154386 experimental discovery set")
    lines.append("=" * 72)
    lines.append("")

    q1 = summary_df.sort_values("eta_timepoint", ascending=False).head(5)
    q3 = summary_df.loc[summary_df["candidate_ifn_pattern"]].head(5)

    lines.append("Q1. Most time-varying patterns")
    for _, row in q1.iterrows():
        genes = top_gene_df.loc[top_gene_df["pattern"] == row["pattern"], "gene"].head(10).tolist()
        lines.append(
            f"- {row['pattern']}: eta_time={row['eta_timepoint']:.3f}, "
            f"eta_cell={row['eta_broad_cell_type']:.3f}, class={row['pattern_class']}, "
            f"peak={row['peak_timepoint']}, q={row['kruskal_time_p_adj']:.3e}"
        )
        lines.append(f"  Top genes: {', '.join(genes)}")

    lines.append("")
    lines.append("Q3. Strongest interferon-associated candidates")
    if len(q3) == 0:
        lines.append("- No pattern met the current interferon heuristic.")
    else:
        for _, row in q3.iterrows():
            genes = top_gene_df.loc[top_gene_df["pattern"] == row["pattern"], "gene"].head(10).tolist()
            lines.append(
                f"- {row['pattern']}: IFN overlap(top15)={row['ifn_top_gene_overlap_top15']}, "
                f"rho(IFN score)={row['spearman_ifn_score_rho']:.3f}, peak={row['peak_timepoint']}"
            )
            lines.append(f"  Top genes: {', '.join(genes)}")

    write_text(out_path, "\n".join(lines))

# -------------------------
# 12) Projection into natural cohort (Q4)
# -------------------------
def choose_projection_genes(gene_patterns_df: pd.DataFrame, top_genes_per_pattern_n: int = 250) -> List[str]:
    genes = set()
    for pattern in gene_patterns_df.columns:
        s = gene_patterns_df[pattern].replace([np.inf, -np.inf], np.nan).fillna(0.0)
        s = s[s > 0]
        genes.update(
            s.sort_values(ascending=False)
            .head(top_genes_per_pattern_n)
            .index
            .tolist()
        )
    return sorted(genes)

def project_patterns_nnls_sparse(
    target_adata: AnnData,
    gene_patterns_df: pd.DataFrame,
    top_genes_per_pattern_n: int = 250,
    progress_every: int = 500,
) -> pd.DataFrame:
    pattern_names = gene_patterns_df.columns.tolist()
    proj_genes = choose_projection_genes(gene_patterns_df, top_genes_per_pattern_n=top_genes_per_pattern_n)
    common_genes = [g for g in proj_genes if g in target_adata.var_names and g in gene_patterns_df.index]

    if len(common_genes) < MIN_COMMON_PROJECTION_GENES:
        raise ValueError(f"Too few common projection genes found ({len(common_genes)}).")

    A = gene_patterns_df.loc[common_genes, pattern_names].to_numpy(dtype=np.float64)
    A = np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0)
    A[A < 0] = 0.0

    X_sub = target_adata[:, common_genes].X
    if sparse.issparse(X_sub):
        X_sub = X_sub.tocsr()
    else:
        X_sub = np.asarray(X_sub, dtype=np.float64)

    scores = np.zeros((target_adata.n_obs, len(pattern_names)), dtype=np.float32)
    print(
        f"[projection] projecting {target_adata.n_obs} cells with {len(common_genes)} genes onto {len(pattern_names)} patterns",
        flush=True,
    )

    for i in range(target_adata.n_obs):
        if sparse.issparse(X_sub):
            x = X_sub.getrow(i).toarray().ravel().astype(np.float64)
        else:
            x = np.asarray(X_sub[i], dtype=np.float64).ravel()

        x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
        x[x < 0] = 0.0

        coeffs, _ = nnls(A, x)
        scores[i, :] = coeffs.astype(np.float32)

        if (i + 1) % progress_every == 0:
            print(f"[projection] processed {i + 1}/{target_adata.n_obs} cells", flush=True)

    return pd.DataFrame(scores, index=target_adata.obs_names, columns=pattern_names)

# -------------------------
# 13) Main
# -------------------------
def main() -> None:
    np.random.seed(RANDOM_SEED)
    sc.settings.figdir = str(FIG_DIR)
    sc.settings.set_figure_params(dpi=120, figsize=(6, 5))

    # Build or load the preprocessed HVG object
    if USE_CACHED_HVG and OUT_ALL_HVG.exists():
        print(f"[cache] loading existing HVG object: {OUT_ALL_HVG}", flush=True)
        adata_hvg = sc.read_h5ad(OUT_ALL_HVG)
        add_merged_timepoints(adata_hvg)
    else:
        adata_hvg = build_preprocessed_hvg()

    # Split experimental / natural
    adata_exp_all = adata_hvg[adata_hvg.obs["cohort"] == "experimental"].copy()
    adata_nat = adata_hvg[adata_hvg.obs["cohort"] == "natural"].copy()

    # Build balanced experimental discovery set
    adata_discovery = balanced_sample_by_group(
        adata_exp_all,
        groupby=DISCOVERY_SAMPLE_GROUP,
        max_cells_per_group=DISCOVERY_MAX_CELLS_PER_SAMPLE,
        seed=RANDOM_SEED,
    )
    add_merged_timepoints(adata_discovery)

    discovery_counts = (
        adata_discovery.obs
        .groupby(["sample_id", "subject", "timepoint_merged", "day_merged_numeric", "broad_cell_type"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["subject", "day_merged_numeric", "broad_cell_type"])
    )
    discovery_counts.to_csv(RESULTS_DIR / "experimental_discovery_cell_counts_nondist.csv", index=False)

    # Keep natural target exactly as before
    adata_nat_target = make_cogaps_ready_sparse(adata_nat)

    # Pretranspose discovery data into dense genes x cells for PyCoGAPS
    adata_discovery_cogaps = make_cogaps_ready_pretransposed_dense(adata_discovery)

    check_matrix("discovery_cogaps_X", adata_discovery_cogaps.X)
    check_matrix("natural_target_X", adata_nat_target.X)

    adata_discovery_cogaps.write_h5ad(OUT_DISCOVERY_COGAPS)
    adata_nat_target.write_h5ad(OUT_NAT_TARGET)

    print(f"[saved discovery CoGAPS input] {OUT_DISCOVERY_COGAPS}", flush=True)
    print(f"[saved natural projection target] {OUT_NAT_TARGET}", flush=True)

    # Run NON-DISTRIBUTED CoGAPS
    cogaps_result = run_cogaps_nondistributed(adata_discovery_cogaps, OUT_COGAPS_RESULT)

    # pass ORIGINAL discovery object (cells x genes), not the pretransposed CoGAPS input
    gene_patterns_df, cell_patterns_df, orientation_msg = extract_cogaps_matrices(
        cogaps_result,
        adata_discovery,
    )
    print(f"[CoGAPS] {orientation_msg}", flush=True)
    write_text(COGAPS_DIR / "cogaps_orientation.txt", orientation_msg)

    # FIX: sanitize CoGAPS outputs before downstream use
    gene_patterns_df = clean_pattern_df(gene_patterns_df, clip_negative=True)
    cell_patterns_df = clean_pattern_df(cell_patterns_df, clip_negative=True)

    pattern_names = gene_patterns_df.columns.tolist()

    gene_patterns_df.to_csv(COGAPS_DIR / "cogaps_gene_loadings.csv")
    cell_patterns_df.to_csv(COGAPS_DIR / "cogaps_cell_scores.csv")

    top_gene_df = top_genes_per_pattern(gene_patterns_df, top_n=50)
    top_gene_df.to_csv(COGAPS_DIR / "cogaps_pattern_top_genes.csv", index=False)

    # Join CoGAPS cell scores to discovery metadata
    cell_meta_patterns = adata_discovery.obs.join(cell_patterns_df, how="inner")
    cell_meta_patterns.to_csv(COGAPS_DIR / "discovery_cells_with_cogaps_patterns.csv")

    # Visualize pattern scores on discovery UMAP
    adata_discovery_patterns = adata_discovery.copy()
    for p in pattern_names:
        adata_discovery_patterns.obs[p] = cell_patterns_df.loc[adata_discovery_patterns.obs_names, p]
    save_umap_panels(
        adata_discovery_patterns,
        colors=pattern_names,
        out_prefix="discovery_cogaps_umap_nondist",
    )

    # Q1-Q3 summaries from direct CoGAPS cell scores
    pattern_summary_df = summarize_cogaps_patterns(cell_meta_patterns, top_gene_df)
    pattern_summary_df.to_csv(COGAPS_DIR / "cogaps_pattern_summary.csv", index=False)

    pattern_summary_df.sort_values("eta_timepoint", ascending=False).to_csv(
        COGAPS_DIR / "RQ1_time_varying_patterns.csv",
        index=False,
    )
    pattern_summary_df[[
        "pattern", "eta_timepoint", "eta_broad_cell_type", "pattern_class", "peak_timepoint"
    ]].sort_values(["pattern_class", "eta_timepoint"], ascending=[True, False]).to_csv(
        COGAPS_DIR / "RQ2_identity_vs_activity_patterns.csv",
        index=False,
    )
    pattern_summary_df.loc[pattern_summary_df["candidate_ifn_pattern"]].to_csv(
        COGAPS_DIR / "RQ3_interferon_candidate_patterns.csv",
        index=False,
    )

    # Experimental pattern means over time / cell type
    pattern_time_means = (
        cell_meta_patterns
        .groupby(["timepoint_merged", "day_merged_numeric"], observed=True)[pattern_names]
        .mean()
        .reset_index()
        .sort_values("day_merged_numeric")
    )
    pattern_time_means.to_csv(COGAPS_DIR / "experimental_pattern_means_by_time.csv", index=False)

    pattern_time_cell_means = (
        cell_meta_patterns
        .groupby(["broad_cell_type", "timepoint_merged", "day_merged_numeric"], observed=True)[pattern_names]
        .mean()
        .reset_index()
        .sort_values(["broad_cell_type", "day_merged_numeric"])
    )
    pattern_time_cell_means.to_csv(COGAPS_DIR / "experimental_pattern_means_by_time_and_celltype.csv", index=False)

    exp_long = pattern_time_means.melt(
        id_vars=["timepoint_merged", "day_merged_numeric"],
        value_vars=pattern_names,
        var_name="pattern",
        value_name="mean_usage",
    )
    save_heatmap(
        df=exp_long,
        index_col="timepoint_merged",
        column_col="pattern",
        value_col="mean_usage",
        day_col="day_merged_numeric",
        out_png=FIG_DIR / "experimental_cogaps_pattern_heatmap_nondist.png",
        title="Experimental CoGAPS pattern usage across time (non-distributed)",
    )

    top_patterns_for_lines = pattern_summary_df.sort_values("eta_timepoint", ascending=False)["pattern"].head(4).tolist()
    for p in top_patterns_for_lines:
        tmp = (
            cell_meta_patterns
            .groupby(["broad_cell_type", "timepoint_merged", "day_merged_numeric"], observed=True)[p]
            .mean()
            .reset_index()
            .sort_values("day_merged_numeric")
        )
        save_pattern_lineplot(
            summary_df=tmp,
            pattern_col=p,
            out_png=FIG_DIR / f"{p}_experimental_by_celltype_nondist.png",
            title=f"{p} over experimental time course (non-distributed)",
        )

    write_pattern_report(
        summary_df=pattern_summary_df,
        top_gene_df=top_gene_df,
        out_path=COGAPS_DIR / "cogaps_preliminary_report.txt",
    )

    # Q4: project learned gene weights into natural cohort
    natural_projected_df = project_patterns_nnls_sparse(
        target_adata=adata_nat_target,
        gene_patterns_df=gene_patterns_df,
        top_genes_per_pattern_n=PROJECTION_TOP_GENES_PER_PATTERN,
        progress_every=PROJECTION_PROGRESS_EVERY,
    )
    natural_meta_patterns = adata_nat.obs.join(natural_projected_df, how="inner")
    natural_meta_patterns.to_csv(COGAPS_DIR / "natural_cells_projected_patterns.csv")

    natural_time_means = (
        natural_meta_patterns
        .groupby(["timepoint_merged", "day_merged_numeric"], observed=True)[pattern_names]
        .mean()
        .reset_index()
        .sort_values("day_merged_numeric")
    )
    natural_time_means.to_csv(COGAPS_DIR / "RQ4_natural_projection_summary.csv", index=False)

    nat_long = natural_time_means.melt(
        id_vars=["timepoint_merged", "day_merged_numeric"],
        value_vars=pattern_names,
        var_name="pattern",
        value_name="mean_usage",
    )
    save_heatmap(
        df=nat_long,
        index_col="timepoint_merged",
        column_col="pattern",
        value_col="mean_usage",
        day_col="day_merged_numeric",
        out_png=FIG_DIR / "natural_projection_pattern_heatmap_nondist.png",
        title="Projection of experimental CoGAPS patterns into natural infection (non-distributed)",
    )

    natural_time_cell_means = (
        natural_meta_patterns
        .groupby(["broad_cell_type", "timepoint_merged", "day_merged_numeric"], observed=True)[pattern_names]
        .mean()
        .reset_index()
        .sort_values(["broad_cell_type", "day_merged_numeric"])
    )
    natural_time_cell_means.to_csv(COGAPS_DIR / "natural_projection_by_time_and_celltype.csv", index=False)

    summary_lines = [
        "GSE154386 NON-DISTRIBUTED CoGAPS analysis summary",
        "=" * 60,
        f"All HVG dataset: {adata_hvg.shape}",
        f"Experimental all-cells dataset: {adata_exp_all.shape}",
        f"Experimental discovery CoGAPS dataset: {adata_discovery_cogaps.shape}",
        f"Natural projection target dataset: {adata_nat_target.shape}",
        "",
        "Key outputs:",
        f"- {OUT_ALL_HVG}",
        f"- {OUT_DISCOVERY_COGAPS}",
        f"- {OUT_COGAPS_RESULT}",
        f"- {OUT_NAT_TARGET}",
        f"- {COGAPS_DIR / 'RQ1_time_varying_patterns.csv'}",
        f"- {COGAPS_DIR / 'RQ2_identity_vs_activity_patterns.csv'}",
        f"- {COGAPS_DIR / 'RQ3_interferon_candidate_patterns.csv'}",
        f"- {COGAPS_DIR / 'RQ4_natural_projection_summary.csv'}",
    ]
    write_text(RESULTS_DIR / "run_summary_nondist.txt", "\n".join(summary_lines))

    print("\n✅ Summary:", flush=True)
    print(" - All HVG dataset:", adata_hvg.shape, flush=True)
    print(" - Experimental discovery CoGAPS input:", adata_discovery_cogaps.shape, flush=True)
    print(" - CoGAPS result saved to:", OUT_COGAPS_RESULT, flush=True)
    print(" - Natural projection target:", adata_nat_target.shape, flush=True)
    print(" - Research-question outputs written to:", COGAPS_DIR, flush=True)
    print(" - Figures written to:", FIG_DIR, flush=True)

if __name__ == "__main__":
    main()