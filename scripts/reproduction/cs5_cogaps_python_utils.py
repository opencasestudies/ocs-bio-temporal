#!/usr/bin/env python3
from __future__ import annotations

import json
import logging
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def setup_logger(log_path: Path, name: str = "cs5_cogaps") -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.INFO)
    logger.addHandler(file_handler)
    return logger


def set_thread_env(*, cogaps_threads: int, blas_threads: int) -> None:
    os.environ["OMP_NUM_THREADS"] = str(int(cogaps_threads))
    os.environ["OPENBLAS_NUM_THREADS"] = str(int(blas_threads))
    os.environ["MKL_NUM_THREADS"] = str(int(blas_threads))
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(int(blas_threads))
    os.environ["NUMEXPR_NUM_THREADS"] = str(int(blas_threads))


def json_ready(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, (np.generic,)):
        return json_ready(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    if pd.isna(value) if not isinstance(value, (str, bytes, list, tuple, dict, set)) else False:
        return None
    return value


def atomic_write_json(payload: Dict[str, Any], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True), encoding="utf-8")
    tmp_path.replace(out_path)


def atomic_write_h5ad(adata: ad.AnnData, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    adata.write_h5ad(tmp_path)
    if tmp_path.stat().st_size == 0:
        raise RuntimeError(f"Temporary h5ad write is empty: {tmp_path}")
    tmp_path.replace(out_path)


def pattern_columns(columns: Iterable[str]) -> List[str]:
    out = []
    for column in columns:
        text = str(column)
        if text.startswith("Pattern_"):
            suffix = text.split("_", 1)[1]
        elif text.startswith("Pattern"):
            suffix = text.replace("Pattern", "", 1)
        else:
            continue
        if suffix.isdigit():
            out.append(text)
    return sorted(out, key=lambda name: int(str(name).replace("Pattern_", "").replace("Pattern", "")))


def safe_shape(value: Any) -> Optional[List[int]]:
    if value is None:
        return None
    try:
        return list(np.asarray(value).shape)
    except Exception:
        return None


def summarize_result_uns(result: ad.AnnData) -> Dict[str, Any]:
    uns = result.uns
    summary: Dict[str, Any] = {}
    scalar_keys = (
        "meanChiSq",
        "totalRunningTime",
        "totalUpdates",
        "seed",
        "averageQueueLengthA",
        "averageQueueLengthP",
    )
    for key in scalar_keys:
        if key in uns:
            summary[key] = json_ready(uns[key])

    for key in ("chisqHistory", "atomhistoryA", "atomhistoryP"):
        if key in uns:
            try:
                summary[f"{key}_length"] = int(len(uns[key]))
            except TypeError:
                summary[f"{key}_length"] = None

    for key in (
        "equilibrationSnapshotsA",
        "equilibrationSnapshotsP",
        "samplingSnapshotsA",
        "samplingSnapshotsP",
        "meanPatternAssignment",
        "pumpMatrix",
        "pumpStat",
    ):
        if key in uns:
            summary[f"{key}_shape"] = safe_shape(uns[key])
            try:
                summary[f"{key}_length"] = int(len(uns[key]))
            except TypeError:
                pass
    return summary


def trace_dataframe(result: ad.AnnData) -> pd.DataFrame:
    uns = result.uns
    chisq = list(uns.get("chisqHistory", []))
    atoms_a = list(uns.get("atomhistoryA", []))
    atoms_p = list(uns.get("atomhistoryP", []))
    n_rows = max(len(chisq), len(atoms_a), len(atoms_p))
    rows = []
    for idx in range(n_rows):
        rows.append(
            {
                "trace_index": idx + 1,
                "chisq": chisq[idx] if idx < len(chisq) else np.nan,
                "atomsA": atoms_a[idx] if idx < len(atoms_a) else np.nan,
                "atomsP": atoms_p[idx] if idx < len(atoms_p) else np.nan,
            }
        )
    return pd.DataFrame(rows, columns=["trace_index", "chisq", "atomsA", "atomsP"])


def effective_pattern_summary(gene_patterns_df: pd.DataFrame, cell_patterns_df: pd.DataFrame) -> Tuple[List[Dict[str, Any]], int, List[str]]:
    rows: List[Dict[str, Any]] = []
    for pattern in gene_patterns_df.columns:
        gene_vec = (
            gene_patterns_df[pattern]
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0.0)
            .to_numpy(dtype=float)
        )
        cell_vec = (
            cell_patterns_df[pattern]
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0.0)
            .to_numpy(dtype=float)
        )
        gene_max = float(np.max(gene_vec)) if gene_vec.size else 0.0
        cell_max = float(np.max(cell_vec)) if cell_vec.size else 0.0
        nonzero_gene_count = int(np.count_nonzero(gene_vec > 0))
        effective = bool(gene_max > 0 and cell_max > 0 and float(np.std(cell_vec)) > 1e-8 and nonzero_gene_count >= 10)
        rows.append(
            {
                "pattern": str(pattern),
                "effective": effective,
                "gene_loading_max": gene_max,
                "cell_score_max": cell_max,
                "nonzero_gene_count": nonzero_gene_count,
            }
        )
    n_effective = sum(1 for row in rows if row["effective"])
    degenerate = [row["pattern"] for row in rows if not row["effective"]]
    return rows, n_effective, degenerate


def redundancy_summary(gene_patterns_df: pd.DataFrame, top_gene_df: pd.DataFrame) -> Dict[str, float]:
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
    jaccards = []
    pattern_names = list(gene_patterns_df.columns)
    for i, left in enumerate(pattern_names):
        for right in pattern_names[i + 1 :]:
            left_set = top_gene_sets.get(str(left), set())
            right_set = top_gene_sets.get(str(right), set())
            union = left_set | right_set
            jaccards.append(float(len(left_set & right_set) / len(union)) if union else 0.0)

    return {
        "within_run_pattern_redundancy_mean": float(np.mean(tri)) if tri.size else 0.0,
        "within_run_pattern_redundancy_max": float(np.max(tri)) if tri.size else 0.0,
        "within_run_top_gene_jaccard_mean": float(np.mean(jaccards)) if jaccards else 0.0,
        "within_run_top_gene_jaccard_max": float(np.max(jaccards)) if jaccards else 0.0,
    }


def top_genes_by_pattern(top_gene_df: pd.DataFrame) -> Dict[str, List[str]]:
    return {
        str(pattern): group["gene"].astype(str).tolist()
        for pattern, group in top_gene_df.groupby("pattern", observed=True)
    }
