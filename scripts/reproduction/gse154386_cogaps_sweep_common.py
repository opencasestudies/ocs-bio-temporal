#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


REPO_ROOT = Path(__file__).resolve().parent
DATA_DIR = REPO_ROOT / "GSE154386"
DEFAULT_OUTDIR = DATA_DIR / "cogaps_sweep_singleprocess_hpc"

CASE_STUDY5_SCRIPT_CANDIDATES = [
    "gse154386_sparse_nondistributed_cogaps.py",
    "gse154386_sparse_distributed_cogaps.py",
]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def parse_int_csv(text: str) -> List[int]:
    values = [chunk.strip() for chunk in str(text).split(",")]
    out = [int(chunk) for chunk in values if chunk]
    if not out:
        raise ValueError(f"Expected at least one integer in CSV string: {text!r}")
    return out


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_path(path: Path | str) -> Path:
    raw = Path(path)
    return raw if raw.is_absolute() else (REPO_ROOT / raw).resolve()


def find_case_study5_script() -> Path:
    checked: List[str] = []
    for name in CASE_STUDY5_SCRIPT_CANDIDATES:
        candidate = REPO_ROOT / name
        checked.append(str(candidate))
        if candidate.exists():
            return candidate.resolve()
    checked_str = "\n".join(f"  - {entry}" for entry in checked)
    raise FileNotFoundError(
        "Could not find the Case Study 5 source script. Checked:\n"
        f"{checked_str}"
    )


def load_reference_module(script_path: Optional[Path | str] = None):
    source_path = resolve_path(script_path) if script_path else find_case_study5_script()
    module_name = f"gse154386_case_study5_source_{abs(hash(str(source_path)))}"
    spec = importlib.util.spec_from_file_location(module_name, source_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create import spec for {source_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def json_ready(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_ready(v) for v in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(json_ready(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def metrics_status_ok(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        payload = read_json(path)
    except Exception:
        return False
    return str(payload.get("status", "")).lower() == "ok"


def run_stem(k: int, seed: int, n_iter: int) -> str:
    return f"cogaps_K{k}_seed{seed}_iter{n_iter}"


def run_artifact_paths(outdir: Path | str, k: int, seed: int, n_iter: int) -> Dict[str, Path]:
    root = resolve_path(outdir)
    runs_dir = ensure_dir(root / "runs")
    logs_dir = ensure_dir(root / "logs")
    figures_dir = ensure_dir(root / "figures")
    stem = run_stem(k=k, seed=seed, n_iter=n_iter)
    return {
        "outdir": root,
        "runs_dir": runs_dir,
        "logs_dir": logs_dir,
        "figures_dir": figures_dir,
        "result_path": runs_dir / f"{stem}.h5ad",
        "metrics_path": runs_dir / f"{stem}.metrics.json",
        "gene_loadings_csv": runs_dir / f"{stem}.gene_loadings.csv",
        "cell_scores_csv": runs_dir / f"{stem}.cell_scores.csv",
        "top_genes_csv": runs_dir / f"{stem}.top_genes.csv",
        "pattern_summary_csv": runs_dir / f"{stem}.pattern_summary.csv",
        "cell_metadata_patterns_csv": runs_dir / f"{stem}.discovery_cells_with_patterns.csv",
        "log_path": logs_dir / f"run_{stem}.log",
    }


def cache_paths(outdir: Path | str) -> Dict[str, Path]:
    root = resolve_path(outdir)
    cache_dir = ensure_dir(root / "cache")
    return {
        "outdir": root,
        "cache_dir": cache_dir,
        "preprocessed_hvg_h5ad": cache_dir / "gse154386_preprocessed_hvg.h5ad",
        "experimental_discovery_h5ad": cache_dir / "gse154386_experimental_discovery_cells_x_genes.h5ad",
        "cogaps_input_h5ad": cache_dir / "gse154386_experimental_discovery_genes_x_cells.h5ad",
        "natural_target_h5ad": cache_dir / "gse154386_natural_projection_target_cells_x_genes.h5ad",
        "discovery_counts_csv": cache_dir / "gse154386_experimental_discovery_cell_counts.csv",
        "prep_manifest_json": cache_dir / "gse154386_prep_manifest.json",
    }


def set_thread_env(n_threads: Optional[int]) -> None:
    if n_threads is None:
        return
    value = str(int(n_threads))
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[key] = value


def summarize_missing_paths(paths: Iterable[Path]) -> List[str]:
    return [str(path) for path in paths if not Path(path).exists()]
