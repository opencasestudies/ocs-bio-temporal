#!/usr/bin/env python3
"""Generate a lightweight pattern-score embedding figure for Case Study 5.

The preferred layout is the transcriptome UMAP stored in the full preprocessed
HVG AnnData object. Only the discovery-cell coordinates and selected-model
pattern scores are exported into the repository. If the large H5AD is not
available, the script falls back to a UMAP of the included CoGAPS pattern-score
matrix so the figure can still be regenerated locally.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "data" / "processed" / "figures"
SOURCE_DIR = FIG_DIR / "source_tables"
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

DISCOVERY_CSV = (
    ROOT
    / "data"
    / "processed"
    / "selected_model_k10"
    / "r"
    / "cogaps_K10_seed2_iter2000.discovery_cells_with_patterns.csv"
)

DEFAULT_FULL_HVG = ROOT / "data" / "external" / "gse154386_preprocessed_hvg.h5ad"
FULL_HVG = Path(os.environ.get("CS5_FULL_HVG_H5AD", DEFAULT_FULL_HVG))

SOURCE_CSV = SOURCE_DIR / "figure_08_transcriptome_umap_pattern_scores_source.csv"
OUT_PNG = FIG_DIR / "figure_08_pattern_scores_on_transcriptome_umap.png"
OUT_SVG = FIG_DIR / "figure_08_pattern_scores_on_transcriptome_umap.svg"
MANIFEST = FIG_DIR / "figure_manifest.json"

PATTERN_COLUMNS = [f"Pattern{i}" for i in range(1, 11)]
PLOT_PATTERNS = [
    ("Pattern3", "Pattern3\nIFN-associated activity"),
    ("Pattern7", "Pattern7\nT-cell identity"),
    ("Pattern10", "Pattern10\nB-cell identity"),
]
TIMEPOINT_ORDER = ["D0", "D2", "D4", "D6", "D8", "D10", "D14/15", "D28"]
CELL_TYPE_ORDER = ["T_cell", "NK_cell", "Monocyte", "B_cell", "Neutrophil", "Plasmablast"]
CELL_TYPE_COLORS = {
    "T_cell": "#66c2a5",
    "NK_cell": "#fc8d62",
    "Monocyte": "#8da0cb",
    "B_cell": "#e78ac3",
    "Neutrophil": "#a6d854",
    "Plasmablast": "#ffd92f",
}


def _read_discovery_scores() -> pd.DataFrame:
    discovery = pd.read_csv(DISCOVERY_CSV, index_col=0)
    discovery.index.name = "cell_id"
    keep_cols = [
        "subject",
        "timepoint_merged",
        "day_merged_numeric",
        "broad_cell_type",
        *PATTERN_COLUMNS,
    ]
    missing = [col for col in keep_cols if col not in discovery.columns]
    if missing:
        raise ValueError(f"Missing expected columns in {DISCOVERY_CSV}: {missing}")
    out = discovery[keep_cols].copy()
    out["cell_id"] = out.index.astype(str)
    return out.reset_index(drop=True)


def _transcriptome_umap(discovery: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    if not FULL_HVG.exists():
        raise FileNotFoundError(FULL_HVG)

    import anndata as ad

    adata = ad.read_h5ad(FULL_HVG, backed="r")
    try:
        if "X_umap" not in adata.obsm:
            raise KeyError("X_umap")
        coords = pd.DataFrame(
            np.asarray(adata.obsm["X_umap"]),
            index=adata.obs_names.astype(str),
            columns=["umap_1", "umap_2"],
        )
    finally:
        adata.file.close()

    merged = discovery.merge(
        coords,
        how="left",
        left_on="cell_id",
        right_index=True,
        validate="one_to_one",
    )
    missing = int(merged["umap_1"].isna().sum())
    if missing:
        raise ValueError(f"{missing} discovery cells were missing transcriptome UMAP coordinates")
    merged["embedding_source"] = "transcriptome_umap_from_full_preprocessed_hvg"
    return merged, "Transcriptome UMAP"


def _pattern_score_umap(discovery: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    try:
        import umap
    except ImportError as exc:
        raise RuntimeError(
            "The full H5AD was unavailable and umap-learn is not installed for fallback."
        ) from exc

    matrix = discovery[PATTERN_COLUMNS].to_numpy(dtype=float)
    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.35,
        metric="euclidean",
        random_state=42,
    )
    coords = reducer.fit_transform(matrix)
    out = discovery.copy()
    out["umap_1"] = coords[:, 0]
    out["umap_2"] = coords[:, 1]
    out["embedding_source"] = "fallback_umap_from_k10_pattern_scores"
    return out, "Pattern-score UMAP"


def make_source_table() -> tuple[pd.DataFrame, str]:
    discovery = _read_discovery_scores()
    try:
        table, label = _transcriptome_umap(discovery)
    except Exception as exc:
        print(f"[warn] using pattern-score UMAP fallback: {exc}")
        table, label = _pattern_score_umap(discovery)

    ordered_cols = [
        "cell_id",
        "umap_1",
        "umap_2",
        "embedding_source",
        "subject",
        "timepoint_merged",
        "day_merged_numeric",
        "broad_cell_type",
        *PATTERN_COLUMNS,
    ]
    table = table[ordered_cols].copy()
    table["timepoint_merged"] = pd.Categorical(
        table["timepoint_merged"],
        categories=TIMEPOINT_ORDER,
        ordered=True,
    )
    table["broad_cell_type"] = pd.Categorical(
        table["broad_cell_type"],
        categories=CELL_TYPE_ORDER,
        ordered=True,
    )
    table.to_csv(SOURCE_CSV, index=False)
    print(f"[write] {SOURCE_CSV} rows={len(table)} source={label}")
    return table, label


def _scatter_base(ax: plt.Axes, table: pd.DataFrame) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1", fontsize=9)
    ax.set_ylabel("UMAP 2", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)


def make_figure(table: pd.DataFrame, embedding_label: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.0), constrained_layout=True)
    fig.suptitle(
        "Selected CoGAPS pattern scores on a transcriptome embedding",
        fontsize=18,
        fontweight="bold",
        y=1.02,
    )

    ax0 = axes[0, 0]
    for cell_type in CELL_TYPE_ORDER:
        sub = table.loc[table["broad_cell_type"].astype(str).eq(cell_type)]
        ax0.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=5,
            alpha=0.72,
            linewidths=0,
            color=CELL_TYPE_COLORS[cell_type],
            label=cell_type,
        )
    ax0.set_title("Broad PBMC identity", loc="left", fontsize=12.5, fontweight="bold")
    _scatter_base(ax0, table)
    handles = [
        Line2D([0], [0], marker="o", color="w", label=ct, markerfacecolor=CELL_TYPE_COLORS[ct], markersize=7)
        for ct in CELL_TYPE_ORDER
    ]
    ax0.legend(
        handles=handles,
        frameon=False,
        fontsize=8,
        loc="upper right",
        bbox_to_anchor=(1.04, 1.02),
        borderaxespad=0,
    )

    for ax, (pattern, title) in zip(axes.ravel()[1:], PLOT_PATTERNS):
        values = table[pattern].astype(float)
        vmax = float(np.nanquantile(values, 0.99))
        scatter = ax.scatter(
            table["umap_1"],
            table["umap_2"],
            c=values,
            s=5,
            alpha=0.82,
            linewidths=0,
            cmap="viridis",
            vmin=0,
            vmax=vmax if vmax > 0 else None,
        )
        ax.set_title(title, loc="left", fontsize=12.5, fontweight="bold")
        _scatter_base(ax, table)
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label("Cell score", fontsize=9)

    fig.text(
        0.01,
        -0.01,
        (
            f"Coordinates: {embedding_label}. Colors in pattern panels show K=10 CoGAPS cell scores. "
            "The embedding is a context view, not a statistical test."
        ),
        ha="left",
        va="bottom",
        fontsize=10,
        color="#486581",
    )
    fig.savefig(OUT_PNG, dpi=220, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {OUT_PNG}")
    print(f"[write] {OUT_SVG}")


def update_manifest(embedding_label: str) -> None:
    if MANIFEST.exists():
        manifest = json.loads(MANIFEST.read_text())
    else:
        manifest = {}
    figures = manifest.setdefault("figures", {})
    figures["figure_08_pattern_scores_on_transcriptome_umap_png"] = str(OUT_PNG.relative_to(ROOT))
    figures["figure_08_pattern_scores_on_transcriptome_umap_svg"] = str(OUT_SVG.relative_to(ROOT))
    figures["figure_08_pattern_scores_on_transcriptome_umap_source_csv"] = str(SOURCE_CSV.relative_to(ROOT))
    manifest["figure_08_embedding_source"] = embedding_label
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[write] {MANIFEST}")


def main() -> None:
    table, embedding_label = make_source_table()
    make_figure(table, embedding_label)
    update_manifest(embedding_label)


if __name__ == "__main__":
    main()
