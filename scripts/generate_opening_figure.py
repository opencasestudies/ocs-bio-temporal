#!/usr/bin/env python3
"""Generate the opening roadmap figure for Case Study 5.

The opening figure is a teaching map for the whole case study. It is generated
from lightweight, GitHub-safe source tables that are also used elsewhere in the
rendered case study.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "data" / "processed" / "figures"
SOURCE_DIR = FIG_DIR / "source_tables"

OUT_PNG = FIG_DIR / "figure_00_opening_summary.png"
OUT_SVG = FIG_DIR / "figure_00_opening_summary.svg"

TIMEPOINTS = ["D0", "D2", "D4", "D6", "D8", "D10", "D14/15", "D28"]
CELL_TYPE_ORDER = ["T_cell", "NK_cell", "Monocyte", "B_cell", "Neutrophil", "Plasmablast"]
CELL_TYPE_COLORS = {
    "T_cell": "#66c2a5",
    "NK_cell": "#fc8d62",
    "Monocyte": "#8da0cb",
    "B_cell": "#e78ac3",
    "Neutrophil": "#a6d854",
    "Plasmablast": "#ffd92f",
}


def fs(size: float) -> float:
    """Scale figure text so it remains legible when rendered inline."""
    return size * 1.25


def _order_timepoints(labels: list[str]) -> list[str]:
    return [x for x in TIMEPOINTS if x in labels]


def _clean_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)


def draw_dataset_panel(ax: plt.Axes) -> None:
    ax.set_title(
        "A. Discovery time course",
        loc="left",
        fontsize=fs(12.2),
        fontweight="bold",
        pad=8,
    )
    positions = list(range(len(TIMEPOINTS)))
    ax.set_xlim(-0.4, len(TIMEPOINTS) - 0.6)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.hlines(0.48, positions[0], positions[-1], color="#53606c", linewidth=2.3, zorder=1)
    late_start, late_end = 4.65, 6.35
    ax.add_patch(
        Rectangle(
            (late_start, 0.25),
            late_end - late_start,
            0.46,
            facecolor="#fee2e2",
            edgecolor="none",
            alpha=0.75,
            zorder=0,
        )
    )
    ax.text(
        (late_start + late_end) / 2,
        0.78,
        "late acute\nwindow",
        ha="center",
        va="bottom",
        fontsize=fs(9.5),
        color="#9f2f20",
        fontweight="bold",
    )

    for label, x_pos in zip(TIMEPOINTS, positions):
        is_late = label in {"D10", "D14/15"}
        face = "#c23b22" if is_late else "#eaf4ff"
        edge = "#8f2a1a" if is_late else "#233f5f"
        text_color = "white" if is_late else "#12263a"
        ax.scatter(x_pos, 0.48, s=540, facecolor=face, edgecolor=edge, linewidth=1.9, zorder=3)
        ax.text(
            x_pos,
            0.48,
            label.replace("/", "/\n") if label == "D14/15" else label,
            ha="center",
            va="center",
            fontsize=fs(9.5),
            fontweight="bold",
            color=text_color,
            zorder=4,
        )

    ax.text(
        positions[0],
        0.96,
        "Three dengue-naive adults, baseline to recovery",
        fontsize=fs(10.8),
        color="#243b53",
        va="top",
    )
    ax.text(
        positions[0],
        0.06,
        "Teaching subset: 600 cells/day; broad PBMC identities retained",
        fontsize=fs(10.2),
        color="#486581",
        va="bottom",
    )


def draw_composition_panel(ax: plt.Axes, composition: pd.DataFrame) -> None:
    ax.set_title(
        "B. Discovery composition",
        loc="left",
        fontsize=fs(12.5),
        fontweight="bold",
        pad=8,
    )

    comp = composition.copy()
    comp["timepoint_merged"] = pd.Categorical(
        comp["timepoint_merged"],
        categories=TIMEPOINTS,
        ordered=True,
    )
    comp["broad_cell_type"] = pd.Categorical(
        comp["broad_cell_type"],
        categories=CELL_TYPE_ORDER,
        ordered=True,
    )
    pivot = (
        comp.pivot_table(
            index="timepoint_merged",
            columns="broad_cell_type",
            values="n_cells",
            aggfunc="sum",
            fill_value=0,
            observed=False,
        )
        .reindex(TIMEPOINTS)
        .loc[:, CELL_TYPE_ORDER]
    )

    x = range(len(pivot.index))
    bottom = [0] * len(pivot.index)
    for cell_type in CELL_TYPE_ORDER:
        values = pivot[cell_type].to_numpy()
        ax.bar(
            x,
            values,
            bottom=bottom,
            width=0.72,
            color=CELL_TYPE_COLORS[cell_type],
            edgecolor="white",
            linewidth=0.8,
            label=cell_type,
        )
        bottom = [b + v for b, v in zip(bottom, values)]

    ax.set_xticks(list(x))
    ax.set_xticklabels(TIMEPOINTS, fontsize=fs(8.7))
    ax.set_ylabel("Discovery cells", fontsize=fs(9.5))
    ax.set_ylim(0, max(bottom) * 1.25)
    ax.grid(axis="y", color="#e6edf3", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.text(
        0.01,
        0.96,
        "About 600 cells per timepoint; colors show retained broad immune-cell identities.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs(8.7),
        color="#486581",
    )
    ax.legend(
        ncol=3,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.995, 0.98),
        fontsize=fs(7.6),
        columnspacing=1.0,
        handlelength=1.25,
    )
    _clean_axis(ax)


def draw_k_temporal_resolution_panel(ax: plt.Axes, k_summary: pd.DataFrame) -> None:
    ax.set_title(
        "C1. K=10 is stable and temporal",
        loc="left",
        fontsize=fs(12.0),
        fontweight="bold",
        pad=8,
    )
    k_summary = k_summary.sort_values("K")

    colors = []
    for _, row in k_summary.iterrows():
        if bool(row["selected_for_case_study"]):
            colors.append("#c44e52")
        elif bool(row["on_stability_plateau"]) and row["max_eta_timepoint"] < 0.12:
            colors.append("#4c72b0")
        else:
            colors.append("#b8b8b8")

    ax.scatter(
        k_summary["stability_core"],
        k_summary["max_eta_timepoint"],
        s=62,
        c=colors,
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
    )

    selected = k_summary.loc[k_summary["selected_for_case_study"].astype(bool)]
    compact = k_summary.loc[k_summary["K"].eq(5)]
    if len(selected) == 1:
        row = selected.iloc[0]
        ax.annotate(
            "K=10:\nstable + temporal",
            xy=(row["stability_core"], row["max_eta_timepoint"]),
            xytext=(row["stability_core"] - 0.032, row["max_eta_timepoint"] - 0.13),
            arrowprops=dict(arrowstyle="->", color="#9f2f20", lw=1.0),
            fontsize=fs(8.0),
            color="#9f2f20",
            fontweight="bold",
        )
    if len(compact) == 1:
        row = compact.iloc[0]
        ax.annotate(
            "K=5:\nstable, compact",
            xy=(row["stability_core"], row["max_eta_timepoint"]),
            xytext=(row["stability_core"] - 0.032, row["max_eta_timepoint"] + 0.17),
            arrowprops=dict(arrowstyle="->", color="#4c72b0", lw=1.0),
            fontsize=fs(8.0),
            color="#2f5f9f",
        )

    for _, row in k_summary.iterrows():
        ax.text(
            row["stability_core"] + 0.002,
            row["max_eta_timepoint"] + 0.004,
            f"K={int(row['K'])}",
            fontsize=fs(7.4),
            color="#263238",
        )

    ax.set_xlabel("Seed-to-seed stability core", fontsize=fs(9.5))
    ax.set_ylabel("Maximum temporal effect size", fontsize=fs(9.5))
    ax.set_xlim(k_summary["stability_core"].min() - 0.01, 1.005)
    ax.set_ylim(-0.015, k_summary["max_eta_timepoint"].max() + 0.055)
    ax.grid(color="#d9d9d9", alpha=0.45, linewidth=0.8)
    ax.tick_params(axis="both", labelsize=fs(8))
    _clean_axis(ax)


def draw_k_ranking_panel(ax: plt.Axes, k_summary: pd.DataFrame) -> None:
    ax.set_title(
        "C2. Revised K ranking",
        loc="left",
        fontsize=fs(12.0),
        fontweight="bold",
        pad=8,
    )
    k_summary = k_summary.sort_values("K").reset_index(drop=True)
    selected = k_summary["selected_for_case_study"].astype(bool)
    colors = ["#c44e52" if x else "#55a868" for x in selected]
    x = range(len(k_summary))

    ax.bar(x, k_summary["goal_aligned_score"], color=colors, edgecolor="white", linewidth=0.8)
    ax.set_ylabel("Goal-aligned score", fontsize=fs(9.5))
    ax.set_xticks(list(x))
    ax.set_xticklabels([f"K={int(k)}" for k in k_summary["K"]], rotation=70, ha="right", fontsize=fs(8))
    ax.tick_params(axis="y", labelsize=fs(8))
    ax.set_ylim(0, max(k_summary["goal_aligned_score"]) * 1.12)
    ax.grid(axis="y", color="#e6edf3", linewidth=0.8)
    ax.set_axisbelow(True)

    if selected.any():
        selected_idx = int(k_summary.index[selected][0])
        selected_score = k_summary.loc[selected_idx, "goal_aligned_score"]
        ax.text(
            selected_idx,
            selected_score + 0.02,
            "selected",
            ha="center",
            va="bottom",
            fontsize=fs(8),
            color="#263238",
        )

    _clean_axis(ax)


def draw_heatmap_panel(
    ax: plt.Axes,
    matrix_df: pd.DataFrame,
    *,
    title: str,
    xlabel: str,
    highlight_pattern: str | None = "Pattern3",
    cmap: str = "viridis",
) -> None:
    ax.set_title(title, loc="left", fontsize=fs(12.5), fontweight="bold", pad=8)

    matrix = matrix_df.set_index("pattern")
    time_cols = _order_timepoints(matrix.columns.tolist())
    if time_cols:
        matrix = matrix.loc[:, time_cols]

    values = matrix.to_numpy(dtype=float)
    im = ax.imshow(values, aspect="auto", cmap=cmap)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=fs(8.5))
    ax.set_yticks(range(matrix.shape[0]))
    y_labels = [
        "Pattern3\n(IFN)" if pattern == "Pattern3" else pattern
        for pattern in matrix.index.tolist()
    ]
    ax.set_yticklabels(y_labels, fontsize=fs(8.5))
    ax.set_xlabel(xlabel, fontsize=fs(10))
    ax.set_ylabel("CoGAPS pattern", fontsize=fs(10))

    if highlight_pattern in matrix.index:
        row = list(matrix.index).index(highlight_pattern)
        ax.add_patch(
            Rectangle(
                (-0.5, row - 0.5),
                matrix.shape[1],
                1,
                fill=False,
                edgecolor="#fff7ed",
                linewidth=2.3,
            )
        )
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.025)
    cbar.set_label("Mean score", fontsize=fs(8.5))
    cbar.ax.tick_params(labelsize=fs(8))


def draw_ifn_panel(
    ax: plt.Axes,
    summary: pd.DataFrame,
    subjects: pd.DataFrame,
) -> None:
    ax.set_title(
        "E. Pattern 3 IFN-like peak",
        loc="left",
        fontsize=fs(12.0),
        fontweight="bold",
        pad=8,
    )

    summary = summary.sort_values("day_merged_numeric")
    subjects = subjects.sort_values("day_merged_numeric")

    for _, sub in subjects.groupby("subject", observed=True):
        ax.plot(
            sub["day_merged_numeric"],
            sub["mean_pattern_score"],
            color="#7b8794",
            lw=1.15,
            alpha=0.7,
            marker="o",
            ms=3.2,
        )

    ax.errorbar(
        summary["day_merged_numeric"],
        summary["mean_subject_score"],
        yerr=summary["se_subject_score"],
        color="#c23b22",
        lw=2.5,
        marker="o",
        ms=5.2,
        capsize=3,
        label="subject mean",
        zorder=5,
    )
    ax.axvspan(9.3, 15.2, color="#fee2e2", alpha=0.65, zorder=0)
    ax.text(
        12.3,
        summary["mean_subject_score"].max() * 1.05,
        "late acute peak",
        ha="center",
        va="bottom",
        fontsize=fs(9),
        color="#9f2f20",
        fontweight="bold",
    )

    labels = summary["timepoint_merged"].tolist()
    days = summary["day_merged_numeric"].tolist()
    ax.set_xticks(days)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=fs(8.5))
    ax.set_xlabel("Experimental day", fontsize=fs(10))
    ax.set_ylabel("Pattern 3 mean score", fontsize=fs(10))
    ax.grid(axis="y", color="#d9e2ec", lw=0.8)
    ax.legend(frameon=False, fontsize=fs(8), loc="upper left")
    _clean_axis(ax)


def main() -> None:
    k_summary = pd.read_csv(SOURCE_DIR / "figure_01_k_selection_source.csv")
    composition = pd.read_csv(SOURCE_DIR / "figure_02_discovery_set_composition_source.csv")
    pattern_by_time = pd.read_csv(SOURCE_DIR / "figure_03_pattern_usage_by_time_matrix.csv")
    ifn_summary = pd.read_csv(SOURCE_DIR / "figure_05_ifn_pattern_summary_trajectory_source.csv")
    ifn_subjects = pd.read_csv(SOURCE_DIR / "figure_05_ifn_pattern_subject_trajectory_source.csv")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.labelcolor": "#12263a",
            "axes.edgecolor": "#829ab1",
            "xtick.color": "#334e68",
            "ytick.color": "#334e68",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.titlesize": 13,
        }
    )

    fig = plt.figure(figsize=(13.8, 17.2), constrained_layout=False)
    gs = gridspec.GridSpec(
        4,
        4,
        figure=fig,
        height_ratios=[0.85, 1.65, 1.85, 2.55],
        hspace=0.36,
        wspace=0.66,
        top=0.885,
        bottom=0.105,
        left=0.075,
        right=0.965,
    )

    fig.suptitle(
        "Temporal immune programs in experimental dengue infection",
        fontsize=fs(18.5),
        fontweight="bold",
        color="#12263a",
        y=0.965,
    )
    fig.text(
        0.5,
        0.925,
        "Roadmap: teaching dataset -> discovery composition -> K sweep and stability -> temporal patterns -> IFN-associated peak.",
        ha="center",
        fontsize=fs(10.8),
        color="#486581",
    )

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    sweep_gs = gs[2, :].subgridspec(1, 2, wspace=0.28)
    ax_b1 = fig.add_subplot(sweep_gs[0, 0])
    ax_b2 = fig.add_subplot(sweep_gs[0, 1])
    result_gs = gs[3, :].subgridspec(1, 2, width_ratios=[1.24, 1.0], wspace=0.34)
    ax_d = fig.add_subplot(result_gs[0, 0])
    ax_e = fig.add_subplot(result_gs[0, 1])

    draw_dataset_panel(ax_a)
    draw_composition_panel(ax_b, composition)
    draw_k_temporal_resolution_panel(ax_b1, k_summary)
    draw_k_ranking_panel(ax_b2, k_summary)
    draw_heatmap_panel(
        ax_d,
        pattern_by_time,
        title="D. K=10 pattern usage across infection days",
        xlabel="Experimental day",
    )
    draw_ifn_panel(ax_e, ifn_summary, ifn_subjects)

    fig.text(
        0.075,
        0.026,
        "PBMC = peripheral blood mononuclear cell. IFN = interferon.",
        fontsize=fs(8.4),
        color="#627d98",
    )

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=190, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
