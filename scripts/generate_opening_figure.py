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


def _order_timepoints(labels: list[str]) -> list[str]:
    return [x for x in TIMEPOINTS if x in labels]


def _clean_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)


def draw_dataset_panel(ax: plt.Axes) -> None:
    ax.set_title(
        "A. Experimental PBMC discovery time course",
        loc="left",
        fontsize=12.5,
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
        fontsize=9.5,
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
            fontsize=9.5,
            fontweight="bold",
            color=text_color,
            zorder=4,
        )

    ax.text(
        positions[0],
        0.96,
        "Three dengue-naive adults sampled from baseline to recovery",
        fontsize=10.8,
        color="#243b53",
        va="top",
    )
    ax.text(
        positions[0],
        0.06,
        "Balanced teaching subset: 600 cells per day, preserving broad PBMC identities",
        fontsize=10.2,
        color="#486581",
        va="bottom",
    )


def draw_k_panel(ax: plt.Axes, k_summary: pd.DataFrame) -> None:
    ax.set_title(
        "B. Sweep evidence for K = 10",
        loc="left",
        fontsize=12.5,
        fontweight="bold",
        pad=8,
    )
    k_summary = k_summary.sort_values("K")
    selected = k_summary.loc[k_summary["selected_for_case_study"].astype(bool)]

    plateau = k_summary.loc[k_summary["on_stability_plateau"].astype(bool)]
    if len(plateau) > 0:
        ax.scatter(
            plateau["K"],
            plateau["stability_core"],
            s=70,
            facecolor="#e0f2fe",
            edgecolor="#0f73a8",
            linewidth=1.4,
            label="stability plateau",
            zorder=2,
        )

    ax.plot(
        k_summary["K"],
        k_summary["stability_core"],
        color="#2878b5",
        marker="o",
        linewidth=2,
        markersize=4,
        label="stability",
        zorder=3,
    )
    ax.plot(
        k_summary["K"],
        k_summary["goal_aligned_score"],
        color="#c23b22",
        marker="s",
        linewidth=2,
        markersize=4,
        label="goal-aligned score",
        zorder=3,
    )

    if len(selected) == 1:
        row = selected.iloc[0]
        ax.axvline(row["K"], color="#1f2933", linestyle="--", linewidth=1.2)
        ax.scatter(
            [row["K"]],
            [row["goal_aligned_score"]],
            s=150,
            facecolor="#c23b22",
            edgecolor="white",
            linewidth=1.8,
            zorder=5,
        )
        ax.text(
            row["K"] + 0.6,
            0.7,
            "selected\nK = 10",
            color="#9f2f20",
            fontsize=10.5,
            fontweight="bold",
            va="bottom",
        )

    ax.text(
        0.02,
        0.04,
        "K = number of CoGAPS patterns\nK=5 is stable but compresses biology\nK=10 stays stable and resolves a temporal program",
        transform=ax.transAxes,
        fontsize=8.7,
        color="#486581",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#bcccdc", alpha=0.92),
    )
    ax.set_xlabel("Model rank K", fontsize=10)
    ax.set_ylabel("Score", fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([5, 10, 16, 24, 32, 40])
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="y", color="#e6edf3", linewidth=0.8)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
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
    ax.set_title(title, loc="left", fontsize=12.5, fontweight="bold", pad=8)

    matrix = matrix_df.set_index("pattern")
    time_cols = _order_timepoints(matrix.columns.tolist())
    if time_cols:
        matrix = matrix.loc[:, time_cols]

    values = matrix.to_numpy(dtype=float)
    im = ax.imshow(values, aspect="auto", cmap=cmap)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8.5)
    ax.set_yticks(range(matrix.shape[0]))
    y_labels = [
        "Pattern3\n(IFN)" if pattern == "Pattern3" else pattern
        for pattern in matrix.index.tolist()
    ]
    ax.set_yticklabels(y_labels, fontsize=8.5)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("CoGAPS pattern", fontsize=10)

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
    cbar.set_label("Mean score", fontsize=8.5)
    cbar.ax.tick_params(labelsize=8)


def draw_ifn_panel(
    ax: plt.Axes,
    summary: pd.DataFrame,
    subjects: pd.DataFrame,
) -> None:
    ax.set_title(
        "C. Pattern 3 late acute trajectory",
        loc="left",
        fontsize=12.5,
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
        fontsize=9,
        color="#9f2f20",
        fontweight="bold",
    )

    labels = summary["timepoint_merged"].tolist()
    days = summary["day_merged_numeric"].tolist()
    ax.set_xticks(days)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8.5)
    ax.set_xlabel("Experimental day", fontsize=10)
    ax.set_ylabel("Pattern 3 mean score", fontsize=10)
    ax.grid(axis="y", color="#d9e2ec", lw=0.8)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    _clean_axis(ax)


def main() -> None:
    k_summary = pd.read_csv(SOURCE_DIR / "figure_01_k_selection_source.csv")
    pattern_by_time = pd.read_csv(SOURCE_DIR / "figure_03_pattern_usage_by_time_matrix.csv")
    pattern_by_cell_type = pd.read_csv(
        SOURCE_DIR / "figure_04_pattern_usage_by_cell_type_matrix.csv"
    )
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

    fig = plt.figure(figsize=(13.8, 12.6), constrained_layout=False)
    gs = gridspec.GridSpec(
        4,
        4,
        figure=fig,
        height_ratios=[0.82, 1.18, 1.3, 1.18],
        hspace=0.78,
        wspace=0.72,
        top=0.88,
        bottom=0.075,
        left=0.075,
        right=0.965,
    )

    fig.suptitle(
        "Temporal immune programs in experimental dengue infection",
        fontsize=18.5,
        fontweight="bold",
        color="#12263a",
        y=0.965,
    )
    fig.text(
        0.5,
        0.925,
        "Roadmap: teaching dataset -> K sweep and stability -> temporal and identity patterns -> IFN-associated peak.",
        ha="center",
        fontsize=10.8,
        color="#486581",
    )

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :2])
    ax_e = fig.add_subplot(gs[1, 2:])
    ax_c = fig.add_subplot(gs[2, :])
    ax_d = fig.add_subplot(gs[3, :])

    draw_dataset_panel(ax_a)
    draw_k_panel(ax_b, k_summary)
    draw_ifn_panel(ax_e, ifn_summary, ifn_subjects)
    draw_heatmap_panel(
        ax_c,
        pattern_by_time,
        title="D. Selected K = 10 patterns across infection days",
        xlabel="Experimental day",
    )
    draw_heatmap_panel(
        ax_d,
        pattern_by_cell_type,
        title="E. The same patterns also track broad PBMC identity",
        xlabel="Broad immune-cell type",
    )

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=190, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
