#!/usr/bin/env python3
"""Generate the opening composite figure for Case Study 5.

The figure is intentionally generated from lightweight, GitHub-safe source
tables that are already part of the rendered case study. It previews the
scientific arc without requiring learners to rerun CoGAPS.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "data" / "processed" / "figures"
SOURCE_DIR = FIG_DIR / "source_tables"

OUT_PNG = FIG_DIR / "figure_00_opening_summary.png"
OUT_SVG = FIG_DIR / "figure_00_opening_summary.svg"


def _order_timepoints(labels: list[str]) -> list[str]:
    preferred = ["D0", "D2", "D4", "D6", "D8", "D10", "D14/15", "D28"]
    return [x for x in preferred if x in labels]


def draw_heatmap_panel(ax: plt.Axes, pattern_matrix: pd.DataFrame) -> None:
    ax.set_title(
        "A. Mean CoGAPS pattern usage across the dengue time course",
        loc="left",
        fontsize=12.5,
        fontweight="bold",
        pad=8,
    )

    matrix = pattern_matrix.set_index("pattern")
    timepoints = _order_timepoints([c for c in matrix.columns])
    matrix = matrix.loc[:, timepoints]

    values = matrix.to_numpy(dtype=float)
    im = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(timepoints)))
    ax.set_xticklabels(timepoints, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(matrix.shape[0]))
    y_labels = [
        "Pattern3\n(IFN)" if pattern == "Pattern3" else pattern
        for pattern in matrix.index.tolist()
    ]
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.set_xlabel("Experimental day", fontsize=10)
    ax.set_ylabel("CoGAPS pattern", fontsize=10)

    if "Pattern3" in matrix.index:
        row = list(matrix.index).index("Pattern3")
        ax.add_patch(
            Rectangle(
                (-0.5, row - 0.5),
                len(timepoints),
                1,
                fill=False,
                edgecolor="#fef3c7",
                linewidth=2.2,
            )
        )
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Mean cell score", fontsize=9)
    cbar.ax.tick_params(labelsize=8)


def draw_ifn_panel(
    ax: plt.Axes,
    summary: pd.DataFrame,
    subjects: pd.DataFrame,
) -> None:
    ax.set_title(
        "B. Pattern 3 peaks during the late acute window",
        loc="left",
        fontsize=12.5,
        fontweight="bold",
        pad=8,
    )

    summary = summary.sort_values("day_merged_numeric")
    subjects = subjects.sort_values("day_merged_numeric")

    for subject, sub in subjects.groupby("subject", observed=True):
        ax.plot(
            sub["day_merged_numeric"],
            sub["mean_pattern_score"],
            color="#7b8794",
            lw=1.3,
            alpha=0.65,
            marker="o",
            ms=3.5,
            label=str(subject),
        )

    ax.errorbar(
        summary["day_merged_numeric"],
        summary["mean_subject_score"],
        yerr=summary["se_subject_score"],
        color="#c23b22",
        lw=2.4,
        marker="o",
        ms=5.5,
        capsize=3,
        label="Subject mean",
        zorder=5,
    )
    ax.axvspan(9.3, 15.2, color="#fee2e2", alpha=0.65, zorder=0)

    labels = summary["timepoint_merged"].tolist()
    days = summary["day_merged_numeric"].tolist()
    ax.set_xticks(days)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("Experimental day", fontsize=10)
    ax.set_ylabel("Pattern 3 mean score", fontsize=10)
    ax.grid(axis="y", color="#d9e2ec", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left")


def main() -> None:
    pattern_matrix = pd.read_csv(SOURCE_DIR / "figure_03_pattern_usage_by_time_matrix.csv")
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
        }
    )

    fig, (ax_a, ax_b) = plt.subplots(
        2,
        1,
        figsize=(11.0, 10.0),
        gridspec_kw={"height_ratios": [1.05, 1.0], "hspace": 0.48},
        constrained_layout=False,
    )

    draw_heatmap_panel(ax_a, pattern_matrix)
    draw_ifn_panel(ax_b, ifn_summary, ifn_subjects)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=220, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
