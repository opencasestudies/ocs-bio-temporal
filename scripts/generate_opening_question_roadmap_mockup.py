#!/usr/bin/env python3
"""Generate a question-centered mockup for the Case Study 5 opening figure.

This mockup keeps the current opening figure untouched. It previews a possible
replacement figure that maps the four research questions to the evidence types
learners will inspect later in the case study.
"""

from __future__ import annotations

from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
from matplotlib.patches import FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "data" / "processed" / "figures"
SOURCE_DIR = FIG_DIR / "source_tables"
OUT_PNG = FIG_DIR / "figure_00_opening_question_roadmap_mockup.png"
OUT_SVG = FIG_DIR / "figure_00_opening_question_roadmap_mockup.svg"

TIMEPOINTS = ["D0", "D2", "D4", "D6", "D8", "D10", "D14/15", "D28"]
CELL_TYPE_ORDER = ["T_cell", "NK_cell", "Monocyte", "B_cell", "Neutrophil", "Plasmablast"]
CELL_TYPE_COLORS = {
    "T_cell": "#5fbf9f",
    "NK_cell": "#f28e63",
    "Monocyte": "#8da0cb",
    "B_cell": "#df83bd",
    "Neutrophil": "#9ccc50",
    "Plasmablast": "#f2c94c",
}

NAVY = "#102a43"
SLATE = "#486581"
LIGHT_GRID = "#e6edf3"
BLUE = "#2f6fd6"
GREEN = "#07866f"
PURPLE = "#7c2bd6"
ORANGE = "#d9541e"
RED = "#c51643"
DARK = "#25364a"


def fs(size: float) -> float:
    return size * 1.32


def wrap(text: str, width: int = 54) -> str:
    return fill(text, width=width)


def draw_card(ax: plt.Axes, border: str, face: str = "white") -> None:
    ax.axis("off")
    ax.add_patch(
        FancyBboxPatch(
            (0, 0),
            1,
            1,
            transform=ax.transAxes,
            boxstyle="round,pad=0.012,rounding_size=0.025",
            facecolor=face,
            edgecolor=border,
            linewidth=1.25,
            clip_on=False,
            zorder=-1,
        )
    )


def style_plot(ax: plt.Axes, border: str, face: str = "white") -> None:
    ax.set_facecolor(face)
    for side in ax.spines:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color(border)
        ax.spines[side].set_linewidth(1.1)
    ax.tick_params(colors=SLATE)


def label_panel(ax: plt.Axes, label: str, title: str, color: str, subtitle: str | None = None) -> None:
    ax.text(
        0.035,
        0.915,
        label,
        transform=ax.transAxes,
        ha="left",
        va="center",
        color="white",
        fontsize=fs(8.8),
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.32,rounding_size=0.12", facecolor=color, edgecolor=color),
        clip_on=False,
    )
    ax.text(
        0.115,
        0.93,
        title,
        transform=ax.transAxes,
        ha="left",
        va="center",
        color=NAVY,
        fontsize=fs(12.4),
        fontweight="bold",
        clip_on=False,
    )
    if subtitle:
        ax.text(
            0.115,
            0.855,
            subtitle,
            transform=ax.transAxes,
            ha="left",
            va="center",
            color=SLATE,
            fontsize=fs(8.7),
            clip_on=False,
        )


def clean_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)


def draw_timeline(ax: plt.Axes) -> None:
    positions = np.arange(len(TIMEPOINTS))
    ax.set_xlim(-0.12, len(TIMEPOINTS) - 0.88)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.hlines(0.52, positions[0], positions[-1], color="#53606c", linewidth=2.2, zorder=1)
    ax.add_patch(Rectangle((4.66, 0.29), 1.70, 0.46, facecolor="#fde7e7", edgecolor="none", zorder=0))
    ax.text(
        5.50,
        0.84,
        "late acute\nwindow",
        ha="center",
        va="bottom",
        fontsize=fs(8.1),
        color="#9f2f20",
        fontweight="bold",
    )

    for label, x_pos in zip(TIMEPOINTS, positions):
        is_late = label in {"D10", "D14/15"}
        face = ORANGE if is_late else "#eaf4ff"
        edge = "#8f2a1a" if is_late else "#233f5f"
        text_color = "white" if is_late else NAVY
        ax.scatter(x_pos, 0.52, s=370, facecolor=face, edgecolor=edge, linewidth=1.7, zorder=3)
        ax.text(
            x_pos,
            0.52,
            label.replace("/", "/\n") if label == "D14/15" else label,
            ha="center",
            va="center",
            fontsize=fs(7.9),
            fontweight="bold",
            color=text_color,
            zorder=4,
        )


def draw_data_setting(ax: plt.Axes, composition: pd.DataFrame) -> None:
    label_panel(
        ax,
        "A",
        "Start with the biological setting",
        BLUE,
        "Use two anchors from the experimental DENV-1 samples: infection day and broad PBMC identity.",
    )
    draw_card(ax, "#bdd0e7", "#fbfdff")

    timeline_ax = ax.inset_axes([0.035, 0.47, 0.565, 0.30])
    draw_timeline(timeline_ax)

    comp_ax = ax.inset_axes([0.625, 0.43, 0.34, 0.34])
    comp = composition.copy()
    comp["timepoint_merged"] = pd.Categorical(comp["timepoint_merged"], categories=TIMEPOINTS, ordered=True)
    comp["broad_cell_type"] = pd.Categorical(comp["broad_cell_type"], categories=CELL_TYPE_ORDER, ordered=True)
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
    x = np.arange(len(pivot.index))
    bottom = np.zeros(len(pivot.index))
    for cell_type in CELL_TYPE_ORDER:
        values = pivot[cell_type].to_numpy()
        comp_ax.bar(
            x,
            values,
            bottom=bottom,
            width=0.74,
            color=CELL_TYPE_COLORS[cell_type],
            edgecolor="white",
            linewidth=0.4,
        )
        bottom += values
    comp_ax.set_xticks([0, 3, 5, 6, 7])
    comp_ax.set_xticklabels(["D0", "D6", "D10", "D14/15", "D28"], rotation=0, ha="center", fontsize=fs(7.1))
    comp_ax.set_yticks([])
    comp_ax.set_title("PBMC mixture retained", fontsize=fs(8.8), loc="left", color=NAVY, pad=2)
    clean_axis(comp_ax)

    legend_ax = ax.inset_axes([0.625, 0.12, 0.34, 0.22])
    legend_ax.axis("off")
    legend_ax.text(
        0.0,
        0.98,
        "Cell identity colors",
        ha="left",
        va="top",
        fontsize=fs(7.7),
        color=NAVY,
        fontweight="bold",
    )
    for i, cell_type in enumerate(CELL_TYPE_ORDER):
        col = i % 2
        row = i // 2
        x0 = col * 0.50
        y0 = 0.58 - row * 0.29
        legend_ax.add_patch(
            Rectangle(
                (x0, y0 - 0.06),
                0.055,
                0.12,
                transform=legend_ax.transAxes,
                facecolor=CELL_TYPE_COLORS[cell_type],
                edgecolor="white",
                linewidth=0.7,
            )
        )
        legend_ax.text(
            x0 + 0.075,
            y0,
            cell_type.replace("_", " "),
            transform=legend_ax.transAxes,
            ha="left",
            va="center",
            fontsize=fs(7.0),
            color=DARK,
        )

    ax.text(
        0.055,
        0.37,
        wrap(
            "Begin with a shared infection clock, then keep PBMC identity visible so a time signal is not mistaken for a changing cell mixture.",
            58,
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs(8.7),
        color=DARK,
    )


def draw_rank_panel(ax: plt.Axes, k_summary: pd.DataFrame) -> None:
    draw_card(ax, "#f0c4aa", "#fffaf5")
    label_panel(
        ax,
        "Q1",
        "Choose a model resolution",
        ORANGE,
        "Check rank evidence before naming patterns.",
    )
    plot_ax = ax.inset_axes([0.085, 0.145, 0.88, 0.65])
    style_plot(plot_ax, "#f0c4aa", "#fffaf5")

    k_summary = k_summary.sort_values("K")
    colors = []
    for _, row in k_summary.iterrows():
        if bool(row["selected_for_case_study"]):
            colors.append(ORANGE)
        elif bool(row["on_stability_plateau"]) and row["max_eta_timepoint"] < 0.12:
            colors.append(BLUE)
        else:
            colors.append("#b7c0ca")

    plot_ax.scatter(
        k_summary["stability_core"],
        k_summary["max_eta_timepoint"],
        s=86,
        c=colors,
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
    )
    plot_ax.set_xlabel("Seed-to-seed stability", fontsize=fs(8.5), labelpad=4)
    plot_ax.set_ylabel("Temporal structure", fontsize=fs(8.5), labelpad=5)
    plot_ax.tick_params(axis="both", labelsize=fs(7.4), pad=2)
    plot_ax.grid(color=LIGHT_GRID, linewidth=0.8)
    plot_ax.set_xlim(k_summary["stability_core"].min() - 0.012, 1.004)
    plot_ax.set_ylim(-0.02, k_summary["max_eta_timepoint"].max() + 0.07)
    clean_axis(plot_ax)


def draw_temporal_panel(ax: plt.Axes, pattern_by_time: pd.DataFrame) -> None:
    draw_card(ax, "#b7dccf", "#f8fffc")
    label_panel(
        ax,
        "Q2",
        "Find programs that change over time",
        GREEN,
        "Summaries show which patterns rise, peak, or fade.",
    )
    plot_ax = ax.inset_axes([0.095, 0.145, 0.82, 0.65])
    style_plot(plot_ax, "#b7dccf", "#f8fffc")

    matrix = pattern_by_time.set_index("pattern").loc[:, TIMEPOINTS]
    values = matrix.to_numpy(dtype=float)
    im = plot_ax.imshow(values, aspect="auto", cmap="viridis")
    plot_ax.set_xticks(range(len(TIMEPOINTS)))
    plot_ax.set_xticklabels(TIMEPOINTS, rotation=45, ha="right", fontsize=fs(7.1))
    plot_ax.set_yticks(range(matrix.shape[0]))
    plot_ax.set_yticklabels([f"P{i + 1}" for i in range(matrix.shape[0])], fontsize=fs(7.1))
    plot_ax.set_xlabel("Experimental day", fontsize=fs(8.5), labelpad=4)
    plot_ax.set_ylabel("Model pattern", fontsize=fs(8.5), labelpad=5)
    plot_ax.add_patch(Rectangle((-0.5, 2 - 0.5), matrix.shape[1], 1, fill=False, edgecolor="white", linewidth=2.4))
    cbar = plt.colorbar(im, ax=plot_ax, fraction=0.046, pad=0.025)
    cbar.set_label("Mean score", fontsize=fs(7.3))
    cbar.ax.tick_params(labelsize=fs(6.9), pad=2)


def draw_identity_activity_panel(ax: plt.Axes, annotation: pd.DataFrame) -> None:
    draw_card(ax, "#d4c2ee", "#fcfaff")
    label_panel(
        ax,
        "Q3",
        "Separate identity from activity",
        PURPLE,
        "Compare infection timing with PBMC identity.",
    )
    plot_ax = ax.inset_axes([0.095, 0.145, 0.87, 0.65])
    style_plot(plot_ax, "#d4c2ee", "#fcfaff")

    classes = annotation["pattern_class"].fillna("mixed").astype(str)
    color_map = {
        "activity-like": GREEN,
        "identity-like": BLUE,
        "mixed": PURPLE,
    }
    colors = [color_map.get(x, "#7b8794") for x in classes]
    plot_ax.scatter(
        annotation["eta_broad_cell_type"],
        annotation["eta_timepoint"],
        s=112,
        c=colors,
        edgecolor="white",
        linewidth=0.8,
        alpha=0.95,
        zorder=3,
    )
    plot_ax.plot([0, 0.52], [0, 0.52], color="#a5b1c2", linewidth=1.2, linestyle="--", zorder=1)
    plot_ax.set_xlim(-0.02, max(0.53, annotation["eta_broad_cell_type"].max() + 0.04))
    plot_ax.set_ylim(-0.02, max(0.52, annotation["eta_timepoint"].max() + 0.04))
    plot_ax.set_xlabel("Broad PBMC identity structure", fontsize=fs(8.5), labelpad=4)
    plot_ax.set_ylabel("Infection-day structure", fontsize=fs(8.5), labelpad=5)
    plot_ax.tick_params(axis="both", labelsize=fs(7.2), pad=2)
    plot_ax.grid(color=LIGHT_GRID, linewidth=0.8)
    clean_axis(plot_ax)
    plot_ax.text(
        0.98,
        0.93,
        "more\nactivity-like",
        transform=plot_ax.transAxes,
        ha="right",
        va="top",
        fontsize=fs(7.6),
        color=GREEN,
        fontweight="bold",
    )
    plot_ax.text(
        0.03,
        0.22,
        "more\nidentity-like",
        transform=plot_ax.transAxes,
        ha="left",
        va="center",
        fontsize=fs(7.6),
        color=BLUE,
        fontweight="bold",
    )


def draw_response_panel(ax: plt.Axes, summary: pd.DataFrame) -> None:
    draw_card(ax, "#f3bdcc", "#fff8fb")
    label_panel(
        ax,
        "Q4",
        "Evaluate the response claim",
        RED,
        "Combine timing, genes, direction, and sensitivity evidence.",
    )

    spark_ax = ax.inset_axes([0.095, 0.30, 0.47, 0.48])
    summary = summary.sort_values("day_merged_numeric")
    spark_ax.plot(
        summary["day_merged_numeric"],
        summary["mean_subject_score"],
        color=RED,
        linewidth=2.7,
        marker="o",
        markersize=5.2,
    )
    spark_ax.fill_between(
        summary["day_merged_numeric"],
        summary["mean_subject_score"] - summary["se_subject_score"],
        summary["mean_subject_score"] + summary["se_subject_score"],
        color=RED,
        alpha=0.14,
    )
    spark_ax.axvspan(9.3, 15.2, color="#fee2e2", alpha=0.7)
    spark_ax.set_xticks([0, 6, 10, 14.5, 28])
    spark_ax.set_xticklabels(["D0", "D6", "D10", "D14/15", "D28"], rotation=45, ha="right", fontsize=fs(7.0))
    spark_ax.set_yticks([])
    spark_ax.set_title("candidate trajectory", fontsize=fs(8.1), loc="left", color=NAVY, pad=3)
    clean_axis(spark_ax)

    checks = [
        ("Timing", "late acute rise"),
        ("Genes", "IFN-associated signal"),
        ("Direction", "higher or lower\nvs baseline"),
        ("Sensitivity", "interpretation is bounded"),
    ]
    y = 0.73
    for label, detail in checks:
        ax.text(
            0.65,
            y,
            label,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=fs(9.0),
            fontweight="bold",
            color=NAVY,
        )
        ax.text(
            0.65,
            y - 0.085,
            detail,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=fs(7.9),
            color=SLATE,
        )
        ax.scatter(0.60, y - 0.012, s=88, transform=ax.transAxes, color=RED, edgecolor="white", linewidth=0.8)
        y -= 0.18

    ax.text(
        0.095,
        0.16,
        wrap("Use all four checks; one plot is not enough.", 34),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs(8.2),
        color=DARK,
    )


def main() -> None:
    k_summary = pd.read_csv(SOURCE_DIR / "figure_01_k_selection_source.csv")
    composition = pd.read_csv(SOURCE_DIR / "figure_02_discovery_set_composition_source.csv")
    pattern_by_time = pd.read_csv(SOURCE_DIR / "figure_03_pattern_usage_by_time_matrix.csv")
    annotation = pd.read_csv(ROOT / "data" / "processed" / "interpretation" / "pattern_annotation_table.csv")
    ifn_summary = pd.read_csv(SOURCE_DIR / "figure_05_ifn_pattern_summary_trajectory_source.csv")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.labelcolor": NAVY,
            "xtick.color": SLATE,
            "ytick.color": SLATE,
        }
    )

    fig = plt.figure(figsize=(12.4, 13.8), constrained_layout=False)
    gs = gridspec.GridSpec(
        3,
        2,
        figure=fig,
        height_ratios=[0.82, 1.30, 1.30],
        hspace=0.14,
        wspace=0.08,
        top=0.905,
        bottom=0.065,
        left=0.055,
        right=0.965,
    )

    fig.suptitle(
        "Roadmap: how the case study turns PBMC data into CoGAPS evidence",
        x=0.5,
        y=0.962,
        fontsize=fs(17.0),
        color=NAVY,
        fontweight="bold",
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_q1 = fig.add_subplot(gs[1, 0])
    ax_q2 = fig.add_subplot(gs[1, 1])
    ax_q3 = fig.add_subplot(gs[2, 0])
    ax_q4 = fig.add_subplot(gs[2, 1])

    draw_data_setting(ax_a, composition)
    draw_rank_panel(ax_q1, k_summary)
    draw_temporal_panel(ax_q2, pattern_by_time)
    draw_identity_activity_panel(ax_q3, annotation)
    draw_response_panel(ax_q4, ifn_summary)

    fig.text(
        0.075,
        0.027,
        "Roadmap key: Q1 rank decision; Q2 temporal programs; Q3 identity versus activity; Q4 response-claim evidence. PBMC = peripheral blood mononuclear cell.",
        ha="left",
        va="bottom",
        fontsize=fs(8.0),
        color="#627d98",
    )

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
