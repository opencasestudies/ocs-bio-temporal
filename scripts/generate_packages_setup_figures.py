from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
IMG_DIR = ROOT / "img"

RESOURCE_PNG = IMG_DIR / "packages_setup_resources.png"
EVIDENCE_PNG = IMG_DIR / "packages_setup_workflow_map_core_path_colors.png"

RESOURCE_W, RESOURCE_H = 1800, 640
EVIDENCE_W, EVIDENCE_H = 1800, 3220
BG = "#F7FAFD"
WHITE = "#FFFFFF"
NAVY = "#08254D"
INK = "#17202A"
MUTED = "#4B6178"
PANEL_LINE = "#CBD7E2"
BLUE = "#1D4ED8"
TEAL = "#0F766E"
PURPLE = "#7E22CE"
ORANGE = "#C2410C"
INDIGO = "#4338CA"
GREEN = "#047857"
RED = "#BE123C"
SLATE = "#475569"
AMBER = "#B45309"

# Colorblind-safer path palette.
# Core lesson steps use saturated colors; upstream/sweep-only steps use quieter
# hues so learners can see which evidence they directly work with.
CORE_ORANGE = "#D55E00"
CORE_BLUE = "#0072B2"
CORE_GREEN = "#009E73"
CORE_GOLD = "#E69F00"
MUTED_ROSE = "#8A737D"
UPSTREAM_BLUE = "#6B7785"
UPSTREAM_GREEN = "#71827A"
SWEEP_PURPLE = "#7D748C"

FONT_REG_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
]
FONT_BOLD_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
]


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = FONT_BOLD_CANDIDATES if bold else FONT_REG_CANDIDATES
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size)
    return ImageFont.load_default(size=size)


H1 = font(61, True)
H2 = font(37, True)
H3 = font(30, True)
SUBTITLE = font(24)
BODY = font(24)
BODY_BOLD = font(24, True)
SMALL = font(20)
SMALL_BOLD = font(20, True)
TINY = font(16)
CHIP = font(18)
CHIP_BOLD = font(18, True)
SCRIPT = font(17, True)

SETUP = (48, 48, 1768, 598)
LEFT_PANEL = (24, 40, 930, 3020)
RIGHT_PANEL = (980, 40, 1772, 3020)
DIVIDER_X = 955
FOOTER = (48, 3062, 1768, 3180)

MAIN_LESSON_MARKER_STEPS = {4, 5, 6, 7, 8}


def rgb(hex_value: str) -> tuple[int, int, int]:
    hex_value = hex_value.lstrip("#")
    return tuple(int(hex_value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(hex_value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*rgb(hex_value), alpha)


def tint(hex_value: str, amount: float = 0.9) -> str:
    r, g, b = rgb(hex_value)
    rr = round(r * (1 - amount) + 255 * amount)
    gg = round(g * (1 - amount) + 255 * amount)
    bb = round(b * (1 - amount) + 255 * amount)
    return f"#{rr:02x}{gg:02x}{bb:02x}"


def text_size(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def wrap(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont, max_width: int) -> list[str]:
    lines: list[str] = []
    for para in text.split("\n"):
        words = para.split()
        if not words:
            lines.append("")
            continue
        line = words[0]
        for word in words[1:]:
            trial = f"{line} {word}"
            if text_size(draw, trial, fnt)[0] <= max_width:
                line = trial
            else:
                lines.append(line)
                line = word
        lines.append(line)
    return lines


class Audit:
    def __init__(self) -> None:
        self.items: list[dict[str, object]] = []
        self.failures: list[dict[str, object]] = []

    def add(self, name: str, box: tuple[int, int, int, int], parent: tuple[int, int, int, int] | None = None) -> None:
        self.items.append({"name": name, "box": list(box)})
        if parent:
            x0, y0, x1, y1 = box
            px0, py0, px1, py1 = parent
            if x0 < px0 or y0 < py0 or x1 > px1 or y1 > py1:
                self.failures.append({"type": "bounds", "name": name, "box": list(box), "parent": list(parent)})

    @staticmethod
    def overlap(a: tuple[int, int, int, int], b: tuple[int, int, int, int], pad: int = 0) -> bool:
        ax0, ay0, ax1, ay1 = a
        bx0, by0, bx1, by1 = b
        return not (ax1 + pad <= bx0 or bx1 + pad <= ax0 or ay1 + pad <= by0 or by1 + pad <= ay0)

    def no_overlap(self, group: str, boxes: list[tuple[str, tuple[int, int, int, int]]], pad: int = 0) -> None:
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                if self.overlap(boxes[i][1], boxes[j][1], pad=pad):
                    self.failures.append(
                        {
                            "type": "overlap",
                            "group": group,
                            "a": boxes[i][0],
                            "a_box": list(boxes[i][1]),
                            "b": boxes[j][0],
                            "b_box": list(boxes[j][1]),
                        }
                    )


def draw_text_block(
    draw: ImageDraw.ImageDraw,
    audit: Audit,
    name: str,
    x: int,
    y: int,
    text: str,
    fnt: ImageFont.FreeTypeFont,
    max_width: int,
    fill: str,
    parent: tuple[int, int, int, int] | None,
    line_gap: int = 7,
) -> tuple[int, tuple[int, int, int, int]]:
    start = y
    max_x = x
    line_boxes: list[tuple[str, tuple[int, int, int, int]]] = []
    for idx, line in enumerate(wrap(draw, text, fnt, max_width)):
        w, h = text_size(draw, line or "Ag", fnt)
        box = (x, y, x + w, y + h)
        draw.text((x, y), line, font=fnt, fill=fill)
        audit.add(f"{name}-line-{idx}", box, parent)
        line_boxes.append((f"{name}-line-{idx}", box))
        max_x = max(max_x, x + w)
        y += h + line_gap
    full = (x, start, max_x, y - line_gap if y > start else start)
    audit.add(name, full, parent)
    audit.no_overlap(name, line_boxes, pad=0)
    return y, full


def shadow_rect(
    image: Image.Image,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str,
    outline: str,
    width: int = 2,
    shadow_alpha: int = 14,
    shadow_blur: int = 14,
    shadow_offset: tuple[int, int] = (0, 8),
) -> None:
    layer = Image.new("RGBA", image.size, (0, 0, 0, 0))
    ld = ImageDraw.Draw(layer)
    x0, y0, x1, y1 = box
    ox, oy = shadow_offset
    ld.rounded_rectangle((x0 + ox, y0 + oy, x1 + ox, y1 + oy), radius=radius, fill=(15, 23, 42, shadow_alpha))
    layer = layer.filter(ImageFilter.GaussianBlur(shadow_blur))
    image.alpha_composite(layer)
    draw = ImageDraw.Draw(image, "RGBA")
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def circle_icon_base(draw: ImageDraw.ImageDraw, x: int, y: int, color: str, r: int = 41) -> tuple[int, int, int, int]:
    box = (x, y, x + 2 * r, y + 2 * r)
    draw.ellipse(box, fill=tint(color, 0.88), outline=color, width=3)
    return box


def draw_setup_icon(draw: ImageDraw.ImageDraw, x: int, y: int, kind: str, color: str) -> None:
    box = (x, y, x + 76, y + 76)
    draw.ellipse(box, fill=color)
    cx, cy = x + 38, y + 38
    if kind == "docker":
        pts = [
            (cx, cy - 24),
            (cx + 24, cy - 10),
            (cx + 24, cy + 16),
            (cx, cy + 30),
            (cx - 24, cy + 16),
            (cx - 24, cy - 10),
        ]
        draw.line(pts + [pts[0]], fill=WHITE, width=4, joint="curve")
        draw.line((cx, cy - 24, cx, cy + 4), fill=WHITE, width=3)
        draw.line((cx - 24, cy - 10, cx, cy + 4, cx + 24, cy - 10), fill=WHITE, width=3, joint="curve")
    elif kind == "repo":
        draw.text((cx - 15, cy - 10), "GH", font=font(20, True), fill=WHITE)
    else:
        draw.rounded_rectangle((cx - 22, cy - 13, cx + 22, cy + 22), radius=3, outline=WHITE, width=4)
        draw.rectangle((cx - 14, cy - 25, cx + 14, cy - 14), outline=WHITE, width=4)
        draw.line((cx - 14, cy - 1, cx + 14, cy - 1), fill=WHITE, width=4)
        draw.line((cx - 12, cy + 10, cx + 12, cy + 10), fill=WHITE, width=4)


def draw_panel_icon(draw: ImageDraw.ImageDraw, x: int, y: int, kind: str, color: str) -> None:
    if kind == "network":
        pts = [(x + 20, y + 9), (x + 7, y + 42), (x + 43, y + 42)]
        draw.line(pts + [pts[0]], fill=color, width=4)
        for px, py in pts:
            draw.ellipse((px - 9, py - 9, px + 9, py + 9), fill=WHITE, outline=color, width=4)
    else:
        cx, cy = x + 25, y + 25
        draw.ellipse((cx - 21, cy - 21, cx + 21, cy + 21), outline=color, width=5)
        draw.ellipse((cx - 8, cy - 8, cx + 8, cy + 8), outline=color, width=4)
        draw.line((cx - 32, cy, cx - 20, cy), fill=color, width=4)
        draw.line((cx + 20, cy, cx + 32, cy), fill=color, width=4)
        draw.line((cx, cy - 32, cx, cy - 20), fill=color, width=4)
        draw.line((cx, cy + 20, cx, cy + 32), fill=color, width=4)


def draw_step_glyph(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], kind: str, color: str) -> None:
    x0, y0, x1, y1 = box
    cx, cy = (x0 + x1) // 2, (y0 + y1) // 2
    draw.ellipse(box, fill=tint(color, 0.88), outline=color, width=2)
    if kind == "file":
        draw.rectangle((cx - 15, cy - 22, cx + 16, cy + 24), outline=color, width=3)
        draw.line((cx + 5, cy - 22, cx + 16, cy - 11), fill=color, width=3)
        draw.line((cx - 8, cy - 4, cx + 9, cy - 4), fill=color, width=3)
        draw.line((cx - 8, cy + 8, cx + 9, cy + 8), fill=color, width=3)
    elif kind == "filter":
        funnel = [(cx - 26, cy - 18), (cx + 26, cy - 18), (cx + 7, cy + 4), (cx + 7, cy + 22), (cx - 7, cy + 29), (cx - 7, cy + 4), (cx - 26, cy - 18)]
        draw.line(funnel, fill=color, width=3, joint="curve")
    elif kind == "bars":
        for dx, h in [(-18, 22), (0, 39), (18, 29)]:
            draw.rounded_rectangle((cx + dx - 5, cy + 24 - h, cx + dx + 5, cy + 24), radius=3, fill=color)
        draw.line((cx - 28, cy + 24, cx + 30, cy + 24), fill=color, width=3)
    elif kind == "scale":
        beam_y = cy - 13
        pan_y = cy + 12
        draw.line((cx, cy - 27, cx, cy + 23), fill=color, width=3)
        draw.ellipse((cx - 5, cy - 31, cx + 5, cy - 21), fill=color)
        draw.line((cx - 30, beam_y, cx + 30, beam_y), fill=color, width=4)
        draw.line((cx - 20, beam_y, cx - 31, pan_y), fill=color, width=3)
        draw.line((cx - 20, beam_y, cx - 9, pan_y), fill=color, width=3)
        draw.line((cx + 20, beam_y, cx + 9, pan_y), fill=color, width=3)
        draw.line((cx + 20, beam_y, cx + 31, pan_y), fill=color, width=3)
        draw.line((cx - 36, pan_y, cx - 4, pan_y), fill=color, width=3)
        draw.line((cx + 4, pan_y, cx + 36, pan_y), fill=color, width=3)
        draw.line((cx - 34, pan_y, cx - 28, cy + 24, cx - 12, cy + 24, cx - 6, pan_y), fill=color, width=3, joint="curve")
        draw.line((cx + 6, pan_y, cx + 12, cy + 24, cx + 28, cy + 24, cx + 34, pan_y), fill=color, width=3, joint="curve")
        draw.line((cx - 17, cy + 25, cx + 17, cy + 25), fill=color, width=4)
        draw.line((cx - 27, cy + 31, cx + 27, cy + 31), fill=color, width=4)
    elif kind == "download":
        draw.line((cx, cy - 26, cx, cy + 9), fill=color, width=4)
        draw.line((cx - 15, cy - 6, cx, cy + 9, cx + 15, cy - 6), fill=color, width=4, joint="curve")
        draw.rectangle((cx - 26, cy + 18, cx + 26, cy + 29), outline=color, width=3)
    elif kind == "pie":
        draw.pieslice((cx - 25, cy - 25, cx + 25, cy + 25), 270, 360, fill=color)
        draw.pieslice((cx - 25, cy - 25, cx + 25, cy + 25), 0, 270, outline=color, width=3)
        draw.line((cx, cy, cx, cy - 25), fill=color, width=3)
        draw.line((cx, cy, cx + 25, cy), fill=color, width=3)
    elif kind == "direction":
        draw.line((cx - 12, cy + 23, cx - 12, cy - 22), fill=color, width=4)
        draw.line((cx - 26, cy - 8, cx - 12, cy - 22, cx + 2, cy - 8), fill=color, width=4, joint="curve")
        draw.line((cx + 18, cy - 23, cx + 18, cy + 22), fill=color, width=4)
        draw.line((cx + 4, cy + 8, cx + 18, cy + 22, cx + 32, cy + 8), fill=color, width=4, joint="curve")
    else:
        draw.rectangle((cx - 24, cy - 20, cx + 24, cy + 22), outline=color, width=3)
        draw.line((cx - 18, cy + 13, cx - 4, cy - 1, cx + 8, cy + 8, cx + 21, cy - 10), fill=color, width=3)
        draw.ellipse((cx + 8, cy - 16, cx + 17, cy - 7), outline=color, width=2)


def draw_card_glyph(draw: ImageDraw.ImageDraw, x: int, y: int, kind: str, color: str) -> tuple[int, int, int, int]:
    box = (x, y, x + 82, y + 82)
    cx, cy = x + 41, y + 41
    draw.ellipse(box, fill=tint(color, 0.88), outline=color, width=2)
    draw.ellipse((x + 8, y + 8, x + 74, y + 74), fill=tint(color, 0.94))
    if kind == "trace":
        pts = [(x + 25, y + 29), (x + 51, y + 25), (x + 58, y + 54), (x + 29, y + 58)]
        draw.line(pts + [pts[0]], fill=color, width=3)
        for px, py in pts:
            draw.ellipse((px - 7, py - 7, px + 7, py + 7), fill=WHITE, outline=color, width=2)
    elif kind == "rank":
        base_y = y + 61
        for bx, bh in [(x + 24, 28), (x + 39, 42), (x + 54, 34)]:
            draw.rounded_rectangle((bx, base_y - bh, bx + 10, base_y), radius=3, fill=color)
        draw.line((x + 21, base_y, x + 66, base_y), fill=color, width=3)
        draw.line((x + 52, y + 24, x + 61, y + 32, x + 72, y + 18), fill=color, width=3, joint="curve")
    elif kind == "temporal":
        draw.line((x + 22, y + 60, x + 68, y + 60), fill=color, width=2)
        draw.line((x + 22, y + 60, x + 22, y + 22), fill=color, width=2)
        curve = [(x + 24, y + 54), (x + 36, y + 47), (x + 46, y + 31), (x + 60, y + 36), (x + 68, y + 25)]
        draw.line(curve, fill=color, width=3, joint="curve")
        for px, py in curve[1:4]:
            draw.ellipse((px - 4, py - 4, px + 4, py + 4), fill=WHITE, outline=color, width=2)
    elif kind == "identity":
        draw.ellipse((x + 20, y + 21, x + 42, y + 43), fill=color)
        draw.ellipse((x + 48, y + 26, x + 66, y + 44), fill=tint(color, 0.25), outline=color, width=2)
        draw.ellipse((x + 34, y + 49, x + 58, y + 66), fill=tint(color, 0.45), outline=color, width=2)
        draw.line((x + 42, y + 33, x + 48, y + 35), fill=color, width=3)
        draw.line((x + 44, y + 44, x + 42, y + 50), fill=color, width=3)
    else:
        draw.line((cx, y + 19, cx, y + 63), fill=color, width=3)
        draw.line((x + 22, cy, x + 60, cy), fill=color, width=3)
        draw.line((x + 28, y + 28, x + 54, y + 54), fill=color, width=3)
        draw.line((x + 54, y + 28, x + 28, y + 54), fill=color, width=3)
        draw.ellipse((cx - 10, cy - 10, cx + 10, cy + 10), fill=WHITE, outline=color, width=3)
    return box


@dataclass(frozen=True)
class Step:
    n: int
    title: str
    body: str
    script: str
    color: str
    fill: str
    glyph: str


LEFT_STEPS = [
    Step(1, "Import and label source files", "The GEO source files are assembled into a single-cell object that connects counts to cells, genes, subjects, cohorts, and infection days.", "prepare_gse154386_subset_r.R and prepare_gse154386_subset_python.py rebuild this object.", UPSTREAM_BLUE, "#EEF1F4", "file"),
    Step(2, "Prepare the discovery matrix", "Filtering, normalization, variable-gene selection, and balanced sampling produce the matrix sent to CoGAPS.", "preprocessing_manifest.json records thresholds, sampling, and matrix orientation.", UPSTREAM_GREEN, "#EEF2EF", "filter"),
    Step(3, "Run candidate-rank models", "The same discovery matrix is fit at several ranks before any biological interpretation is attached to a pattern.", "gse154386_make_cogaps_jobs_tsv.py writes the sweep job table.", SWEEP_PURPLE, "#F0EEF3", "bars"),
    Step(4, "Compare rank evidence", "Rank summaries compare stability, separation, and interpretability so the analysis does not depend on a convenient rank alone.", "cs5_generate_revised_k_selection_report.py gathers the comparison tables.", CORE_ORANGE, "#FFF7ED", "scale"),
    Step(5, "Export model outputs", "The selected model run is saved as cell scores, gene loadings, top genes, diagnostics, trace, and metadata.", "cs5_run_selected_model_r.R and cs5_run_selected_model_python.py export model outputs.", CORE_BLUE, "#EEF7FC", "download"),
    Step(6, "Summarize scores by biology", "Cell scores are summarized by subject, infection day, and broad PBMC identity before they are read as patterns.", "cs5_build_k10_interpretation_layer.py builds the interpretation tables.", CORE_GREEN, "#ECFDF5", "pie"),
    Step(7, "Add expression direction", "CoGAPS connects genes to patterns, but a separate baseline comparison shows whether top genes increase or decrease after infection.", "gse154386_pattern_directionality.py builds the pseudobulk comparison tables.", MUTED_ROSE, "#F2EEF0", "direction"),
    Step(8, "Assemble figures and tables", "The final figures and tables bring together timing, cell context, embeddings, trajectories, and gene-direction summaries.", "Figure builders include generate_pattern_time_small_multiples.R and generate_pattern_embedding_figure.py.", CORE_GOLD, "#FFF8E5", "image"),
]


STEP_CHIPS = {
    1: ("Source labels", UPSTREAM_BLUE),
    2: ("Discovery matrix", UPSTREAM_GREEN),
    3: ("Rank sweep", SWEEP_PURPLE),
    4: ("Rank evidence", CORE_ORANGE),
    5: ("Model outputs", CORE_BLUE),
    6: ("Score summaries", CORE_GREEN),
    7: ("Directionality", MUTED_ROSE),
    8: ("Figures + tables", CORE_GOLD),
}


@dataclass(frozen=True)
class RightCard:
    label: str
    glyph: str
    title: str
    body: str
    objective: str
    steps: tuple[int, ...]
    y: int
    h: int
    color: str


RIGHT_CARDS = [
    RightCard("Trace", "trace", "Keep the evidence chain visible", "Before interpretation begins, source labels and preprocessing choices stay attached to the matrix CoGAPS used.", "Objective: trace where the included analysis files came from.", (1, 2), 310, 390, TEAL),
    RightCard("Q1", "rank", "Choose a model before naming patterns", "Rank evidence and run diagnostics decide which selected CoGAPS output is ready to inspect.", "Objective: evaluate K with stability, resolution, diagnostics, and sensitivity evidence.", (3, 4, 5), 790, 455, ORANGE),
    RightCard("Q2", "temporal", "Look for temporal programs", "Subject-timepoint summaries and trajectories ask which patterns change across infection days.", "Objective: summarize scores by subject and timepoint before interpreting timing.", (5, 6, 8), 1335, 410, GREEN),
    RightCard("Q3", "identity", "Separate identity from activity", "PBMC identity context stays beside timing evidence, so cell-type structure is not mistaken for infection dynamics.", "Objective: distinguish identity-like, activity-like, and mixed signals.", (2, 6, 8), 1840, 420, TEAL),
    RightCard("Q4", "ifn", "Weigh interferon evidence", "Timing, top genes, IFN signal, and directionality are combined before making the late acute response interpretation.", "Objective: connect pattern evidence without overclaiming.", (5, 6, 7, 8), 2350, 430, RED),
]


@dataclass(frozen=True)
class RouteCard:
    label: str
    glyph: str
    title: str
    subtitle: str
    body: str
    callout: str
    steps: tuple[int, ...]
    y: int
    h: int
    color: str
    fill: str


ROUTE_CARDS = [
    RouteCard(
        "Main",
        "temporal",
        "Main lesson path",
        "Work with included evidence files.",
        "Begin with included rank evidence and selected-model files, then follow summaries, expression direction, and final figures into the interpretation sections.",
        "Use this path to complete the case study.",
        (4, 5, 6, 7, 8),
        300,
        820,
        CORE_GREEN,
        "#ECFDF5",
    ),
    RouteCard(
        "Rebuild",
        "trace",
        "Full reproduction path",
        "Start from source files, then rejoin the main lesson.",
        "Start with GEO source files, rebuild the discovery matrix, then fit and export the selected model with a preselected rank. This path skips the candidate-rank sweep.",
        "Use this when you want to recreate the included files.",
        (1, 2, 4, 5, 6, 7, 8),
        1210,
        1060,
        UPSTREAM_BLUE,
        "#EEF1F4",
    ),
    RouteCard(
        "Sweep",
        "rank",
        "High-compute rank sweep",
        "Separate candidate-rank search.",
        "Rerun candidate-rank models to recreate the evidence behind the rank choice. This path is separate from the main lesson and from full reproduction because it needs high-performance computing.",
        "Do not run this on a typical laptop.",
        (3,),
        2370,
        520,
        SWEEP_PURPLE,
        "#F0EEF3",
    ),
]

ROUTE_SEQUENCES = {
    "Main": (
        "Main lesson sequence",
        [
            (4, "Compare included rank evidence"),
            (5, "Load selected-model outputs"),
            (6, "Summarize scores by subject, day, and identity"),
            (7, "Check expression direction"),
            (8, "Read final figures and tables"),
        ],
    ),
    "Rebuild": (
        "Full reproduction sequence",
        [
            (1, "Import and label source files"),
            (2, "Prepare the discovery matrix"),
            (4, "Use included rank evidence"),
            (5, "Fit and export preselected model"),
            (6, "Summarize scores by subject, day, and identity"),
            (7, "Add expression direction"),
            (8, "Assemble figures and tables"),
        ],
    ),
    "Sweep": (
        "High-compute sweep",
        [
            (3, "Run candidate-rank models"),
        ],
    ),
}


def draw_chip(draw: ImageDraw.ImageDraw, audit: Audit, name: str, x: int, y: int, step: int, parent: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    label, color = STEP_CHIPS[step]
    label_w, _ = text_size(draw, label, CHIP)
    chip_w = 52 + label_w + 20
    box = (x, y, x + chip_w, y + 40)
    draw.rounded_rectangle(box, radius=20, fill=tint(color, 0.9), outline=color, width=2)
    draw.ellipse((x + 7, y + 7, x + 33, y + 33), fill=color)
    n = str(step)
    nw, nh = text_size(draw, n, CHIP_BOLD)
    draw.text((x + 7 + (26 - nw) / 2, y + 7 + (26 - nh) / 2 - 1), n, font=CHIP_BOLD, fill=WHITE)
    draw.text((x + 42, y + 10), label, font=CHIP, fill=INK)
    audit.add(name, box, parent)
    return box


def draw_chips(draw: ImageDraw.ImageDraw, audit: Audit, name: str, x: int, y: int, max_width: int, steps: tuple[int, ...], parent: tuple[int, int, int, int]) -> tuple[int, tuple[int, int, int, int]]:
    boxes: list[tuple[str, tuple[int, int, int, int]]] = []
    row_x = x
    row_y = y
    for step in steps:
        label_w, _ = text_size(draw, STEP_CHIPS[step][0], CHIP)
        width = 52 + label_w + 20
        if row_x + width > x + max_width:
            row_x = x
            row_y += 50
        box = draw_chip(draw, audit, f"{name}-{step}", row_x, row_y, step, parent)
        boxes.append((f"{name}-{step}", box))
        row_x = box[2] + 10
    full = (x, y, max(b[1][2] for b in boxes), max(b[1][3] for b in boxes))
    audit.add(name, full, parent)
    audit.no_overlap(name, boxes, pad=0)
    return full[3], full


def draw_dashed_ring(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], color: str, width: int = 7) -> None:
    for start in range(0, 360, 48):
        draw.arc(box, start=start, end=start + 28, fill=color, width=width)


def draw_step_marker(draw: ImageDraw.ImageDraw, x: int, y: int, step: Step) -> None:
    if step.n in MAIN_LESSON_MARKER_STEPS:
        draw.ellipse((x - 52, y - 52, x + 52, y + 52), outline=CORE_GOLD, width=8)
        draw.ellipse((x - 45, y - 45, x + 45, y + 45), outline=WHITE, width=6)
    draw.ellipse((x - 40, y - 40, x + 40, y + 40), fill=step.color)
    n = str(step.n)
    nw, nh = text_size(draw, n, font(34, True))
    draw.text((x - nw / 2, y - nh / 2 - 2), n, font=font(34, True), fill=WHITE)


def draw_ring_key(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    box = (x, y, x + 590, y + 58)
    draw.rounded_rectangle(box, radius=18, fill=WHITE, outline="#D8E4EF", width=2)
    items = [
        ("Main evidence", GREEN, "solid"),
        ("Prepared upstream", SLATE, "dashed"),
        ("Used later", RED, "dashed"),
    ]
    legend_font = font(17)
    cursor = x + 30
    for label, color, style in items:
        label_w, label_h = text_size(draw, label, legend_font)
        icon_center_y = y + 29
        icon_box = (cursor, icon_center_y - 11, cursor + 22, icon_center_y + 11)
        if style == "solid":
            draw.ellipse(icon_box, outline=color, width=5)
        else:
            draw_dashed_ring(draw, icon_box, color, width=4)
        draw.text((cursor + 38, y + (58 - label_h) / 2 - 3), label, font=legend_font, fill=INK)
        cursor += 22 + 38 + label_w + 24


def draw_setup_card(image: Image.Image, audit: Audit, box: tuple[int, int, int, int], title: str, sub: str, body: str, color: str, kind: str) -> None:
    draw = ImageDraw.Draw(image, "RGBA")
    shadow_rect(image, box, 16, WHITE, color, width=3, shadow_alpha=9, shadow_blur=10, shadow_offset=(0, 5))
    draw = ImageDraw.Draw(image, "RGBA")
    icon_box = (box[0] + 32, box[1] + 42, box[0] + 108, box[1] + 118)
    draw_setup_icon(draw, icon_box[0], icon_box[1], kind, color)
    text_x = box[0] + 138
    inner = (box[0] + 24, box[1] + 24, box[2] - 24, box[3] - 24)
    _, tbox = draw_text_block(draw, audit, f"{title}-setup-title", text_x, box[1] + 42, title, H3, box[2] - text_x - 24, color, inner, 6)
    _, sbox = draw_text_block(draw, audit, f"{title}-setup-sub", text_x, tbox[3] + 12, sub, SMALL_BOLD, box[2] - text_x - 24, color, inner, 4)
    draw_text_block(draw, audit, f"{title}-setup-body", box[0] + 32, box[1] + 134, body, BODY, box[2] - box[0] - 64, INK, inner, 7)
    audit.add(f"{title}-setup-card", box)


def draw_step_card(image: Image.Image, audit: Audit, step: Step, x: int, y: int, w: int, h: int) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(image, "RGBA")
    box = (x, y, x + w, y + h)
    inner = (x + 32, y + 24, x + w - 32, y + h - 24)
    shadow_rect(image, box, 18, step.fill, step.color, width=3, shadow_alpha=16, shadow_blur=12, shadow_offset=(0, 7))
    draw = ImageDraw.Draw(image, "RGBA")
    draw.rounded_rectangle((x, y, x + 18, y + h), radius=18, fill=step.color)
    draw.rectangle((x + 9, y, x + 26, y + h), fill=step.color)
    glyph_box = (x + w - 126, y + 48, x + w - 44, y + 130)
    draw_step_glyph(draw, glyph_box, step.glyph, step.color)
    title_x = x + 54
    title_w = w - 210
    cursor, title_box = draw_text_block(draw, audit, f"step-{step.n}-title", title_x, y + 32, step.title, H3, title_w, step.color, inner, 7)
    cursor += 12
    cursor, body_box = draw_text_block(draw, audit, f"step-{step.n}-body", title_x, cursor, step.body, BODY, title_w, INK, inner, 8)
    script_y = y + h - 74
    script_box = (title_x, script_y, x + w - 34, y + h - 24)
    draw.rounded_rectangle(script_box, radius=12, fill=tint(step.color, 0.92))
    audit.add(f"step-{step.n}-script-box", script_box, inner)
    draw_text_block(draw, audit, f"step-{step.n}-script", title_x + 12, script_y + 11, step.script, SCRIPT, script_box[2] - script_box[0] - 24, step.color, script_box, 5)
    audit.no_overlap(
        f"step-{step.n}-major",
        [("title", title_box), ("body", body_box), ("script", script_box), ("glyph", glyph_box)],
        pad=6,
    )
    return box


def draw_right_card(image: Image.Image, audit: Audit, card: RightCard, x: int, w: int) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(image, "RGBA")
    box = (x, card.y, x + w, card.y + card.h)
    inner = (x + 28, card.y + 26, x + w - 28, card.y + card.h - 26)
    shadow_rect(image, box, 18, WHITE, card.color, width=3, shadow_alpha=13, shadow_blur=12, shadow_offset=(0, 7))
    draw = ImageDraw.Draw(image, "RGBA")
    draw.rounded_rectangle((x, card.y, x + 18, card.y + card.h), radius=18, fill=card.color)
    draw.rectangle((x + 9, card.y, x + 26, card.y + card.h), fill=card.color)
    glyph_box = draw_card_glyph(draw, x + 42, card.y + 34, card.glyph, card.color)
    audit.add(f"{card.label}-glyph", glyph_box, inner)
    label_w, _ = text_size(draw, card.label, SMALL_BOLD)
    label_box = (x + 43, card.y + 128, x + 43 + label_w + 34, card.y + 164)
    draw.rounded_rectangle(label_box, radius=18, fill=card.color)
    draw.text((label_box[0] + 17, label_box[1] + 7), card.label, font=SMALL_BOLD, fill=WHITE)
    audit.add(f"{card.label}-label", label_box, inner)
    tx = x + 158
    tw = w - 196
    cursor, title_box = draw_text_block(draw, audit, f"{card.label}-title", tx, card.y + 32, card.title, H3, tw, INK, inner, 6)
    cursor += 12
    cursor, body_box = draw_text_block(draw, audit, f"{card.label}-body", tx, cursor, card.body, BODY, tw, INK, inner, 8)
    cursor += 18
    objective_box = (tx, cursor, tx + tw, cursor + 62)
    draw.rounded_rectangle(objective_box, radius=14, fill=tint(card.color, 0.92))
    audit.add(f"{card.label}-objective-box", objective_box, inner)
    draw_text_block(draw, audit, f"{card.label}-objective-text", tx + 16, cursor + 12, card.objective, SMALL_BOLD, tw - 32, card.color, objective_box, 5)
    cursor += 80
    label_text = "Evidence labels"
    draw.text((tx, cursor), label_text, font=TINY, fill=MUTED)
    chip_label_box = (tx, cursor, tx + text_size(draw, label_text, TINY)[0], cursor + text_size(draw, label_text, TINY)[1])
    audit.add(f"{card.label}-chip-label", chip_label_box, inner)
    cursor += 25
    chip_bottom, chips_box = draw_chips(draw, audit, f"{card.label}-chips", tx, cursor, tw, card.steps, inner)
    audit.no_overlap(
        f"{card.label}-major",
        [("glyph", glyph_box), ("label", label_box), ("title", title_box), ("body", body_box), ("objective", objective_box), ("chips", chips_box)],
        pad=8,
    )
    if chip_bottom > inner[3]:
        audit.failures.append({"type": "card-overflow", "card": card.label, "bottom": chip_bottom, "limit": inner[3]})
    return box


def draw_route_card(image: Image.Image, audit: Audit, card: RouteCard, x: int, w: int) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(image, "RGBA")
    box = (x, card.y, x + w, card.y + card.h)
    inner = (x + 30, card.y + 28, x + w - 30, card.y + card.h - 28)
    shadow_rect(image, box, 18, card.fill, card.color, width=3, shadow_alpha=13, shadow_blur=12, shadow_offset=(0, 7))
    draw = ImageDraw.Draw(image, "RGBA")
    draw.rounded_rectangle((x, card.y, x + 18, card.y + card.h), radius=18, fill=card.color)
    draw.rectangle((x + 9, card.y, x + 28, card.y + card.h), fill=card.color)

    glyph_box = draw_card_glyph(draw, x + 44, card.y + 38, card.glyph, card.color)
    audit.add(f"{card.label}-route-glyph", glyph_box, inner)
    label_w, _ = text_size(draw, card.label, SMALL_BOLD)
    label_box = (x + 44, card.y + 134, x + 44 + label_w + 34, card.y + 170)
    draw.rounded_rectangle(label_box, radius=18, fill=card.color)
    draw.text((label_box[0] + 17, label_box[1] + 7), card.label, font=SMALL_BOLD, fill=WHITE)
    audit.add(f"{card.label}-route-label", label_box, inner)

    tx = x + 164
    tw = w - 206
    cursor, title_box = draw_text_block(draw, audit, f"{card.label}-route-title", tx, card.y + 34, card.title, H3, tw, INK, inner, 6)
    cursor += 10
    cursor, subtitle_box = draw_text_block(draw, audit, f"{card.label}-route-subtitle", tx, cursor, card.subtitle, BODY_BOLD, tw, card.color, inner, 6)
    cursor += 12
    cursor, body_box = draw_text_block(draw, audit, f"{card.label}-route-body", tx, cursor, card.body, BODY, tw, INK, inner, 8)
    cursor += 20

    callout_box = (tx, cursor, tx + tw, cursor + 66)
    draw.rounded_rectangle(callout_box, radius=14, fill=tint(card.color, 0.92))
    audit.add(f"{card.label}-route-callout-box", callout_box, inner)
    draw_text_block(draw, audit, f"{card.label}-route-callout", tx + 16, cursor + 12, card.callout, SMALL_BOLD, tw - 32, card.color, callout_box, 5)
    cursor += 84

    sequence_boxes: list[tuple[str, tuple[int, int, int, int]]] = []
    sequence_label, sequence_rows = ROUTE_SEQUENCES[card.label]
    label_w, label_h = text_size(draw, sequence_label, TINY)
    draw.text((tx, cursor), sequence_label, font=TINY, fill=MUTED)
    sequence_boxes.append(("sequence-label", (tx, cursor, tx + label_w, cursor + label_h)))
    audit.add(f"{card.label}-route-sequence-label", sequence_boxes[-1][1], inner)
    cursor += 28
    for step_number, label in sequence_rows:
        _, color = STEP_CHIPS[step_number]
        text_x = tx + 54
        text_max = tw - 72
        row_lines = wrap(draw, label, SMALL_BOLD, text_max)[:2]
        row_h = 48 if len(row_lines) == 1 else 66
        row_box = (tx, cursor, tx + tw, cursor + row_h)
        draw.rounded_rectangle(row_box, radius=14, fill=tint(color, 0.93), outline=tint(color, 0.65), width=1)
        number_y = cursor + (row_h - 28) / 2
        draw.ellipse((tx + 12, number_y, tx + 40, number_y + 28), fill=color)
        n = str(step_number)
        nw, nh = text_size(draw, n, CHIP_BOLD)
        draw.text((tx + 12 + (28 - nw) / 2, number_y + (28 - nh) / 2 - 1), n, font=CHIP_BOLD, fill=WHITE)
        total_text_h = sum(text_size(draw, line, SMALL_BOLD)[1] for line in row_lines) + 4 * (len(row_lines) - 1)
        line_y = cursor + (row_h - total_text_h) / 2 - 1
        for idx, line in enumerate(row_lines):
            draw.text((text_x, line_y), line, font=SMALL_BOLD, fill=INK)
            lw, lh = text_size(draw, line, SMALL_BOLD)
            audit.add(f"{card.label}-route-sequence-{step_number}-line-{idx}", (text_x, line_y, text_x + lw, line_y + lh), row_box)
            line_y += lh + 4
        sequence_boxes.append((f"sequence-{step_number}", row_box))
        audit.add(f"{card.label}-route-sequence-{step_number}", row_box, inner)
        cursor += row_h + 10
    audit.no_overlap(
        f"{card.label}-route-major",
        [
            ("glyph", glyph_box),
            ("label", label_box),
            ("title", title_box),
            ("subtitle", subtitle_box),
            ("body", body_box),
            ("callout", callout_box),
        ]
        + sequence_boxes,
        pad=8,
    )
    route_bottom = max(box[3] for _, box in sequence_boxes)
    if route_bottom > inner[3]:
        audit.failures.append({"type": "route-card-overflow", "card": card.label, "bottom": route_bottom, "limit": inner[3]})
    return box


def step_center_y(step_number: int) -> int:
    step_h = 310
    step_gap = 35
    return 245 + (step_number - 1) * (step_h + step_gap) + step_h // 2


def draw_route_connectors(draw: ImageDraw.ImageDraw) -> None:
    x0 = 895
    x1 = 1040
    for card in ROUTE_CARDS:
        for step_number in card.steps:
            y0 = step_center_y(step_number)
            y1 = ROUTE_DESTINATIONS[card.label][step_number]
            color = STEP_CHIPS[step_number][1]
            draw.line((x0, y0, x1, y1), fill=rgba(color, 135), width=6)
            draw.ellipse((x0 - 7, y0 - 7, x0 + 7, y0 + 7), fill=color)
            draw.ellipse((x1 - 7, y1 - 7, x1 + 7, y1 + 7), fill=color)


def save_with_qa(image: Image.Image, audit: Audit, png_path: Path, rendering: str) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    image.convert("RGB").save(png_path, quality=95)
    qa = {
        "output": str(png_path),
        "size": list(image.size),
        "rendering": rendering,
        "audit_item_count": len(audit.items),
        "audit_failures": audit.failures,
        "status": "pass" if not audit.failures else "review",
    }
    print(json.dumps(qa, indent=2))


def generate_resource_figure() -> None:
    image = Image.new("RGBA", (RESOURCE_W, RESOURCE_H), rgb(BG) + (255,))
    audit = Audit()
    draw = ImageDraw.Draw(image, "RGBA")

    shadow_rect(image, SETUP, 24, WHITE, PANEL_LINE, width=2, shadow_alpha=10, shadow_blur=14, shadow_offset=(0, 8))
    draw = ImageDraw.Draw(image, "RGBA")
    draw.text((86, 104), "What setup provides", font=H2, fill=NAVY)
    draw.text((535, 115), "Docker supplies software; the repository supplies included evidence; the Full Reproduction Guide covers rebuild and sweep details.", font=SMALL, fill=INK)
    draw_setup_card(image, audit, (95, 190, 600, 468), "Docker image", "software environment", "R, Python, CoGAPS, PyCoGAPS, and plotting packages are already installed.", BLUE, "docker")
    draw_setup_card(image, audit, (655, 190, 1160, 468), "GitHub repository", "case-study files", "Pages, scripts, manifests, figures, and included evidence tables are in the project folder.", GREEN, "repo")
    draw_setup_card(image, audit, (1215, 190, 1720, 468), "External archive", "rebuild and sweep-support files", "Larger files support source rebuilding, selected-model reruns, summaries, directionality, and rank-sweep checks.", ORANGE, "archive")
    draw.text((86, 516), "Most sections use Docker plus the repository. The separate high-compute rank sweep is described later.", font=SMALL, fill=MUTED)
    audit.add("setup-panel", SETUP)
    save_with_qa(image, audit, RESOURCE_PNG, "separate resource figure, full redraw")


def generate_evidence_figure() -> None:
    image = Image.new("RGBA", (EVIDENCE_W, EVIDENCE_H), rgb(BG) + (255,))
    audit = Audit()
    draw = ImageDraw.Draw(image, "RGBA")

    shadow_rect(image, LEFT_PANEL, 24, "#FAFCFF", PANEL_LINE, width=2, shadow_alpha=10, shadow_blur=14, shadow_offset=(0, 8))
    shadow_rect(image, RIGHT_PANEL, 24, "#FBFEFD", PANEL_LINE, width=2, shadow_alpha=10, shadow_blur=14, shadow_offset=(0, 8))
    draw = ImageDraw.Draw(image, "RGBA")
    draw.line((DIVIDER_X, LEFT_PANEL[1] + 10, DIVIDER_X, LEFT_PANEL[3] - 10), fill=PANEL_LINE, width=5)

    draw_panel_icon(draw, 82, 84, "network", BLUE)
    draw.text((155, 85), "Evidence-building path", font=H2, fill=NAVY)
    draw_panel_icon(draw, 1025, 84, "target", GREEN)
    draw.text((1094, 80), "Choose the path you need", font=H2, fill=NAVY)
    draw.text((1094, 128), "main lesson, full reproduction, or high-compute sweep", font=SUBTITLE, fill=MUTED)
    draw_text_block(draw, audit, "right-intro", 1036, 208, "Use these paths to see what you do in the main lesson, what you can rebuild from source files, and what requires high-performance computing.", SMALL, 680, MUTED, RIGHT_PANEL, 6)

    line_x = 112
    draw.line((line_x, 250, line_x, 2910), fill="#9DB0C5", width=3)
    step_x = 180
    step_w = 700
    step_h = 310
    step_gap = 35
    step_boxes = []
    for idx, step in enumerate(LEFT_STEPS):
        y = 245 + idx * (step_h + step_gap)
        cy = y + 62
        draw_step_marker(draw, line_x, cy, step)
        step_boxes.append((f"step-{step.n}", draw_step_card(image, audit, step, step_x, y, step_w, step_h)))
    audit.no_overlap("left-step-cards", step_boxes, pad=10)

    card_boxes = []
    for card in ROUTE_CARDS:
        card_boxes.append((card.label, draw_route_card(image, audit, card, 1040, 688)))
    audit.no_overlap("route-cards", card_boxes, pad=10)

    shadow_rect(image, FOOTER, 14, WHITE, "#D8E4EF", width=2, shadow_alpha=8, shadow_blur=12, shadow_offset=(0, 6))
    draw = ImageDraw.Draw(image, "RGBA")
    draw.ellipse((74, 3108, 102, 3136), outline=BLUE, width=3)
    draw.text((84, 3110), "i", font=font(20, True), fill=BLUE)
    draw.text((118, 3106), "How to read this map:", font=SMALL_BOLD, fill=NAVY)
    draw.text((338, 3106), "Read the paths from top to bottom: main lesson, full reproduction without the sweep, then the separate high-compute rank sweep.", font=SMALL, fill=MUTED)

    audit.add("left-panel", LEFT_PANEL)
    audit.add("right-panel", RIGHT_PANEL)
    audit.add("footer", FOOTER)
    save_with_qa(image, audit, EVIDENCE_PNG, "separate evidence-path figure, full redraw")


def main() -> None:
    generate_resource_figure()
    generate_evidence_figure()


if __name__ == "__main__":
    main()
