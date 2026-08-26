"""JHEP-sized, color-accessible plotting style shared by every figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl


SINGLE_COLUMN_WIDTH = 3.35
DOUBLE_COLUMN_WIDTH = 7.10
WIDE_FIGURE_WIDTH = 7.20

COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#111111",
    "gray": "#6B7280",
    "light_gray": "#D7DCE0",
}

METADATA = {
    "Author": "PHA QNM collaboration",
    "Creator": "Matplotlib",
}


def use_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8.3,
        "axes.labelsize": 8.8,
        "axes.titlesize": 8.8,
        "legend.fontsize": 7.3,
        "xtick.labelsize": 7.7,
        "ytick.labelsize": 7.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.75,
        "lines.linewidth": 1.35,
        "lines.markersize": 3.4,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.transparent": False,
    })


def add_panel_label(ax, label: str, *, x: float = 0.02, y: float = 0.98) -> None:
    ax.text(x, y, label, transform=ax.transAxes, ha="left", va="top",
            fontweight="semibold", zorder=20)


def save_figure(fig, base_path: Path, *, title: str, subject: str,
                dpi: int = 320) -> None:
    """Write matching vector PDF and inspection PNG with stable metadata."""
    base_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {**METADATA, "Title": title, "Subject": subject}
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight",
                pad_inches=0.035, metadata=metadata)
    fig.savefig(base_path.with_suffix(".png"), dpi=dpi, bbox_inches="tight",
                pad_inches=0.035)

