"""Shared publication-figure infrastructure for the PHA-QNM paper."""

from .style import (COLORS, DOUBLE_COLUMN_WIDTH, SINGLE_COLUMN_WIDTH,
                    WIDE_FIGURE_WIDTH, add_panel_label, save_figure,
                    use_publication_style)

__all__ = [
    "COLORS", "DOUBLE_COLUMN_WIDTH", "SINGLE_COLUMN_WIDTH",
    "WIDE_FIGURE_WIDTH", "add_panel_label", "save_figure",
    "use_publication_style",
]
