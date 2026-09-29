#!/usr/bin/env python
"""Shared figure style for the genes/ layer — the GrENE-Net paper look.

Matched to the CAM5 panel figure the user set as the reference. The recipe, in four
rules:

  1. WHITE panel, LIGHT GREY GRID, NO FRAME. The grid carries the structure, so the box
     comes off entirely. Grid sits under the data.
  2. GREY AXES. Tick marks and tick labels in grey, axis titles only slightly darker.
     Nothing on the page is pure black.
  3. NO BOLD ANYWHERE. Panel letters are separated from the text by size and colour, not
     weight.
  4. DESATURATED PALETTE. Blue->red diverging for anything climate-ordered (cold to warm),
     a green family for statistical series, grey for context/background points, and a
     single warm accent reserved for "this is the variant we are talking about".

Nothing here sets a title: the house rule is no chart or subplot titles (see CLAUDE.md).
Panel identity is an in-panel corner annotation, drawn by `panel()` below.

Usage:
    import plot_theme as TH
    TH.apply()
    ...
    TH.panel(ax, "A", "climate association")
"""
from __future__ import annotations
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# -- ink ------------------------------------------------------------------------------------
GRID = "#E6E6E6"        # background grid
AXIS = "#C2C2C2"        # tick marks / any spine left on
TICK = "#5E5E5E"        # tick labels
LABEL = "#333333"       # axis titles
MUTED = "#8A8A8A"       # secondary annotation
FAINT = "#BFBFBF"       # context points

# -- palettes -------------------------------------------------------------------------------
# climate-ordered, cold -> warm (the reference's Temp. scale)
COLD, WARM = "#2166AC", "#B2182B"
TEMP = LinearSegmentedColormap.from_list("temp", ["#2166AC", "#89B3D4", "#F4E7E3",
                                                  "#D98C7A", "#B2182B"])
# statistical series, as the reference's LFMM / Binomial GLM / Kendall greens
GREENS = ["#8CC98C", "#3F9B49", "#1B5E20"]
# variant classes: SNPs are context, the non-SNP layer is the subject
CLASS = {"snp": FAINT, "smallindel": "#3F9B49", "sv": "#1B5E20"}
# the one warm accent, reserved for the focal variant
ACCENT = "#B2182B"
# a second, cooler accent for the SNP the eQTL is actually read at -- it has to be
# distinguishable from the focal variant (red), the variant classes (greens) and the
# tissue rows, so it is the only purple in the figure
TAG = "#6A4C93"
# tissues (ATAC): muted qualitative, no primary colours
TISSUE = {"flower": "#B48EAD", "leaf": "#8FBF7F", "root": "#C0A080", "shoot": "#8AAFCB"}
# two-way comparisons (candidates vs background)
FOCUS, BACKDROP = "#3F9B49", "#D4D4D4"


def apply(base: float = 9.5) -> None:
    """Install the style globally. Call once, before creating any figure."""
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",

        # no frame -- the grid is the structure
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": False, "axes.spines.bottom": False,

        "axes.grid": True, "axes.axisbelow": True,
        "grid.color": GRID, "grid.linewidth": 0.7, "grid.alpha": 1.0,

        "xtick.color": AXIS, "ytick.color": AXIS,
        "xtick.labelcolor": TICK, "ytick.labelcolor": TICK,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "xtick.direction": "out", "ytick.direction": "out",

        "axes.labelcolor": LABEL, "text.color": LABEL,
        "axes.labelsize": base + 0.5, "xtick.labelsize": base - 0.5,
        "ytick.labelsize": base - 0.5, "font.size": base,
        "axes.labelweight": "normal", "font.weight": "normal",

        "legend.frameon": False, "legend.fontsize": base - 1.5,
        "lines.solid_capstyle": "round",
    })


def panel(ax, letter: str, text: str = "", y: float = 0.98, x: float = 0.0) -> None:
    """Corner panel label: letter then description, distinguished by size not weight.

    The description is offset from the letter in POINTS, not axes fraction -- a fraction
    scales with panel width, so the same value that looks right on a tall narrow panel
    collides with the letter on a wide short one.
    """
    ax.annotate(letter, (x, y), xycoords="axes fraction", va="top", ha="left",
                fontsize=mpl.rcParams["font.size"] + 3.5, color=LABEL)
    if text:
        ax.annotate(text, (x, y), xycoords="axes fraction",
                    xytext=(mpl.rcParams["font.size"] + 6, 0), textcoords="offset points",
                    va="top", ha="left", fontsize=mpl.rcParams["font.size"], color=MUTED)


def note(ax, text: str, y: float = 0.86, x: float = 0.0, color: str | None = None,
         ha: str = "left") -> None:
    """A short in-panel note, in the secondary ink unless a colour is given."""
    ax.annotate(text, (x, y), xycoords="axes fraction", va="top", ha=ha,
                fontsize=mpl.rcParams["font.size"] - 1.5, color=color or MUTED)


def grid_only(ax, axis: str = "both") -> None:
    """Restrict the grid to one direction (e.g. horizontal rules on a categorical panel)."""
    ax.grid(False)
    ax.grid(True, axis=axis, color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
