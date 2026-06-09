"""
plot_style.py
-------------
Master matplotlib/plotly style loader.
Usage:
    from plot_style import *
    from plot_style import get_n_colors, MidFactor
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import ScalarFormatter, LogLocator
from matplotlib.colors import to_rgba
from matplotlib.gridspec import GridSpec
import matplotlib.patheffects as path_effects
import plotly.graph_objects as go


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_n_colors(n, cmap_name="viridis"):
    """Sample n evenly-spaced colours from a matplotlib colormap."""
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


def MidFactor(n):
    """Return the middle divisor of n (useful for subplot grid layout)."""
    divisors = [i for i in range(1, n + 1) if n % i == 0]
    return divisors[len(divisors) // 2]


# ── RC Params ─────────────────────────────────────────────────────────────────

def load_style(usetex=True, fontsize=24, tick_labelsize=20):
    """Apply master rc params. Call once at the top of a script."""

    params = {
        # Font
        "font.family":          "serif",
        "font.size":            fontsize,
        "text.usetex":          usetex,

        # Tick labels
        "xtick.labelsize":      tick_labelsize,
        "ytick.labelsize":      tick_labelsize,

        # Tick direction & size
        "xtick.direction":      "in",
        "ytick.direction":      "in",
        "xtick.major.size":     7,
        "ytick.major.size":     7,
        "xtick.minor.size":     3.5,
        "ytick.minor.size":     3.5,
        "xtick.minor.visible":  False,
        "ytick.minor.visible":  False,

        # Tick padding & axes
        "xtick.major.pad":      7,
        "axes.linewidth":       2.0,

        # tight layout
        "figure.constrained_layout.use": True
    }

    plt.rcParams.update(params)

# Apply immediately on import
load_style()


# ── Re-exports (so callers can just `from plot_style import *`) ───────────────

__all__ = [
    "plt", "cm", "go",
    "ScalarFormatter", "LogLocator",
    "to_rgba", "GridSpec", "path_effects",
    "get_n_colors", "MidFactor", "load_style",
]
