"""Optional plotting helpers for MCB results."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .core import MCBResult

__all__ = ["plot_simultaneous_intervals"]


def plot_simultaneous_intervals(
    result: MCBResult,
    *,
    ax: plt.Axes | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (8.5, 4.8),
) -> plt.Axes:
    """Plot simultaneous confidence intervals horizontally."""
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    labels = list(result.labels)
    lower = np.asarray(result.lower, dtype=float)
    upper = np.asarray(result.upper, dtype=float)
    centers = (lower + upper) / 2.0
    xerr = np.vstack([centers - lower, upper - centers])
    y = np.arange(len(labels))[::-1]

    ax.errorbar(centers, y, xerr=xerr, fmt="o", capsize=4)
    ax.axvline(0.0, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(result.parameter)
    ax.set_ylabel("Variant")
    if title is None:
        title = f"{result.method} simultaneous {100 * (1 - result.alpha):.1f}% CIs"
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.3)
    return ax
