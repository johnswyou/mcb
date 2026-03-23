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
    """Plot simultaneous confidence intervals horizontally.

    Parameters
    ----------
    result:
        Result returned by a constrained or unconstrained MCB procedure.
    ax:
        Existing Matplotlib axes to draw on. A new figure and axes are created
        when omitted.
    title:
        Optional plot title. A default title is derived from the result when
        omitted.
    figsize:
        Figure size used only when a new figure is created.

    Returns
    -------
    matplotlib.axes.Axes
        Axes containing the rendered interval plot.
    """
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
