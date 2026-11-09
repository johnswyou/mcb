from __future__ import annotations

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mcb import MCBResult
from mcb.plotting import plot_simultaneous_intervals


def test_plot_simultaneous_intervals_returns_axes() -> None:
    result = MCBResult(
        method="Example method",
        parameter="theta",
        alpha=0.05,
        labels=("A", "B"),
        estimates=(0.0, -0.5),
        lower=(-0.1, -0.8),
        upper=(0.2, 0.0),
        critical_values={"A": 2.0, "B": 2.0},
        candidate_best_set=("A",),
    )

    ax = plot_simultaneous_intervals(result, title="Example plot")

    assert ax.get_title() == "Example plot"
    assert ax.get_xlabel() == "theta"
    assert ax.get_ylabel() == "Arm"
    # Markers sit at the estimates, not at the interval midpoints (0.05, -0.4).
    assert list(ax.containers[0][0].get_xdata()) == [0.0, -0.5]
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["A", "B"]
    plt.close(ax.figure)
