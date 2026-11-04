#!/usr/bin/env python3
"""Draw the figures in docs/sample_size.md from simulated tests analyzed with mcb.

The SVGs are committed, so building the docs does not run this script.
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from mcb import constrained_mcb, summarize_from_means_ns_and_mse

DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parents[1] / "docs" / "figures"

K, SIGMA, DELTA_STAR, ALPHA, PLANNED_N = 5, 40.0, 2.0, 0.05, 3797
NU = K * (PLANNED_N - 1)
ARMS = ("Control", "A", "B", "C", "D")
TRUE_MEANS = np.array([50.0, 50.0, 51.0, 53.5, 54.5])
BEST = TRUE_MEANS.max()
INNER, OUTER = BEST - DELTA_STAR, BEST - 2 * DELTA_STAR

INK, INK_2, MUTED, AXIS = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
RED, C_BLUE, D_BLUE = "#e34948", "#6da7ec", "#184f95"
RED_TINT, GRAY_TINT, BLUE_TINT = "#fbe0df", "#f0efec", "#d6e6fa"
BAR, HAIRLINE = "#ecebe5", "#d6d4cd"

STYLE = {
    # Carlito has Lato's metrics, the docs theme's body font, plus Greek. DejaVu is the fallback.
    "font.family": ["Carlito", "DejaVu Sans"],
    "font.size": 12.5,
    "svg.fonttype": "path",
    "svg.hashsalt": "mcb",  # keeps regenerated SVGs byte-stable
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.spines.left": False,
    "xtick.color": INK_2,
    "ytick.color": INK,
    "xtick.major.size": 0,
    "ytick.major.size": 0,
}


class SimulatedTest(NamedTuple):
    means: np.ndarray
    w: float
    leader: float
    cut: float
    winner: str
    ruled_out: set[str]


def analyze(means: np.ndarray, n: int, sigma_hat: float) -> SimulatedTest:
    summary = summarize_from_means_ns_and_mse(means, [n] * K, sigma_hat**2, labels=ARMS)
    result = constrained_mcb(summary, alpha=ALPHA)
    w = result.critical_values["D"] * sigma_hat * np.sqrt(2 / n)  # balanced, so every arm has the same d
    leader = float(means.max())
    cut = leader - w
    ruled_out = set(ARMS) - set(result.candidate_best_set)
    assert ruled_out == {arm for arm, mean in zip(ARMS, means) if mean <= cut}, "the cut line must match constrained_mcb"
    return SimulatedTest(means, w, leader, cut, ARMS[int(means.argmax())], ruled_out)


def followed_test(n: int) -> SimulatedTest:
    # Seed chosen so the test shows the guide's story: the true best loses at n = 3797,
    # a clearly worse arm survives at n = 1000, and the gray arm goes at n = 15,000.
    rng = np.random.default_rng(35972)
    errors = rng.normal(size=K)
    sigma_hat = SIGMA * np.sqrt(rng.chisquare(NU) / NU)
    return analyze(TRUE_MEANS + errors * SIGMA / np.sqrt(n), n, sigma_hat)


def consecutive_tests(reps: int = 16) -> list[SimulatedTest]:
    rng = np.random.default_rng(1)
    tests = []
    for _ in range(reps):
        means = TRUE_MEANS + rng.normal(0, SIGMA / np.sqrt(PLANNED_N), K)
        tests.append(analyze(means, PLANNED_N, SIGMA * np.sqrt(rng.chisquare(NU) / NU)))
    return sorted(tests, key=lambda test: test.leader)


def draw_zones(ax: plt.Axes, left: float) -> None:
    for start, end, color in ((left, OUTER, RED_TINT), (OUTER, INNER, GRAY_TINT), (INNER, BEST, BLUE_TINT)):
        ax.axvspan(start, end, color=color, lw=0, zorder=0)
    for x in (OUTER, INNER):
        ax.axvline(x, color="white", lw=2, zorder=0.5)


def draw_arms(ax: plt.Axes, test: SimulatedTest, size: float) -> None:
    for y, (arm, true_mean, mean) in enumerate(zip(ARMS, TRUE_MEANS, test.means)):
        ax.plot([true_mean, mean], [y, y], color=INK_2, lw=1, zorder=2)
        ax.plot(true_mean, y, "o", ms=size, mfc="white", mec=INK_2, mew=1.4, zorder=3)
        ax.plot(mean, y, "o", ms=size, mfc=MUTED if arm in test.ruled_out else INK, mec="white", mew=1.2, zorder=4)


def draw_lines(ax: plt.Axes, test: SimulatedTest, top: float) -> None:
    ax.plot([test.leader] * 2, [-0.45, top], color=INK_2, lw=1, zorder=1.5)
    ax.plot([test.cut] * 2, [-0.45, top], color=INK, lw=2, zorder=1.5)


def draw_bracket(ax: plt.Axes, start: float, end: float, y: float) -> None:
    style = dict(arrowstyle="|-|,widthA=.35,widthB=.35", color=INK, lw=1.6, shrinkA=0, shrinkB=0)
    ax.annotate("", (start, y), (end, y), arrowprops=style)


def set_minutes_axis(ax: plt.Axes) -> None:
    ax.set_xticks([OUTER, INNER, BEST])
    ax.set_xticklabels([f"best − 2δ*\n{OUTER}", f"best − δ*\n{INNER}", f"true best\n{BEST}"])
    ax.set_xlabel("Minutes watched per user", color=INK_2, labelpad=8)


def one_test_figure() -> plt.Figure:
    test = followed_test(PLANNED_N)
    assert test.winner == "C" and test.ruled_out == {"Control", "A"}, "the guide narrates this outcome"

    fig, ax = plt.subplots(figsize=(9.2, 4.9))
    left, right = 48.6, 56.5
    draw_zones(ax, left)
    draw_arms(ax, test, 10)
    draw_lines(ax, test, 4.62)
    draw_bracket(ax, test.cut, test.leader, 4.62)
    ax.text((test.cut + test.leader) / 2, 4.75, f"w = {test.w:.2f}", ha="center", va="bottom", color=INK)
    zones = (
        (49.55, "clearly worse", "more than 2δ* below best"),
        (51.5, "gray zone", "δ* to 2δ* below best"),
        (53.5, "good enough", "within δ* of best"),
    )
    for x, name, meaning in zones:
        ax.text(x, 5.62, name, ha="center", va="center", color=INK, fontweight="bold")
        ax.text(x, 5.28, meaning, ha="center", va="center", color=INK_2, fontsize=11)
    for y, arm in enumerate(ARMS):
        status = "observed winner" if arm == test.winner else "ruled out" if arm in test.ruled_out else "kept"
        ax.text(
            54.82,
            y,
            status,
            va="center",
            color=INK_2 if arm in test.ruled_out else INK,
            fontweight="bold" if arm == test.winner else "normal",
        )
    ax.text(test.cut - 0.08, -0.42, "ruled out", ha="right", va="center", color=INK, fontsize=11.5)
    ax.text(test.cut + 0.08, -0.42, "kept", ha="left", va="center", color=INK, fontsize=11.5)
    ax.text(test.leader + 0.08, -0.42, "observed leader", ha="left", va="center", color=INK_2, fontsize=11.5)
    ax.set_xlim(left, right)
    ax.set_ylim(-0.7, 5.9)
    ax.set_yticks(range(K))
    ax.set_yticklabels(ARMS)
    set_minutes_axis(ax)
    ax.plot([], [], "o", ms=10, mfc="white", mec=INK_2, mew=1.4, ls="", label="True mean (hidden from you)")
    ax.plot([], [], "o", ms=10, mfc=INK, mec="white", ls="", label="Sample mean (what the test shows)")
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=2,
        frameon=False,
        fontsize=12,
        handletextpad=0.3,
        columnspacing=1.6,
    )
    fig.tight_layout()
    return fig


def more_users_figure() -> plt.Figure:
    notes = {1000: "wider than δ*", PLANNED_N: "just inside δ*", 15000: "well inside δ*"}
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 4.4), sharey=True)
    left, right = 47.6, 55.2
    for ax, (n, note) in zip(axes, notes.items()):
        test = followed_test(n)
        draw_zones(ax, left)
        draw_arms(ax, test, 7.5)
        draw_lines(ax, test, 4.45)
        # Pink and gray nearly merge under red-green colorblindness, so name every zone.
        for x, name in ((49.05, "clearly worse"), (51.5, "gray zone"), (53.5, "good enough")):
            ax.text(x, 4.86, name, ha="center", va="center", fontsize=9.5, color=INK)
        ax.set_xlim(left, right)
        ax.set_ylim(-1.75, 5.2)
        label = str(n) if n < 10000 else f"{n:,}"
        ax.set_title(f"n = {label} per arm\nw = {test.w:.2f}, {note}", fontsize=12.5, color=INK, loc="left", pad=8)
        ax.set_xticks([OUTER, INNER, BEST])
        ax.set_xticklabels(["best−2δ*", "best−δ*", "best"], fontsize=11)
        ax.set_yticks(range(K))
        ax.set_yticklabels(ARMS)
        ax.barh(-0.92, DELTA_STAR, left=test.leader - DELTA_STAR, height=0.3, color=INK_2, zorder=3)
        ax.barh(-1.38, test.w, left=test.leader - test.w, height=0.3, color=INK, zorder=3)
        ax.text(test.leader - DELTA_STAR - 0.1, -0.92, "δ*", ha="right", va="center", fontsize=12, color=INK)
        ax.text(test.leader - test.w - 0.1, -1.38, "w", ha="right", va="center", fontsize=12, color=INK)
        if n == 1000:
            assert test.w > DELTA_STAR and "A" not in test.ruled_out, "the guide narrates this outcome"
            ax.text(test.means[1] + 0.3, 1, "still kept", va="center", fontsize=11.5, color=INK)
        if n == 15000:
            assert "B" in test.ruled_out, "the guide narrates this outcome"
            ax.text(TRUE_MEANS[2] - 0.3, 2, "now out", ha="right", va="center", fontsize=11.5, color=INK)
    fig.tight_layout(w_pad=1.2)
    return fig


def many_tests_figure() -> plt.Figure:
    tests = consecutive_tests()
    rows = len(tests)
    assert all({"Control", "A"} <= test.ruled_out for test in tests), "the guide narrates this outcome"
    assert sum(test.winner == "C" for test in tests) == 2, "the guide narrates this outcome"
    assert sum("B" not in test.ruled_out for test in tests) == 2, "the guide narrates this outcome"

    fig, ax = plt.subplots(figsize=(9.2, 5.9))
    left, right = 47.9, 56.55
    y_truth, y_strip = rows + 0.95, rows + 2.35
    for x in (OUTER, INNER, BEST):
        ax.plot([x, x], [-0.5, y_strip], color=HAIRLINE, lw=0.8, zorder=1.5)
    strips = (
        (left, OUTER, RED_TINT, "clearly worse", 49.75),
        (OUTER, INNER, GRAY_TINT, "gray zone", 51.5),
        (INNER, BEST, BLUE_TINT, "good enough", 53.5),
    )
    for start, end, color, name, x in strips:
        ax.add_patch(Rectangle((start + 0.015, y_strip), end - start - 0.03, 1.0, fc=color, lw=0, zorder=2))
        ax.text(x, y_strip + 0.5, name, ha="center", va="center", color=INK, fontweight="bold", zorder=3)

    # Identity never rests on hue alone: D has its own shape and C and D differ in lightness.
    color = {"Control": RED, "A": RED, "B": INK, "C": C_BLUE, "D": D_BLUE}
    for x, name, arm in ((50, "Control, A", "A"), (51, "B", "B"), (53.5, "C", "C"), (54.5, "D", "D")):
        square = arm == "D"
        ax.plot(x, y_truth, "s" if square else "o", ms=8 if square else 8.5, mfc="white", mec=color[arm], mew=1.8, zorder=3)
        before = name == "Control, A"
        ax.text(x + (-0.17 if before else 0.17), y_truth, name, ha="right" if before else "left", va="center", fontsize=12, color=INK)

    for row, test in enumerate(tests):
        y = rows - 1 - row
        ax.barh(y, test.cut - left, left=left, height=0.84, color=BAR, lw=0, zorder=1)
        ax.plot([test.cut] * 2, [y - 0.42, y + 0.42], color=INK, lw=2, solid_capstyle="butt", zorder=2)
        for arm, mean in zip(ARMS, test.means):
            square = arm == "D"
            ax.plot(mean, y, "s" if square else "o", ms=6.2 if square else 6.8, mfc=color[arm], mec="white", mew=0.9, zorder=4)
        notes = ([f"{test.winner} leads"] if test.winner != "D" else []) + (["B kept"] if "B" not in test.ruled_out else [])
        if notes:
            ax.text(55.45, y, ", ".join(notes), va="center", fontsize=11.5, color=INK_2)

    last = tests[-1]
    draw_bracket(ax, last.cut, last.leader, -1.15)
    ax.text((last.cut + last.leader) / 2, -1.3, "w: from cut line to leader", ha="center", va="top", fontsize=11.5, color=INK)
    ax.text(left + 0.04, -1.15, "Gray bar: ruled out in that test", va="center", fontsize=11.5, color=INK_2)
    ax.set_xlim(left, right)
    ax.set_ylim(-2.1, y_strip + 1.05)
    ax.set_yticks([y_truth, (rows - 1) / 2])
    ax.set_yticklabels(["True means", f"{rows} simulated\ntests"])
    set_minutes_axis(ax)
    ax.xaxis.labelpad = 6
    ax.spines["bottom"].set_visible(False)
    fig.tight_layout()
    return fig


def main(output_dir: Path = DEFAULT_OUTPUT_DIR) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = {"one_test": one_test_figure, "more_users": more_users_figure, "many_tests": many_tests_figure}
    with plt.rc_context(STYLE):
        for name, make_figure in figures.items():
            fig = make_figure()
            fig.savefig(output_dir / f"sample_size_{name}.svg", facecolor="white", metadata={"Date": None})
            plt.close(fig)


if __name__ == "__main__":
    main()
