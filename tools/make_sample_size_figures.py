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
LUCK_BAR, RULED_OUT_DOT = "#b9b7ae", "#b5b3aa"
ON_WHITE = dict(fc="white", ec="none", pad=1.5)

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
    rng = np.random.default_rng(245003)
    errors = rng.normal(size=K)
    sigma_hat = SIGMA * np.sqrt(rng.chisquare(NU) / NU)
    test = analyze(TRUE_MEANS + errors * SIGMA / np.sqrt(n), n, sigma_hat)
    # The guide teaches the correctness event from this test, so the test must satisfy it.
    true_gaps, gaps_in_test = BEST - TRUE_MEANS, test.means[TRUE_MEANS.argmax()] - test.means
    assert np.all(true_gaps - gaps_in_test < test.w), "a chosen test must satisfy the correctness event"
    return test


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


def draw_zone_strip(ax: plt.Axes, left: float, y: float, height: float, clearly_worse_x: float) -> None:
    strips = (
        (left, OUTER, RED_TINT, "clearly worse", clearly_worse_x),
        (OUTER, INNER, GRAY_TINT, "gray zone", 51.5),
        (INNER, BEST, BLUE_TINT, "good enough", 53.5),
    )
    for start, end, color, name, x in strips:
        ax.add_patch(Rectangle((start + 0.015, y), end - start - 0.03, height, fc=color, lw=0, zorder=2))
        ax.text(x, y + height / 2, name, ha="center", va="center", color=INK, fontweight="bold", zorder=3)


def draw_ring(ax: plt.Axes, x: float, y: float, size: float = 11) -> None:
    ax.plot(x, y, "o", ms=size, mfc="white", mec=INK_2, mew=1.5, zorder=5)


def draw_dot(ax: plt.Axes, x: float, y: float, size: float = 11) -> None:
    ax.plot(x, y, "o", ms=size, mfc=INK, mec="white", mew=1.2, zorder=5)


def draw_gap(ax: plt.Axes, start: float, end: float, y: float, color: str = INK_2) -> None:
    ax.plot([start, end], [y, y], color=color, lw=5, solid_capstyle="butt", zorder=3)


def draw_links(ax: plt.Axes, pairs: tuple, top: float, bottom: float) -> None:
    for x_true, x_test in pairs:
        ax.plot([x_true, x_test], [top, bottom], color=MUTED, lw=0.9, ls=(0, (3, 2)), zorder=2)


def centered(ax: plt.Axes, x: float, y: float, text: str, **style) -> None:
    ax.text(x, y, text, ha="center", va="center", color=INK, **style)


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
    draw_zone_strip(ax, left, y_strip, 1.0, 49.75)

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


def luck_in_one_test_figure() -> plt.Figure:
    test = followed_test(PLANNED_N)
    a_true, c_true = TRUE_MEANS[1], TRUE_MEANS[3]
    a, c, d = test.means[1], test.means[3], test.means[4]
    a_true_gap, c_true_gap = BEST - a_true, BEST - c_true
    a_gap_in_test = d - a
    a_luck, c_luck = a_true_gap - a_gap_in_test, c_true_gap + c - d
    largest_luck = np.max((BEST - TRUE_MEANS) - (d - test.means))
    narrated = [test.w, 2 * test.w, a_gap_in_test, a_luck, a - a_true, BEST - d, c - d, c_luck, largest_luck]
    assert test.winner == "C" and [round(float(x), 2) for x in narrated] == [1.97, 3.93, 3.10, 1.40, 0.42, 0.98, 0.21, 1.21, 1.87], (
        "the guide narrates this outcome"
    )

    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    left, right = 49.1, 58.25
    y_strip = 9.15
    for x in (OUTER, INNER, BEST):
        ax.plot([x, x], [-0.15, y_strip], color=HAIRLINE, lw=0.8, zorder=1)
    draw_zone_strip(ax, left, y_strip, 0.62, (left + OUTER) / 2)

    a_true_row, a_test_row, a_leader_row, c_true_row, c_test_row = 7.7, 6.45, 5.2, 2.6, 1.35
    heading = dict(color=INK, fontweight="bold", va="center", bbox=ON_WHITE, zorder=4)
    ax.text(left + 0.02, 8.6, "A is more than 2w below D, so it is ruled out", **heading)
    ax.text(left + 0.02, 3.5, "C won, so its true gap to D is less than w", **heading)

    draw_gap(ax, a_true, BEST, a_true_row)
    draw_links(ax, ((a_true, a), (BEST, d)), a_true_row - 0.2, a_test_row + 0.2)
    draw_gap(ax, a_true, a, a_test_row, LUCK_BAR)
    draw_gap(ax, d, BEST, a_test_row, LUCK_BAR)
    draw_gap(ax, a, d, a_test_row)
    draw_gap(ax, a, c, a_leader_row)
    draw_gap(ax, c_true, BEST, c_true_row)
    draw_links(ax, ((c_true, c), (BEST, d)), c_true_row - 0.2, c_test_row + 0.2)
    draw_gap(ax, d, c, c_test_row)
    for x, y in ((a_true, a_true_row), (BEST, a_true_row), (c_true, c_true_row), (BEST, c_true_row)):
        draw_ring(ax, x, y)
    for x, y in ((a, a_test_row), (d, a_test_row), (a, a_leader_row), (c, a_leader_row), (d, c_test_row), (c, c_test_row)):
        draw_dot(ax, x, y)

    above, below = 0.42, -0.42
    arm_names = (
        (a_true, a_true_row + above, "A"),
        (BEST, a_true_row + above, "D"),
        (a, a_test_row + below, "A"),
        (d, a_test_row + below, "D"),
        (a, a_leader_row + below, "A"),
        (c, a_leader_row + below, "C leads"),
        (c_true, c_true_row + above, "C"),
        (BEST, c_true_row + above, "D"),
        (d - 0.02, c_test_row + below, "D"),
        (c + 0.02, c_test_row + below, "C"),
    )
    for x, y, name in arm_names:
        centered(ax, x, y, name)
    centered(ax, 52.25, a_true_row + 0.3, f"true gap {a_true_gap:.2f}", fontsize=12)
    centered(ax, (a + d) / 2, a_test_row + 0.3, f"gap in the test {a_gap_in_test:.2f}", fontsize=12)
    ax.plot([test.cut] * 2, [a_leader_row - 0.38, a_leader_row + 0.38], color=INK, lw=2.6, zorder=6, solid_capstyle="butt")
    centered(ax, test.cut, a_leader_row - 0.62, "cut line", fontsize=11.5)
    draw_bracket(ax, test.cut, c, a_leader_row + 0.42)
    centered(ax, (test.cut + c) / 2, a_leader_row + 0.68, f"w = {test.w:.2f}", fontsize=12)
    side_label = dict(ha="right", va="center", color=INK, fontsize=12, bbox=ON_WHITE, zorder=4)
    ax.text(53.2, c_true_row, f"true gap {c_true_gap:.2f}", **side_label)
    ax.text(d - 0.3, c_test_row, f"C ahead by {c - d:.2f}", **side_label)

    note = dict(va="center", color=INK, fontsize=12)
    x_note = 54.95
    ax.text(x_note, a_true_row, f"more than 2w = {2 * test.w:.2f}", **note)
    ax.text(
        x_note,
        a_test_row,
        f"luck = {a_true_gap:.2f} − {a_gap_in_test:.2f} = {a_luck:.2f}, less than w,\nso the gap in the test is more than w",
        **note,
    )
    ax.text(x_note, a_leader_row, "the leader is at or ahead of D, so A is\nmore than w behind it and is ruled out", **note)
    ax.text(x_note, c_true_row, "less than w, and w is less than δ*", **note)
    ax.text(
        x_note,
        c_test_row,
        f"luck = {c_true_gap:.2f} + {c - d:.2f} = {c_luck:.2f}, less than w.\nA winner's luck is at least its true gap.",
        **note,
    )

    y_key = y_strip + 0.31
    key = dict(va="center", color=INK, fontsize=11.5)
    draw_ring(ax, 54.98, y_key, 9)
    ax.text(55.1, y_key, "true mean", **key)
    draw_dot(ax, 56.02, y_key, 9)
    ax.text(56.14, y_key, "sample mean", **key)
    draw_gap(ax, 54.9, 55.14, y_key - 0.62, LUCK_BAR)
    ax.text(55.2, y_key - 0.62, "luck", **key)

    ax.set_xlim(left, right)
    ax.set_ylim(0.1, y_strip + 0.7)
    ax.set_yticks([a_true_row, a_test_row, a_leader_row, c_true_row, c_test_row])
    ax.set_yticklabels(["True means", "In the test", "Against\nthe leader", "True means", "In the test"])
    set_minutes_axis(ax)
    ax.set_xlabel("Minutes watched per user", color=INK_2, labelpad=6, x=0.33)
    ax.spines["bottom"].set_bounds(left, BEST + 0.3)
    fig.tight_layout()
    return fig


def dot_columns(x: np.ndarray, width: float, anchor: float, rank: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Column edges fall on the row's thresholds, so a dot never straddles a dashed line.
    columns = np.floor((x - anchor) / width).astype(int)
    level = np.zeros(len(x))
    for column in np.unique(columns):
        members = np.where(columns == column)[0]
        level[members[np.argsort(rank[members], kind="stable")]] = np.arange(len(members))
    return anchor + (columns + 0.5) * width, level


def luck_needed_figure() -> plt.Figure:
    tests = consecutive_tests(200)
    means = np.array([test.means for test in tests])
    each_w = np.array([test.w for test in tests])
    w = float(each_w.mean())
    arms = ARMS[:-1]  # D is the arm that luck is measured against
    true_gaps = BEST - TRUE_MEANS[:-1]
    luck = true_gaps - (means[:, [-1]] - means[:, :-1])
    correct = (luck < each_w[:, None]).all(axis=1)
    kept = {arm: np.array([arm not in test.ruled_out for test in tests]) for arm in arms}
    won = {arm: np.array([test.winner == arm for test in tests]) for arm in arms}
    gray_past_b = (luck[:, 2] > true_gaps[2] - w) & ~kept["B"]
    black_before_c = (luck[:, 3] < true_gaps[3] - w) & kept["C"]
    assert correct.sum() == 191 and gray_past_b.sum() == 1 and black_before_c.sum() == 1, "the guide narrates this outcome"
    assert [round(float(x), 2) for x in (each_w.min(), w, each_w.max())] == [1.95, 1.98, 2.02], "the guide narrates this outcome"
    lead_over_d = means[gray_past_b, 3] - means[gray_past_b, 4]
    assert round(float(lead_over_d[0]), 2) == 1.52, "the guide narrates this outcome"
    conclusions_held = [{"Control", "A"} <= test.ruled_out and test.winner in ("C", "D") for test in tests]
    assert np.all(conclusions_held), "the guide narrates this outcome"

    fig, ax = plt.subplots(figsize=(10.4, 5.9))
    left, right = -3.0, 4.85
    dot_size, level_height = 4.7, 0.083
    baseline = {"C": 3.72, "B": 2.48, "A": 1.24, "Control": 0.0}
    verdicts = {
        "C": ("good enough", BLUE_TINT, "can stay in,\ncan win"),
        "B": ("gray zone", GRAY_TINT, "can stay in,\ncannot win"),
        "A": ("clearly worse", RED_TINT, "cannot stay in"),
        "Control": ("clearly worse", RED_TINT, "cannot stay in"),
    }

    top = baseline["C"] + 1.2
    ax.add_patch(Rectangle((w, -0.1), right - w, top + 0.1, fc=BAR, lw=0, zorder=0))
    ax.plot([w, w], [-0.1, top], color=INK, lw=2.2, zorder=4, solid_capstyle="butt")
    ax.text(w, top + 0.05, f"w, about {w:.2f}", ha="center", va="bottom", color=INK, fontweight="bold")
    above_rows = dict(va="center", color=INK_2, fontsize=11.5)
    ax.text(w + 0.1, top - 0.12, f"{len(tests) - correct.sum()} of the {len(tests)} tests have a dot here", ha="left", **above_rows)
    ax.text(w - 0.1, top - 0.12, "in a correct test, every dot is left of w", ha="right", **above_rows)

    dashes = (0, (4, 2))
    for i, arm in enumerate(arms):
        base = baseline[arm]
        x, level = dot_columns(luck[:, i], w / 24, true_gaps[i], kept[arm].astype(int) + won[arm])
        y = base + level * level_height
        out, stayed, winner = ~kept[arm], kept[arm] & ~won[arm], won[arm]
        ax.plot([left, right], [base - 0.06, base - 0.06], color=AXIS, lw=0.8, zorder=1)
        ax.plot(x[out], y[out], "o", ms=dot_size, mfc=RULED_OUT_DOT, mec="white", mew=0.3, ls="", zorder=3)
        ax.plot(x[stayed], y[stayed], "o", ms=dot_size, mfc=INK, mec="white", mew=0.3, ls="", zorder=3)
        ax.plot(x[winner], y[winner], "D", ms=dot_size * 0.95, mfc=D_BLUE, mec="white", mew=0.3, ls="", zorder=3)
        line_top = base + 0.72
        for needed, label in ((true_gaps[i] - w, "to stay in"), (true_gaps[i], "to win")):
            ax.plot([needed, needed], [base - 0.06, line_top], color=INK, lw=1.3, ls=dashes, zorder=2.5)
            ax.text(
                needed,
                line_top + 0.03,
                f"{label} {needed:.2f}".replace("-", "−"),
                ha="center",
                va="bottom",
                color=INK,
                fontsize=10.5,
                bbox=dict(fc="white" if needed < w else BAR, ec="none", pad=0.8),
                zorder=2.6,
            )
        zone, tint, verdict = verdicts[arm]
        ax.text(
            right + 0.14,
            base + 0.62,
            zone,
            ha="left",
            va="center",
            color=INK,
            fontweight="bold",
            fontsize=11.5,
            bbox=dict(fc=tint, ec="none", pad=2.5),
            clip_on=False,
        )
        ax.text(right + 0.14, base + 0.42, verdict, ha="left", va="top", color=INK, fontsize=11.5, clip_on=False)

    ax.set_xlim(left, right)
    ax.set_ylim(-0.2, top + 0.34)
    ax.set_yticks([baseline[arm] + 0.32 for arm in arms])
    ax.set_yticklabels([f"{arm}\ntrue gap {gap:.2f}" for arm, gap in zip(arms, true_gaps)], fontsize=12)
    ax.set_xticks(range(-3, 5))
    ax.set_xticklabels([str(tick).replace("-", "−") for tick in range(-3, 5)])
    ax.set_xlabel("Luck: true gap to D minus gap in the test, in minutes", color=INK_2, labelpad=6)
    ax.plot([], [], "o", ms=6, mfc=RULED_OUT_DOT, mec="white", ls="", label="ruled out")
    ax.plot([], [], "o", ms=6, mfc=INK, mec="white", ls="", label="kept")
    ax.plot([], [], "D", ms=5.5, mfc=D_BLUE, mec="white", ls="", label="won the test")
    ax.plot([], [], color=INK, lw=1.3, ls=dashes, label="luck needed")
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(-0.01, 0.955),
        ncol=4,
        frameon=False,
        fontsize=11.5,
        handletextpad=0.1,
        columnspacing=1.0,
    )
    ax.spines["bottom"].set_visible(False)
    fig.tight_layout()
    return fig


def main(output_dir: Path = DEFAULT_OUTPUT_DIR) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = {
        "one_test": one_test_figure,
        "more_users": more_users_figure,
        "many_tests": many_tests_figure,
        "luck_in_one_test": luck_in_one_test_figure,
        "luck_needed": luck_needed_figure,
    }
    with plt.rc_context(STYLE):
        for name, make_figure in figures.items():
            fig = make_figure()
            fig.savefig(output_dir / f"sample_size_{name}.svg", facecolor="white", metadata={"Date": None})
            plt.close(fig)


if __name__ == "__main__":
    main()
