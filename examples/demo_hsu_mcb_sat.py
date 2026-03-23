"""Demo of Hsu's unbalanced MCB procedures using Hsu's SAT-score example.

This script reproduces the example from Chapter 4 using summary statistics:
- Variant 1: Computer science
- Variant 2: Engineering
- Variant 3: Other
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from mcb import constrained_mcb, print_result, summarize_from_means_ns_and_mse, unconstrained_mcb_edwards_hsu
from mcb.plotting import plot_simultaneous_intervals

OUTDIR = Path(__file__).resolve().parent / "output"


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)

    summary = summarize_from_means_ns_and_mse(
        means=[619.0, 629.0, 575.0],
        ns=[103, 31, 122],
        mse=82.52**2,
        labels=["Computer science", "Engineering", "Other"],
    )

    alpha = 0.10

    constrained = constrained_mcb(summary, alpha=alpha)
    unconstrained = unconstrained_mcb_edwards_hsu(summary, alpha=alpha)

    print("=" * 80)
    print("CONSTRAINED MCB")
    print_result(constrained, digits=3)
    print()
    print("=" * 80)
    print("UNCONSTRAINED EDWARDS-HSU MCB")
    print_result(unconstrained, digits=3)

    fig1, ax1 = plt.subplots(figsize=(9, 4.8))
    plot_simultaneous_intervals(
        constrained,
        ax=ax1,
        title="SAT example: constrained MCB (mu_i - max_{j != i} mu_j)",
    )
    fig1.tight_layout()
    fig1.savefig(OUTDIR / "sat_constrained_mcb.png", dpi=160)

    fig2, ax2 = plt.subplots(figsize=(9, 4.8))
    plot_simultaneous_intervals(
        unconstrained,
        ax=ax2,
        title="SAT example: Edwards-Hsu unconstrained MCB (mu_i - max_j mu_j)",
    )
    fig2.tight_layout()
    fig2.savefig(OUTDIR / "sat_unconstrained_edwards_hsu.png", dpi=160)

    print()
    print(f"Saved plots to: {OUTDIR / 'sat_constrained_mcb.png'}")
    print(f"Saved plots to: {OUTDIR / 'sat_unconstrained_edwards_hsu.png'}")


if __name__ == "__main__":
    main()
