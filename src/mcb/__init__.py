"""Public package interface for the MCB library."""

from .core import (
    MCBResult,
    OneWaySummary,
    constrained_critical_value,
    constrained_mcb,
    pooled_mse_from_group_stats,
    pooled_mse_from_samples,
    print_result,
    summarize_from_means_ns_and_mse,
    summarize_from_samples,
    unconstrained_edwards_hsu_critical_value,
    unconstrained_mcb_edwards_hsu,
)

__all__ = [
    "MCBResult",
    "OneWaySummary",
    "constrained_critical_value",
    "constrained_mcb",
    "pooled_mse_from_group_stats",
    "pooled_mse_from_samples",
    "print_result",
    "summarize_from_means_ns_and_mse",
    "summarize_from_samples",
    "unconstrained_edwards_hsu_critical_value",
    "unconstrained_mcb_edwards_hsu",
]

__version__ = "0.1.0"
