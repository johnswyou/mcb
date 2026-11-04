"""Public package interface for the MCB library."""

from .core import (
    MCBResult,
    OneWaySummary,
    constrained_critical_value,
    constrained_mcb,
    constrained_mcb_correct_useful_probability,
    constrained_mcb_sample_size,
    pooled_mse_from_group_stats,
    pooled_mse_from_samples,
    print_result,
    summarize_from_means_ns_and_mse,
    summarize_from_samples,
    unconstrained_edwards_hsu_critical_value,
    unconstrained_mcb_edwards_hsu,
    unconstrained_mcb_edwards_hsu_correct_useful_probability,
    unconstrained_mcb_edwards_hsu_sample_size,
)

__all__ = [
    "MCBResult",
    "OneWaySummary",
    "constrained_critical_value",
    "constrained_mcb",
    "constrained_mcb_correct_useful_probability",
    "constrained_mcb_sample_size",
    "pooled_mse_from_group_stats",
    "pooled_mse_from_samples",
    "print_result",
    "summarize_from_means_ns_and_mse",
    "summarize_from_samples",
    "unconstrained_edwards_hsu_critical_value",
    "unconstrained_mcb_edwards_hsu",
    "unconstrained_mcb_edwards_hsu_correct_useful_probability",
    "unconstrained_mcb_edwards_hsu_sample_size",
]

__version__ = "0.1.0"
