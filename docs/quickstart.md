# Quickstart

## Install from a local clone

Install the core package in editable mode:

```bash
python3 -m pip install -e .
```

Install the optional plotting helpers as well when you want figures:

```bash
python3 -m pip install -e '.[plot]'
```

## Start from summary statistics

If you already know the per-arm means, sample sizes, and pooled within-arm variance estimate, start with `summarize_from_means_ns_and_mse`:

```python
from mcb import constrained_mcb, pooled_mse_from_group_stats, summarize_from_means_ns_and_mse

mse = pooled_mse_from_group_stats(
    ns=[103, 31, 122],
    variances=[82.52**2, 82.52**2, 82.52**2],
)

summary = summarize_from_means_ns_and_mse(
    means=[619.0, 629.0, 575.0],
    ns=[103, 31, 122],
    mse=mse,
    labels=["Computer science", "Engineering", "Other"],
)

result = constrained_mcb(summary, alpha=0.10)

print(result.candidate_best_set)
print(result.as_dict())
```

Use `result.candidate_best_set` to see which variants still have a strictly positive constrained upper bound. The per-arm estimates and intervals are available on the dataclass fields or through `result.as_dict()`.

## Start from raw samples

If you have raw observations instead, let `mcb` compute the pooled summary for you:

```python
from mcb import constrained_mcb, summarize_from_samples

summary = summarize_from_samples(
    {
        "Computer science": [612.0, 625.0, 620.0, 619.0],
        "Engineering": [635.0, 621.0, 628.0, 632.0],
        "Other": [571.0, 580.0, 576.0, 573.0],
    }
)

result = constrained_mcb(summary, alpha=0.10)
```

## Choose the procedure

- Use `constrained_mcb` when you want intervals for $\theta_i = \mu_i - \max_{j \neq i} \mu_j$ and a candidate-best set.
- Use `unconstrained_mcb_edwards_hsu` when you want intervals for $\eta_i = \mu_i - \max_j \mu_j$, the distance below the global best.

Both procedures return the same `MCBResult` dataclass, so you can inspect or display them the same way.
