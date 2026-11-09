# Quickstart

## Install

From a local clone:

```bash
python3 -m pip install -e .
```

Add the `plot` extra for the plotting helper:

```bash
python3 -m pip install -e '.[plot]'
```

## Analyze summary statistics

`mcb` needs each arm's sample mean and sample size, and the pooled within-arm variance (MSE). Larger responses are treated as better.

This example uses the SAT mathematics scores from Hsu (1996, Section 4.1.5): three groups of students with sample sizes 103, 31 and 122. The source reports a standard deviation per group, so the pooled MSE is computed from those first.

```python
from mcb import (
    constrained_mcb,
    pooled_mse_from_group_stats,
    print_result,
    summarize_from_means_ns_and_mse,
)

mse = pooled_mse_from_group_stats(
    ns=[103, 31, 122],
    variances=[86.0**2, 67.0**2, 83.0**2],
)

summary = summarize_from_means_ns_and_mse(
    means=[619.0, 629.0, 575.0],
    ns=[103, 31, 122],
    mse=mse,
    labels=["Computer science", "Engineering", "Other"],
)

result = constrained_mcb(summary, alpha=0.10)
print_result(result, digits=2)
```

```text
Hsu constrained MCB (unbalanced one-way)
Parameter: mu_i - max_{j != i} mu_j
Familywise alpha: 0.1
Candidate-best set: ('Computer science', 'Engineering')
Critical values:
  Computer science: 1.61
  Engineering: 1.50
  Other: 1.61
Intervals:
  Computer science   estimate=-10.00  CI=[-35.43, 17.14]
  Engineering        estimate=10.00  CI=[-17.14, 35.43]
  Other              estimate=-54.00  CI=[-78.97, 0.00]
```

Each interval is for an arm's true mean minus the largest of the other arms' true means. `Other` has an upper bound of 0, so it is ruled out as the best. The other two arms could each be the best: Engineering is at most 17.14 points below the better of the other two and at most 35.43 points above it. All three intervals contain their true values at once with probability at least 0.90.

The same values are fields of the returned `MCBResult` (`labels`, `estimates`, `lower`, `upper`, `critical_values`, `candidate_best_set`), and `result.as_dict()` returns the intervals as a dictionary keyed by label.

## Analyze raw observations

Given the observations of each arm, `summarize_from_samples` computes the means, sample sizes and pooled MSE:

```python
from mcb import constrained_mcb, summarize_from_samples

summary = summarize_from_samples(
    {
        "A": [612.0, 625.0, 620.0, 619.0],
        "B": [635.0, 621.0, 628.0, 632.0],
        "C": [571.0, 580.0, 576.0, 573.0],
    }
)

result = constrained_mcb(summary, alpha=0.10)
```

## Next steps

- [Overview](mcb.md) explains what the intervals mean and when to use `unconstrained_mcb_edwards_hsu` instead of `constrained_mcb`.
- [SAT worked example](sat_example.md) traces the numbers above back to the formulas.
- [Sample size](sample_size.md) chooses the number of observations per arm before an experiment runs.
