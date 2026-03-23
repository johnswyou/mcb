# mcb

`mcb` packages Hsu's Multiple Comparisons with the Best procedures for the unbalanced one-way normal model.

## Install

Install the library in editable mode:

```bash
python3 -m pip install -e .
```

Install the plotting and test extras when needed:

```bash
python3 -m pip install -e '.[plot,test]'
```

## Usage

When you already have per-arm means, sample sizes, and within-arm variances:

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
```

When you have raw sample-level data instead:

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

Plotting is available through the optional module:

```python
from mcb.plotting import plot_simultaneous_intervals
```

See `examples/demo_hsu_mcb_sat.py` for the SAT-score example.

## Tests

The regression suite uses `pytest` and validates the critical value implementation against Appendix E tables extracted from Hsu (1996).

```bash
python3 -m pip install -e '.[test]'
python3 -m pytest
```
