# mcb

[![Documentation Status](https://readthedocs.org/projects/mcb/badge/?version=latest)](https://mcb.readthedocs.io/en/latest/?badge=latest)

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
    variances=[86.0**2, 67.0**2, 83.0**2],
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

## Planning sample size

Choose the per-arm sample size of a balanced experiment with Hsu's Appendix C computation: the smallest `n` such that, with probability at least `target`, every arm more than `2 * delta_star` below the true best is ruled out and the sample-best arm is correctly shown to be within `delta_star` of the true best.

```python
from mcb import constrained_mcb_sample_size, unconstrained_mcb_edwards_hsu_sample_size

# Hsu (1996) Appendix C.4: 8 alloys, planning sigma 3000 psi, margin delta* = 5000 psi.
constrained_mcb_sample_size(k=8, delta_star=5000, sigma=3000, alpha=0.05, target=0.90)              # 6
unconstrained_mcb_edwards_hsu_sample_size(k=8, delta_star=5000, sigma=3000, alpha=0.05, target=0.90)  # 8
```

See the [sample-size guide](https://mcb.readthedocs.io/en/latest/sample_size.html) for what the guarantee means.

## Documentation

The Sphinx documentation source lives in `docs/` and is configured for Read the Docs through `.readthedocs.yaml`.
The source repository is https://github.com/johnswyou/mcb.
The hosted documentation is available at https://mcb.readthedocs.io/en/latest/.
Build it locally with:

```bash
python3 -m pip install -e '.[docs,plot]'
python3 -m sphinx -b html docs docs/_build/html
```

## Tests

The regression suite uses `pytest` and validates the critical value implementation against Appendix E tables extracted from Hsu (1996).

```bash
python3 -m pip install -e '.[test]'
python3 -m pytest
```
