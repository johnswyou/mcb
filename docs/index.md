# `mcb`

`mcb` packages Hsu's Multiple Comparisons with the Best (MCB) procedures for the unbalanced one-way normal model.

Source repository: <https://github.com/johnswyou/mcb>

The library focuses on two complementary workflows:

- `constrained_mcb` for simultaneous intervals on $\mu_i - \max_{j \neq i} \mu_j$ and the candidate-best set.
- `unconstrained_mcb_edwards_hsu` for simultaneous intervals on $\mu_i - \max_j \mu_j$, the gap to the global best.

```{toctree}
:maxdepth: 2
:caption: Guide

quickstart
example_workflow
api
mcb
mcb_critical_values
```

## Highlights

- Accept either raw per-arm samples or summary statistics.
- Compute Hsu's constrained MCB intervals for ranking plausible winners.
- Compute the unbalanced Edwards-Hsu unconstrained intervals for gaps to the global best.
- Plot simultaneous intervals with the optional `mcb.plotting` helper.

## Local documentation build

From a local clone, install the docs and plotting extras and build the HTML output:

```bash
python3 -m pip install -e '.[docs,plot]'
python3 -m sphinx -b html docs docs/_build/html
```

The generated site will be written to `docs/_build/html`.
