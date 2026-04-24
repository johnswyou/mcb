# `mcb`

`mcb` implements Hsu's Multiple Comparisons with the Best (MCB) procedures for the unbalanced one-way normal model.

The package is organized around two inferential tasks:

- `constrained_mcb`: simultaneous intervals for $\mu_i - \max_{j \neq i} \mu_j$ and a candidate-best set.
- `unconstrained_mcb_edwards_hsu`: Edwards-Hsu simultaneous intervals for $\mu_i - \max_j \mu_j$, the gap from arm $i$ to the global best.

The statistical reference is Hsu, *Multiple Comparisons: Theory and Methods* (1996), Chapter 4, with critical-value computation following Appendix D.2 and regression checks against Appendix E.

Source repository: <https://github.com/johnswyou/mcb>

```{toctree}
:maxdepth: 2
:caption: Guide

quickstart
mcb
model_notation
constrained_mcb
edwards_hsu
sat_example
mcb_critical_values
assumptions
example_workflow
api
```

## Highlights

- Accept either raw per-arm samples or summary statistics.
- Compute Hsu's constrained MCB intervals for ranking plausible winners.
- Compute unbalanced Edwards-Hsu unconstrained intervals for gaps to the global best.
- Plot simultaneous intervals with the optional `mcb.plotting` helper.
- Build documentation locally with the same Sphinx configuration used by Read the Docs.

## Local documentation build

From a local clone, install the docs and plotting extras and build the HTML output:

```bash
python3 -m pip install -e '.[docs,plot,test]'
python3 -m sphinx -W -b html docs docs/_build/html
```

The generated site will be written to `docs/_build/html`.
