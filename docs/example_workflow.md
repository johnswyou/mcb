# Example workflow

The repository includes a worked SAT-score example in `examples/demo_hsu_mcb_sat.py`. It mirrors the summary-statistics workflow used in Hsu's chapter and shows both supported MCB procedures.

## Reproduce the SAT example

```python
from mcb import constrained_mcb, summarize_from_means_ns_and_mse, unconstrained_mcb_edwards_hsu

summary = summarize_from_means_ns_and_mse(
    means=[619.0, 629.0, 575.0],
    ns=[103, 31, 122],
    mse=82.52**2,
    labels=["Computer science", "Engineering", "Other"],
)

constrained = constrained_mcb(summary, alpha=0.10)
unconstrained = unconstrained_mcb_edwards_hsu(summary, alpha=0.10)
```

The constrained result answers, "which variants can still plausibly be the best?" The unconstrained Edwards-Hsu result answers, "how far below the global best is each variant?"

## Plot the simultaneous intervals

Install the plotting extra first:

```bash
python3 -m pip install -e '.[plot]'
```

Then render the intervals with `plot_simultaneous_intervals`:

```python
import matplotlib.pyplot as plt

from mcb.plotting import plot_simultaneous_intervals

fig, ax = plt.subplots(figsize=(9, 4.8))
plot_simultaneous_intervals(constrained, ax=ax, title="SAT example: constrained MCB")
fig.tight_layout()
fig.savefig("sat_constrained_mcb.png", dpi=160)
```

Running the example script directly will save both constrained and unconstrained plots under `examples/output/`:

```bash
python3 examples/demo_hsu_mcb_sat.py
```
