# SAT worked example

This page reproduces the SAT mathematics example in Hsu (1996), Sections 4.1.5 and 4.2.5, at $\alpha = 0.10$. Students are grouped by their major at the start of their second year:

| Major | $n_i$ | $\hat{\mu}_i$ | Sample SD |
| --- | ---: | ---: | ---: |
| Computer science | 103 | 619 | 86 |
| Engineering | 31 | 629 | 67 |
| Other | 122 | 575 | 83 |

The error degrees of freedom are $\nu = 102 + 30 + 121 = 253$, and the pooled standard deviation is $\hat{\sigma} = 82.52$.

## Build the summary

```python
from mcb import (
    constrained_mcb,
    pooled_mse_from_group_stats,
    summarize_from_means_ns_and_mse,
    unconstrained_mcb_edwards_hsu,
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
```

## Constrained MCB

```python
constrained = constrained_mcb(summary, alpha=0.10)
```

The critical values $d^i$ and the intervals for $\theta_i = \mu_i - \max_{j \neq i}\mu_j$ are:

| Arm | $d^i$ | Estimate | Lower | Upper |
| --- | ---: | ---: | ---: | ---: |
| Computer science | 1.605 | -10.00 | -35.43 | 17.14 |
| Engineering | 1.505 | 10.00 | -17.14 | 35.43 |
| Other | 1.612 | -54.00 | -78.97 | 0.00 |

The bounds for `Other` (arm 3) show how the formulas in [Constrained MCB](constrained_mcb.md) are applied. The upper bound compares `Other` with each of the other arms, using its own critical value $d^3 = 1.612$:

$$
D_3^+ =
\left[\min\left\{
\begin{array}{l}
575 - 619 + 1.612\,(82.52)\sqrt{1/122 + 1/103} \\
575 - 629 + 1.612\,(82.52)\sqrt{1/122 + 1/31}
\end{array}
\right\}\right]^+
= \left[\min\{-26.21,\ -27.25\}\right]^+ = 0.
$$

`Other` is ruled out as the best. The same calculation gives positive upper bounds for the other two arms, so $G$ contains Computer science and Engineering. The lower bound compares `Other` with each arm in $G$, using that arm's critical value ($d^1 = 1.605$, $d^2 = 1.505$):

$$
D_3^- =
\min\left\{
\begin{array}{l}
575 - 619 - 1.605\,(82.52)\sqrt{1/122 + 1/103} \\
575 - 629 - 1.505\,(82.52)\sqrt{1/122 + 1/31}
\end{array}
\right\}
= \min\{-61.72,\ -78.97\} = -78.97.
$$

The conclusions at familywise error rate 0.10 are:

- `Other` is not the best. It is at most 78.97 points below the best.
- Computer science and Engineering could each be the best.
- Engineering is at most 17.14 points below the better of the other two arms, and Computer science at most 35.43 points below.

## Edwards-Hsu MCB

```python
unconstrained = unconstrained_mcb_edwards_hsu(summary, alpha=0.10)
```

The critical values $|d|^i$ and the intervals for $\eta_i = \mu_i - \max_j\mu_j$ are:

| Arm | $\lvert d\rvert^i$ | Estimate | Lower | Upper |
| --- | ---: | ---: | ---: | ---: |
| Computer science | 1.942 | -10.00 | -41.44 | 0.00 |
| Engineering | 1.860 | 0.00 | -22.82 | 0.00 |
| Other | 1.946 | -54.00 | -84.87 | -22.56 |

The screening set $S$ again contains Computer science and Engineering. `Other` is now shown to be between 22.56 and 84.87 points below the best, where constrained MCB showed only that it is at most 78.97 below. The larger critical values move the other lower bounds away from 0: Computer science is at most 41.44 points below the best, against 35.43 under constrained MCB.

## Plot the intervals

`plot_simultaneous_intervals` draws the intervals of either result. It needs the `plot` extra (see [Quickstart](quickstart.md)).

```python
import matplotlib.pyplot as plt

from mcb.plotting import plot_simultaneous_intervals

fig, ax = plt.subplots(figsize=(9, 4.8))
plot_simultaneous_intervals(constrained, ax=ax, title="SAT example: constrained MCB")
fig.tight_layout()
fig.savefig("sat_constrained_mcb.png", dpi=160)
```

The script `examples/demo_hsu_mcb_sat.py` in the repository runs this example and saves a plot for each procedure under `examples/output/`.

## Differences from the printed textbook values

The constrained results above agree with Section 4.1.5 of the book to within 0.01, the rounding in the book's arithmetic.

The Edwards-Hsu critical values differ from Section 4.2.5 in the second decimal place:

| Arm | Book $\lvert d\rvert^i$ | Package $\lvert d\rvert^i$ |
| --- | ---: | ---: |
| Computer science | 1.933 | 1.942 |
| Engineering | 1.852 | 1.860 |
| Other | 1.937 | 1.946 |

Hsu's programs replaced $\nu$ by $\infty$ whenever $\nu > 240$ (Appendix D.2), and the book's Edwards-Hsu values are the $\nu = \infty$ values. The package uses the actual $\nu = 253$, which gives slightly larger critical values and slightly wider intervals. Calling `unconstrained_edwards_hsu_critical_value` with `nu=math.inf` reproduces the book's values.

The book says the same replacement was made in the constrained example, but the constrained values it prints (1.605, 1.505, 1.612) are the $\nu = 253$ values. With $\nu = \infty$ they would be 1.600, 1.500 and 1.606.
