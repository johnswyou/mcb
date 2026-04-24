# SAT Worked Example

This page reproduces the SAT mathematics example from Hsu (1996), Sections 4.1.5 and 4.2.5, using the package API.

The data summarize SAT mathematics scores by second-year major:

| Major | Label | $n_i$ | $\hat{\mu}_i$ | Sample SD |
| --- | ---: | ---: | ---: | ---: |
| Computer science | 1 | 103 | 619 | 86 |
| Engineering | 2 | 31 | 629 | 67 |
| Other | 3 | 122 | 575 | 83 |

The pooled standard deviation from the three sample variances is approximately 82.52, with

$$
\nu = (103-1) + (31-1) + (122-1) = 253.
$$

## Build the Summary

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

The finite-df constrained critical values are:

| Arm | $d_i$ |
| --- | ---: |
| Computer science | 1.605 |
| Engineering | 1.505 |
| Other | 1.612 |

The simultaneous intervals for

$$
\theta_i = \mu_i - \max_{j \neq i}\mu_j
$$

are:

| Arm | Estimate | Lower | Upper |
| --- | ---: | ---: | ---: |
| Computer science | -10.00 | -35.43 | 17.14 |
| Engineering | 10.00 | -17.14 | 35.43 |
| Other | -54.00 | -78.97 | 0.00 |

The constrained candidate-best set is:

```python
("Computer science", "Engineering")
```

`Other` has upper bound 0, so it is ruled out as the true best at familywise level 0.10. Computer science and Engineering remain compatible with being best.

## Edwards-Hsu Unconstrained MCB

```python
unconstrained = unconstrained_mcb_edwards_hsu(summary, alpha=0.10)
```

The package computes finite-df critical values because $\nu = 253 < 10^6$:

| Arm | $|d|_i$ |
| --- | ---: |
| Computer science | 1.942 |
| Engineering | 1.860 |
| Other | 1.946 |

The direct Edwards-Hsu intervals for

$$
\eta_i = \mu_i - \max_j\mu_j
$$

are:

| Arm | Estimate | Lower | Upper |
| --- | ---: | ---: | ---: |
| Computer science | -10.00 | -41.44 | 0.00 |
| Engineering | 0.00 | -22.83 | 0.00 |
| Other | -54.00 | -84.87 | -22.56 |

The Edwards-Hsu candidate-best set is again:

```python
("Computer science", "Engineering")
```

`Other` is inferred to be below the global best because its upper bound is strictly negative.

## About the Textbook's Infinite-df Shortcut

Appendix D.2 describes Hsu's original computer code as switching to the infinite-degrees-of-freedom calculation when $\nu > 240$. The package does not use that historical shortcut for $\nu = 253$; it computes the finite-df values directly.

For the Edwards-Hsu SAT example, the infinite-df shortcut gives the textbook's displayed critical values:

| Arm | Historical $\nu=\infty$ $|d|_i$ | Package finite-df $|d|_i$ |
| --- | ---: | ---: |
| Computer science | 1.933 | 1.942 |
| Engineering | 1.852 | 1.860 |
| Other | 1.937 | 1.946 |

The finite-df intervals are slightly wider. This is intentional: finite $\nu$ below the package's high-df stability switch is not collapsed to Hsu's `dinfnu = 240` approximation.
