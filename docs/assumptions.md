# Assumptions and scope

## Model assumptions

The coverage guarantee, that all intervals contain their true values simultaneously with probability at least $1-\alpha$ (Hsu's Theorems 4.1.2 and 4.2.6), is exact under the [one-way normal model](model_notation.md):

- the observations are independent, within and across arms;
- every arm has the same variance $\sigma^2$;
- the errors are normally distributed.

The guarantee holds at every sample size. It is not a large-sample approximation.

## Use with A/B/n tests

In an A/B/n test the model usually holds only approximately, and the coverage is then approximate too.

- Normality. Per-user metrics such as revenue are often skewed or zero-inflated. With large samples the arm means are approximately normal, which is enough when the other two assumptions hold.
- Common variance. Large samples do not make up for unequal variances. Compare the per-arm sample variances before relying on the pooled MSE.
- Independence. Large samples do not make up for dependent observations either. Analyze the unit that was randomized: if users were randomized, use one value per user, not one per session.

## Ties

Hsu's proofs hold whether or not two arms can have exactly equal means. What changes is the reading of a constrained MCB bound of 0 (Hsu 1996, Sections 4.1 and 4.1.7):

- If exact ties are considered impossible, a lower bound of 0 means the arm is the best, and an upper bound of 0 means it is not the best.
- If ties are possible, a lower bound of 0 means the arm is one of the best, and an upper bound of 0 means it is at most tied for the best.

## Not implemented

- Smaller-is-better analyses. Multiply the response by $-1$ instead.
- Unconstrained MCB intervals deduced from Tukey-Kramer all-pairwise comparisons (Hsu 1996, Section 4.2.4.1).
- All-pairwise (Tukey, Tukey-Kramer) and comparison-with-a-control (Dunnett) procedures as analyses in their own right.
- Nonparametric MCB, and MCB in the general linear model.
- Sample-size computation for unbalanced designs.
