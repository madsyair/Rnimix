# Potts prior kernels for NIMBLE

The Potts prior on a labelling, \\p(z \mid \beta) \propto \exp(\beta
\sum\_{i \sim j} I(z_i = z_j))\\ (Potts, 1952; Besag, 1974). The
normalising constant is intractable but depends only on \\\beta\\, so
with \\\beta\\ fixed it cancels and the unnormalised density gives exact
MCMC for the labels; when \\\beta\\ is estimated, a pseudo-likelihood
sampler is used instead (see
[`MRFEngine`](https://madsyair.github.io/Rnimix/reference/MRFEngine-class.md)).

## Usage

``` r
dPottsNimix(x, beta, e1, e2, log = 0)

rPottsNimix(n, beta, e1, e2)
```

## Arguments

- x:

  Vector of labels.

- beta:

  Interaction strength.

- e1, e2:

  Edge endpoints of the neighbourhood graph.

- log:

  Return the log density.

- n:

  Number of draws (required by NIMBLE).

## Value

The (log) unnormalised Potts density, or a stub draw.

## Details

The paired `r` function is a stub: labels are always supplied as initial
values and updated by the Gibbs sweep, so exact Potts simulation is
never needed for inference. It exists because `registerDistributions`
requires one.

Exported because NIMBLE's code generator resolves distribution names on
the search path; not intended to be called directly.

## References

Besag, J. (1974). Spatial interaction and the statistical analysis of
lattice systems. *JRSS B*, 36(2), 192–236.
[doi:10.1111/j.2517-6161.1974.tb00999.x](https://doi.org/10.1111/j.2517-6161.1974.tb00999.x)
