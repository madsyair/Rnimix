# Scalar neo-normal densities for NIMBLE

Custom NIMBLE distributions for the neo-normal families, used both as
mixture component densities and as the emission primitives called by the
hidden Markov forward kernels.

## Usage

``` r
dMSNBurr_k(x, mu, sigma, alpha, log = 0)

rMSNBurr_k(n, mu, sigma, alpha)

dMSNBurr2a_k(x, mu, sigma, alpha, log = 0)

rMSNBurr2a_k(n, mu, sigma, alpha)

dGMSNBurr_k(x, mu, sigma, alpha, theta, log = 0)

rGMSNBurr_k(n, mu, sigma, alpha, theta)

dSEP_k(x, mu, sigma, nu, log = 0)

rSEP_k(n, mu, sigma, nu)

dLEP_k(x, mu, sigma, nu, log = 0)

rLEP_k(n, mu, sigma, nu)

dFSSN_k(x, mu, sigma, alpha, log = 0)

rFSSN_k(n, mu, sigma, alpha)

dFOSSEP_k(x, mu, sigma, alpha, theta, log = 0)

rFOSSEP_k(n, mu, sigma, alpha, theta)

dFSST_k(x, mu, sigma, alpha, nu, log = 0)

rFSST_k(n, mu, sigma, alpha, nu)

dJFST_k(x, mu, sigma, alpha, theta, log = 0)

rJFST_k(n, mu, sigma, alpha, theta)

dSkewMvN_k(x, mu, Sigma, gam, log = 0)

rSkewMvN_k(n, mu, Sigma, gam)

dSkewMvNO_k(x, mu, Sigma, gam, theta, log = 0)

rSkewMvNO_k(n, mu, Sigma, gam, theta)

dSkewMvNOG_k(x, mu, Sigma, gam, theta, log = 0)

dSkewMvIT_k(x, mu, Sigma, gam, nu, log = 0)

rSkewMvIT_k(n, mu, Sigma, gam, nu)

dSkewMvITO_k(x, mu, Sigma, gam, nu, theta, log = 0)

rSkewMvITO_k(n, mu, Sigma, gam, nu, theta)

dSkewMvITOG_k(x, mu, Sigma, gam, nu, theta, log = 0)
```

## Arguments

- x:

  Observation.

- mu:

  Location.

- sigma:

  Scale.

- alpha:

  Shape.

- log:

  Return the log density.

- n:

  Number of draws (required by NIMBLE; always 1).

## Value

The (log) density, or a draw.

## Details

`dMSNBurr_k` implements the MSNBurr density (Iriawan, 2000; Choir,
2020), parameterised by location `mu`, scale `sigma` and shape `alpha`,
where `alpha = 1` recovers symmetry and the normalising constant `omega`
is computed on the log scale for stability at extreme `alpha`.

These are exported because NIMBLE resolves user-defined distributions by
name on the search path when generating C++ code. They are not intended
to be called directly.

## References

Iriawan, N. (2000). *Computationally intensive approaches to inference
in neo-normal linear models*. PhD thesis, Curtin University.

Choir, A.S. (2020). *The new neo-normal distributions and their
properties*. Dissertation, Institut Teknologi Sepuluh Nopember.
ISBN:978-602-5539-45-4
