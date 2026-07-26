# Forward-algorithm kernels for hidden Markov mixtures

Custom NIMBLE distributions implementing the scaled forward algorithm
(Rabiner, 1989) for a hidden Markov mixture: the discrete state path is
marginalised out of the likelihood analytically, leaving a density in
the component parameters and the transition matrix. The paired `r`
function exists because `registerDistributions` requires one; the
marginalised model never simulates from it during MCMC.

## Usage

``` r
dRegimeHMMNorm_k(x, mu, s2, P, init, len, log = 0)

rRegimeHMMNorm_k(n, mu, s2, P, init, len)

dRegimeHMMBinom_k(x, prob, size, P, init, len, log = 0)

rRegimeHMMBinom_k(n, prob, size, P, init, len)

dRegimeHMMBinomReg_k(x, X, beta, size, P, init, log = 0)

rRegimeHMMBinomReg_k(n, X, beta, size, P, init)

dRegimeHMMFOSSEP_k(x, mu, sigma, alpha, theta, P, init, len, log = 0)

rRegimeHMMFOSSEP_k(n, mu, sigma, alpha, theta, P, init, len)

dRegimeHMMFSSN_k(x, mu, sigma, alpha, P, init, len, log = 0)

rRegimeHMMFSSN_k(n, mu, sigma, alpha, P, init, len)

dRegimeHMMFSST_k(x, mu, sigma, alpha, nu, P, init, len, log = 0)

rRegimeHMMFSST_k(n, mu, sigma, alpha, nu, P, init, len)

dRegimeHMMGMSNB_k(x, mu, sigma, alpha, theta, P, init, len, log = 0)

rRegimeHMMGMSNB_k(n, mu, sigma, alpha, theta, P, init, len)

dRegimeHMMNormReg_k(x, X, beta, s2, P, init, log = 0)

rRegimeHMMNormReg_k(n, X, beta, s2, P, init)

dRegimeHMMJFST_k(x, mu, sigma, alpha, theta, P, init, len, log = 0)

rRegimeHMMJFST_k(n, mu, sigma, alpha, theta, P, init, len)

dRegimeHMMLEP_k(x, mu, sigma, nu, P, init, len, log = 0)

rRegimeHMMLEP_k(n, mu, sigma, nu, P, init, len)

dRegimeHMMMSNB_k(x, mu, sigma, alpha, P, init, len, log = 0)

rRegimeHMMMSNB_k(n, mu, sigma, alpha, P, init, len)

dRegimeHMMMSNB2a_k(x, mu, sigma, alpha, P, init, len, log = 0)

rRegimeHMMMSNB2a_k(n, mu, sigma, alpha, P, init, len)

dRegimeHMMPois_k(x, lambda, P, init, len, log = 0)

rRegimeHMMPois_k(n, lambda, P, init, len)

dRegimeHMMPoisReg_k(x, X, beta, P, init, log = 0)

rRegimeHMMPoisReg_k(n, X, beta, P, init)

dRegimeHMMSEP_k(x, mu, sigma, nu, P, init, len, log = 0)

rRegimeHMMSEP_k(n, mu, sigma, nu, P, init, len)

dRegimeHMMT_k(x, mu, tau, df, P, init, len, log = 0)

rRegimeHMMT_k(n, mu, tau, df, P, init, len)

dRegimeHMMStudentTReg_k(x, X, beta, s2, df, P, init, log = 0)

rRegimeHMMStudentTReg_k(n, X, beta, s2, df, P, init)

dRegimeHMMMSNBReg_k(x, X, beta, sigma, alpha, P, init, log = 0)

rRegimeHMMMSNBReg_k(n, X, beta, sigma, alpha, P, init)

dRegimeHMMMSNB2aReg_k(x, X, beta, sigma, alpha, P, init, log = 0)

rRegimeHMMMSNB2aReg_k(n, X, beta, sigma, alpha, P, init)

dRegimeHMMFSSNReg_k(x, X, beta, sigma, alpha, P, init, log = 0)

rRegimeHMMFSSNReg_k(n, X, beta, sigma, alpha, P, init)

dRegimeHMMSEPReg_k(x, X, beta, sigma, nu, P, init, log = 0)

rRegimeHMMSEPReg_k(n, X, beta, sigma, nu, P, init)

dRegimeHMMLEPReg_k(x, X, beta, sigma, nu, P, init, log = 0)

rRegimeHMMLEPReg_k(n, X, beta, sigma, nu, P, init)

dRegimeHMMGMSNBReg_k(x, X, beta, sigma, alpha, theta, P, init, log = 0)

rRegimeHMMGMSNBReg_k(n, X, beta, sigma, alpha, theta, P, init)

dRegimeHMMFSSTReg_k(x, X, beta, sigma, alpha, nu, P, init, log = 0)

rRegimeHMMFSSTReg_k(n, X, beta, sigma, alpha, nu, P, init)

dRegimeHMMFOSSEPReg_k(x, X, beta, sigma, alpha, theta, P, init, log = 0)

rRegimeHMMFOSSEPReg_k(n, X, beta, sigma, alpha, theta, P, init)

dRegimeHMMJFSTReg_k(x, X, beta, sigma, alpha, theta, P, init, log = 0)

rRegimeHMMJFSTReg_k(n, X, beta, sigma, alpha, theta, P, init)
```

## Arguments

- x:

  Observed series (vector).

- mu:

  Per-state location.

- s2:

  Per-state variance.

- P:

  Transition matrix, row-stochastic.

- init:

  Initial state distribution.

- len:

  Length of the series (needed by the `r` function, which cannot infer
  it from the other arguments).

- log:

  Return the log density.

- n:

  Number of draws (required by NIMBLE; always 1).

## Value

`dRegimeHMMNorm_k` returns the (log) marginal likelihood;
`rRegimeHMMNorm_k` returns a simulated series.

## Details

These are exported because NIMBLE resolves user-defined distributions by
name on the search path when generating C++ code. They are not intended
to be called directly.

## References

Rabiner, L.R. (1989). A tutorial on hidden Markov models and selected
applications in speech recognition. *Proceedings of the IEEE*, 77(2),
257–286. [doi:10.1109/5.18626](https://doi.org/10.1109/5.18626)
