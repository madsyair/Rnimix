# kernels-hmm.R
#
# HMM forward kernels defined at BUILD time and exported.
#
# Why not created lazily and assigned into globalenv(), as the rest of the HMM
# kernels still are? Because writing to the global environment is a CRAN policy
# violation. The lazy pattern existed because NIMBLE has to resolve these
# functions by name at three separate moments -- registerDistributions(),
# nimbleModel(), and the C++ code generation in compileNimble() -- which happen
# in different frames at different times, and globalenv() is the only scope
# that survives across all three.
#
# Defining them at build time solves that, but only if they are also EXPORTED:
# a namespace-private kernel is found when the model is built and then fails
# during compilation, because the code generator resolves names on the search
# path. Exporting them (with nimble attached via Depends) satisfies all three
# lookups with no global assignment anywhere. Verified end to end: build,
# install, fresh session, nimbleModel, compileNimble, calculate.
#
# These are exported for NIMBLE's benefit, not as a user-facing API; they are
# documented together and flagged internal.

#' Forward-algorithm kernels for hidden Markov mixtures
#'
#' Custom NIMBLE distributions implementing the scaled forward algorithm
#' (Rabiner, 1989) for a hidden Markov mixture: the discrete state path is
#' marginalised out of the likelihood analytically, leaving a density in the
#' component parameters and the transition matrix. The paired \code{r} function
#' exists because \code{registerDistributions} requires one; the marginalised
#' model never simulates from it during MCMC.
#'
#' These are exported because NIMBLE resolves user-defined distributions by
#' name on the search path when generating C++ code. They are not intended to
#' be called directly.
#'
#' @param x Observed series (vector).
#' @param n Number of draws (required by NIMBLE; always 1).
#' @param mu Per-state location.
#' @param s2 Per-state variance.
#' @param P Transition matrix, row-stochastic.
#' @param init Initial state distribution.
#' @param len Length of the series (needed by the \code{r} function, which
#'   cannot infer it from the other arguments).
#' @param log Return the log density.
#' @return \code{dRegimeHMMNorm_k} returns the (log) marginal likelihood;
#'   \code{rRegimeHMMNorm_k} returns a simulated series.
#' @references Rabiner, L.R. (1989). A tutorial on hidden Markov models and
#'   selected applications in speech recognition. \emph{Proceedings of the
#'   IEEE}, 77(2), 257--286. \doi{10.1109/5.18626}
#' @keywords internal
#' @name RegimeHMMKernels
NULL

#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMNorm_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), s2 = double(1),
                 P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * dnorm(x[1], mu[s], sqrt(s2[s]))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * dnorm(x[t], mu[sp], sqrt(s2[sp])) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })

#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMNorm_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), s2 = double(1),
                 P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sqrt(s2[z])); z <- rcat(1, P[z, ]) }
    return(out) })
## --- Binomial ------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMBinom_k <- nimble::nimbleFunction(
  run = function(x = double(1), prob = double(1), size = double(0),
                 P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(prob)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dbinom(x[1], size, prob[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dbinom(x[t], size, prob[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMBinom_k <- nimble::nimbleFunction(
  run = function(n = integer(0), prob = double(1), size = double(0),
                 P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rbinom(1, size, prob[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- Binomial regression -------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMBinomReg_k <- nimble::nimbleFunction(
  run = function(x = double(1), X = double(2), beta = double(2),
                 size = double(0), P = double(2), init = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    T <- length(x); S <- dim(beta)[1]; p <- dim(X)[2]
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S)
      a[s] <- init[s] * dbinom(x[1], size,
                               ilogit(sum(X[1, 1:p] * beta[s, 1:p])))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * dbinom(x[t], size,
                               ilogit(sum(X[t, 1:p] * beta[sp, 1:p]))) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMBinomReg_k <- nimble::nimbleFunction(
  run = function(n = integer(0), X = double(2), beta = double(2),
                 size = double(0), P = double(2), init = double(1)) {
    returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rbinom(1, size, 0.5); z <- rcat(1, P[z, ]) }
    return(out) })

## --- FOSSEP --------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMFOSSEP_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 alpha = double(1), theta = double(1),
                 P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dFOSSEP_k(x[1], mu[s], sigma[s], alpha[s], theta[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dFOSSEP_k(x[t], mu[sp], sigma[sp], alpha[sp], theta[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMFOSSEP_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 alpha = double(1), theta = double(1),
                 P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- FSSN ----------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMFSSN_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 alpha = double(1), P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dFSSN_k(x[1], mu[s], sigma[s], alpha[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dFSSN_k(x[t], mu[sp], sigma[sp], alpha[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMFSSN_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 alpha = double(1), P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- FSST ----------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMFSST_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 alpha = double(1), nu = double(1),
                 P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dFSST_k(x[1], mu[s], sigma[s], alpha[s], nu[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dFSST_k(x[t], mu[sp], sigma[sp], alpha[sp], nu[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMFSST_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 alpha = double(1), nu = double(1),
                 P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- GMSNBurr ------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMGMSNB_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 alpha = double(1), theta = double(1),
                 P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dGMSNBurr_k(x[1], mu[s], sigma[s], alpha[s], theta[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dGMSNBurr_k(x[t], mu[sp], sigma[sp], alpha[sp], theta[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMGMSNB_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 alpha = double(1), theta = double(1),
                 P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- Gaussian regression -------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMNormReg_k <- nimble::nimbleFunction(
  run = function(x = double(1), X = double(2), beta = double(2),
                 s2 = double(1), P = double(2), init = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    T <- length(x); S <- length(s2); p <- dim(X)[2]
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S)
      a[s] <- init[s] * dnorm(x[1], sum(X[1, 1:p] * beta[s, 1:p]),
                              sqrt(s2[s]))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * dnorm(x[t], sum(X[t, 1:p] * beta[sp, 1:p]),
                              sqrt(s2[sp])) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMNormReg_k <- nimble::nimbleFunction(
  run = function(n = integer(0), X = double(2), beta = double(2),
                 s2 = double(1), P = double(2), init = double(1)) {
    returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, 0, 1); z <- rcat(1, P[z, ]) }
    return(out) })

## --- JFST ----------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMJFST_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 alpha = double(1), theta = double(1),
                 P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dJFST_k(x[1], mu[s], sigma[s], alpha[s], theta[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dJFST_k(x[t], mu[sp], sigma[sp], alpha[sp], theta[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMJFST_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 alpha = double(1), theta = double(1),
                 P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- LEP -----------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMLEP_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 nu = double(1), P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dLEP_k(x[1], mu[s], sigma[s], nu[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dLEP_k(x[t], mu[sp], sigma[sp], nu[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMLEP_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 nu = double(1), P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- MSNBurr -------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMMSNB_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 alpha = double(1), P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dMSNBurr_k(x[1], mu[s], sigma[s], alpha[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dMSNBurr_k(x[t], mu[sp], sigma[sp], alpha[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMMSNB_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 alpha = double(1), P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- MSNBurr2a -----------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMMSNB2a_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 alpha = double(1), P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dMSNBurr2a_k(x[1], mu[s], sigma[s], alpha[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dMSNBurr2a_k(x[t], mu[sp], sigma[sp], alpha[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMMSNB2a_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 alpha = double(1), P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- Poisson -------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMPois_k <- nimble::nimbleFunction(
  run = function(x = double(1), lambda = double(1),
                 P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(lambda)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * dpois(x[1], lambda[s])
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * dpois(x[t], lambda[sp]) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMPois_k <- nimble::nimbleFunction(
  run = function(n = integer(0), lambda = double(1),
                 P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rpois(1, lambda[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- Poisson regression --------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMPoisReg_k <- nimble::nimbleFunction(
  run = function(x = double(1), X = double(2), beta = double(2),
                 P = double(2), init = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    T <- length(x); S <- dim(beta)[1]; p <- dim(X)[2]
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S)
      a[s] <- init[s] * dpois(x[1], exp(sum(X[1, 1:p] * beta[s, 1:p])))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * dpois(x[t], exp(sum(X[t, 1:p] * beta[sp, 1:p]))) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMPoisReg_k <- nimble::nimbleFunction(
  run = function(n = integer(0), X = double(2), beta = double(2),
                 P = double(2), init = double(1)) {
    returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rpois(1, 1); z <- rcat(1, P[z, ]) }
    return(out) })

## --- SEP -----------------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMSEP_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), sigma = double(1),
                 nu = double(1), P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * exp(dSEP_k(x[1], mu[s], sigma[s], nu[s], 1))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * exp(dSEP_k(x[t], mu[sp], sigma[sp], nu[sp], 1)) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMSEP_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), sigma = double(1),
                 nu = double(1), P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, mu[z], sigma[z]); z <- rcat(1, P[z, ]) }
    return(out) })

## --- Student-t -----------------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMT_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), tau = double(1),
                 df = double(0), P = double(2), init = double(1),
                 len = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0)); T <- length(x); S <- length(mu)
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S) a[s] <- init[s] * dt_nonstandard(x[1], df, mu[s], 1/sqrt(tau[s]))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * dt_nonstandard(x[t], df, mu[sp], 1/sqrt(tau[sp])) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMT_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), tau = double(1),
                 df = double(0), P = double(2), init = double(1), len = double(0)) {
    returnType(double(1)); Tlen <- len; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rt_nonstandard(1, df, mu[z], 1/sqrt(tau[z])); z <- rcat(1, P[z, ]) }
    return(out) })

## --- Student-t regression ------------------------------------------------
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMStudentTReg_k <- nimble::nimbleFunction(
  run = function(x = double(1), X = double(2), beta = double(2),
                 s2 = double(1), df = double(0), P = double(2),
                 init = double(1), log = integer(0, default = 0)) {
    returnType(double(0))
    T <- length(x); S <- length(s2); p <- dim(X)[2]
    a <- numeric(S); an <- numeric(S); ll <- 0
    for (s in 1:S)
      a[s] <- init[s] * dt_nonstandard(x[1], df,
                 sum(X[1, 1:p] * beta[s, 1:p]), sqrt(s2[s]))
    c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
    ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
    if (T >= 2) for (t in 2:T) {
      for (sp in 1:S) { acc <- 0
        for (s in 1:S) acc <- acc + a[s] * P[s, sp]
        an[sp] <- acc * dt_nonstandard(x[t], df,
                    sum(X[t, 1:p] * beta[sp, 1:p]), sqrt(s2[sp])) }
      ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
    if (log) return(ll) else return(exp(ll))
  })
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMStudentTReg_k <- nimble::nimbleFunction(
  run = function(n = integer(0), X = double(2), beta = double(2),
                 s2 = double(1), df = double(0), P = double(2),
                 init = double(1)) {
    returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
    z <- rcat(1, init)
    for (t in 1:Tlen) { out[t] <- rnorm(1, 0, 1); z <- rcat(1, P[z, ]) }
    return(out) })

# --- Neo-normal Markov-switching regression kernels -------------------------
#
# These nine families differ only in which primitive density they call and how
# many shape parameters they carry, and the shape count changes the run()
# signature -- so one static kernel cannot cover them. The generator below
# writes each source and evaluates it; running it here, at build time, keeps
# the results out of the global environment. Each pair is exported for the
# same reason as the rest: NIMBLE's code generator resolves names on the
# search path.


.neoHMMRegKernelPair <- function(densName, shapeNames) {
  shapeSig <- paste(sprintf("%s = double(1)", shapeNames), collapse = ", ")
  shapeArg1 <- paste(sprintf("%s[s]", shapeNames), collapse = ", ")
  shapeArgT <- paste(sprintf("%s[sp]", shapeNames), collapse = ", ")
  firstShape <- shapeNames[1]
  dSrc <- sprintf('
    nimble::nimbleFunction(run = function(x = double(1), X = double(2),
        beta = double(2), %s, P = double(2), init = double(1),
        log = integer(0, default = 0)) {
      returnType(double(0))
      T <- length(x); S <- length(%s); p <- dim(X)[2]
      a <- numeric(S); an <- numeric(S); ll <- 0
      for (s in 1:S)
        a[s] <- init[s] * exp(%s(x[1],
                  sum(X[1, 1:p] * beta[s, 1:p]), %s, 1))
      c1 <- sum(a); if (c1 <= 0) { if (log) return(-Inf) else return(0) }
      ll <- ll + log(c1); for (s in 1:S) a[s] <- a[s] / c1
      if (T >= 2) for (t in 2:T) {
        for (sp in 1:S) { acc <- 0
          for (s in 1:S) acc <- acc + a[s] * P[s, sp]
          an[sp] <- acc * exp(%s(x[t],
                      sum(X[t, 1:p] * beta[sp, 1:p]), %s, 1)) }
        ct <- sum(an); if (ct <= 0) { if (log) return(-Inf) else return(0) }
        ll <- ll + log(ct); for (s in 1:S) a[s] <- an[s] / ct }
      if (log) return(ll) else return(exp(ll))
    })',
    shapeSig, firstShape, densName, shapeArg1, densName, shapeArgT)
  dK <- eval(parse(text = dSrc))
  rSrc <- sprintf('
    nimble::nimbleFunction(run = function(n = integer(0), X = double(2),
        beta = double(2), %s, P = double(2), init = double(1)) {
      returnType(double(1)); Tlen <- dim(X)[1]; out <- numeric(Tlen)
      z <- rcat(1, init)
      for (t in 1:Tlen) { out[t] <- rnorm(1, 0, 1); z <- rcat(1, P[z, ]) }
      return(out) })', shapeSig)
  rK <- eval(parse(text = rSrc))
  list(d = dK, r = rK)
}


.kMSNBReg <- .neoHMMRegKernelPair("dMSNBurr_k", c("sigma", "alpha"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMMSNBReg_k <- .kMSNBReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMMSNBReg_k <- .kMSNBReg$r

.kMSNB2aReg <- .neoHMMRegKernelPair("dMSNBurr2a_k", c("sigma", "alpha"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMMSNB2aReg_k <- .kMSNB2aReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMMSNB2aReg_k <- .kMSNB2aReg$r

.kFSSNReg <- .neoHMMRegKernelPair("dFSSN_k", c("sigma", "alpha"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMFSSNReg_k <- .kFSSNReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMFSSNReg_k <- .kFSSNReg$r

.kSEPReg <- .neoHMMRegKernelPair("dSEP_k", c("sigma", "nu"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMSEPReg_k <- .kSEPReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMSEPReg_k <- .kSEPReg$r

.kLEPReg <- .neoHMMRegKernelPair("dLEP_k", c("sigma", "nu"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMLEPReg_k <- .kLEPReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMLEPReg_k <- .kLEPReg$r

.kGMSNBReg <- .neoHMMRegKernelPair("dGMSNBurr_k", c("sigma", "alpha", "theta"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMGMSNBReg_k <- .kGMSNBReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMGMSNBReg_k <- .kGMSNBReg$r

.kFSSTReg <- .neoHMMRegKernelPair("dFSST_k", c("sigma", "alpha", "nu"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMFSSTReg_k <- .kFSSTReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMFSSTReg_k <- .kFSSTReg$r

.kFOSSEPReg <- .neoHMMRegKernelPair("dFOSSEP_k", c("sigma", "alpha", "theta"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMFOSSEPReg_k <- .kFOSSEPReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMFOSSEPReg_k <- .kFOSSEPReg$r

.kJFSTReg <- .neoHMMRegKernelPair("dJFST_k", c("sigma", "alpha", "theta"))
#' @rdname RegimeHMMKernels
#' @export
dRegimeHMMJFSTReg_k <- .kJFSTReg$d
#' @rdname RegimeHMMKernels
#' @export
rRegimeHMMJFSTReg_k <- .kJFSTReg$r
