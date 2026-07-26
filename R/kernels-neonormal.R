# kernels-neonormal.R
#
# Scalar neo-normal densities defined at BUILD time and exported.
#
# These are the primitives the HMM forward kernels call: dRegimeHMMMSNB_k
# evaluates dMSNBurr_k once per observation per state. They used to be built
# inside .nimixEnsureMSNBurr() with their enclosure forced to globalenv(),
# because a namespace-private enclosure fails NIMBLE's C++ code generation for
# these scalar densities -- writing to the global environment is a CRAN policy
# violation, so that could not stay.
#
# Exporting them resolves it. The code generator looks names up on the search
# path, so a build-time definition that is exported satisfies it, while a
# build-time definition that is merely in the namespace does not. Verified for
# exactly the nesting used here: an exported scalar density called from inside
# an exported vector kernel compiles and calculates correctly.

#' Scalar neo-normal densities for NIMBLE
#'
#' Custom NIMBLE distributions for the neo-normal families, used both as
#' mixture component densities and as the emission primitives called by the
#' hidden Markov forward kernels.
#'
#' \code{dMSNBurr_k} implements the MSNBurr density (Iriawan, 2000; Choir,
#' 2020), parameterised by location \code{mu}, scale \code{sigma} and shape
#' \code{alpha}, where \code{alpha = 1} recovers symmetry and the normalising
#' constant \code{omega} is computed on the log scale for stability at extreme
#' \code{alpha}.
#'
#' These are exported because NIMBLE resolves user-defined distributions by
#' name on the search path when generating C++ code. They are not intended to
#' be called directly.
#'
#' @param x Observation.
#' @param n Number of draws (required by NIMBLE; always 1).
#' @param mu Location.
#' @param sigma Scale.
#' @param alpha Shape.
#' @param log Return the log density.
#' @return The (log) density, or a draw.
#' @references Iriawan, N. (2000). \emph{Computationally intensive approaches
#'   to inference in neo-normal linear models}. PhD thesis, Curtin University.
#'
#'   Choir, A.S. (2020). \emph{The new neo-normal distributions and their
#'   properties}. Dissertation, Institut Teknologi Sepuluh Nopember.
#'   ISBN:978-602-5539-45-4
#' @keywords internal
#' @name NeoNormalKernels
NULL

#' @rdname NeoNormalKernels
#' @export
dMSNBurr_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 alpha = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if (alpha < 1e-300) {
      lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
    } else {
      lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
    }
    omega <- exp(lomega); zo <- -omega * ((x - mu) / sigma)
    u <- zo - log(alpha); sp <- max(u, 0) + log1p(exp(-abs(u)))
    lp <- lomega - log(sigma) + zo - (alpha + 1) * sp
    if (log) return(lp) else return(exp(lp))
  })

#' @rdname NeoNormalKernels
#' @export
rMSNBurr_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 alpha = double(0)) {
    returnType(double(0))
    if (alpha < 1e-300) {
      lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
    } else {
      lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
    }
    omega <- exp(lomega); p <- runif(1)
    lt <- log(exp(-log(p) / alpha) - 1)
    return(mu - (sigma / omega) * (log(alpha) + lt))
  })
## --- MSNBurr2a ---------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dMSNBurr2a_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 alpha = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if (alpha < 1e-300) {
      lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
    } else {
      lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
    }
    omega <- exp(lomega); zt <- omega * ((x - mu) / sigma)
    u <- log(alpha) - zt; sp <- max(u, 0) + log1p(exp(-abs(u)))
    lp <- lomega - log(sigma) + (alpha + 1) * log(alpha) - alpha * zt -
      (alpha + 1) * sp
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rMSNBurr2a_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 alpha = double(0)) {
    returnType(double(0))
    if (alpha < 1e-300) {
      lomega <- -(alpha + 1) * log(alpha) - 0.9189385332046727
    } else {
      lomega <- (alpha + 1) * log1p(1 / alpha) - 0.9189385332046727
    }
    omega <- exp(lomega); q <- 1 - runif(1)
    lt <- log(exp(-log(q) / alpha) - 1)
    return(mu + (sigma / omega) * (log(alpha) + lt))
  })

## --- GMSNBurr ----------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dGMSNBurr_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 alpha = double(0), theta = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    lbeta_at <- lgamma(alpha) + lgamma(theta) - lgamma(alpha + theta)
    if (alpha < 1e-300) {
      lomega <- -0.9189385332046727 + lbeta_at +
        alpha * (log(theta) - log(alpha))
    } else {
      lomega <- -0.9189385332046727 + lbeta_at -
        theta * (log(theta) - log(alpha)) +
        (alpha + theta) * log1p(theta / alpha)
    }
    omega <- exp(lomega)
    zo <- -omega * ((x - mu) / sigma)
    zoa <- zo + log(theta) - log(alpha)
    sp <- max(zoa, 0) + log1p(exp(-abs(zoa)))
    lp <- lomega - log(sigma) + theta * (log(theta) - log(alpha)) +
      theta * zo - (alpha + theta) * sp - lbeta_at
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rGMSNBurr_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 alpha = double(0), theta = double(0)) {
    returnType(double(0))
    lbeta_at <- lgamma(alpha) + lgamma(theta) - lgamma(alpha + theta)
    if (alpha < 1e-300) {
      lomega <- -0.9189385332046727 + lbeta_at +
        alpha * (log(theta) - log(alpha))
    } else {
      lomega <- -0.9189385332046727 + lbeta_at -
        theta * (log(theta) - log(alpha)) +
        (alpha + theta) * log1p(theta / alpha)
    }
    omega <- exp(lomega)
    Xg <- rgamma(1, shape = alpha, rate = 1)
    Yg <- rgamma(1, shape = theta, rate = 1)
    log_ratio <- log(Yg) - log(Xg) + log(alpha) - log(theta)
    return(mu - (sigma / omega) * log_ratio)
  })

## --- SEP ---------------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dSEP_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 nu = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    z <- abs(x - mu) / sigma
    lp <- -log(2) - (1 / nu) * log(2) - lgamma(1 + 1 / nu) -
      log(sigma) - 0.5 * z^nu
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rSEP_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 nu = double(0)) {
    returnType(double(0))
    W <- rgamma(1, shape = 1 / nu, rate = 0.5)
    za <- W^(1 / nu)
    s <- 2 * (runif(1, 0, 1) < 0.5) - 1
    return(mu + sigma * s * za)
  })

## --- LEP ---------------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dLEP_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 nu = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    z <- abs(x - mu) / sigma
    lp <- -log(2) - (1 / nu) * log(nu) - lgamma(1 + 1 / nu) -
      log(sigma) - z^nu / nu
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rLEP_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 nu = double(0)) {
    returnType(double(0))
    W <- rgamma(1, shape = 1 / nu, rate = 1 / nu)
    za <- W^(1 / nu)
    s <- 2 * (runif(1, 0, 1) < 0.5) - 1
    return(mu + sigma * s * za)
  })

## --- FSSN --------------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dFSSN_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 alpha = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    z <- (x - mu) / sigma
    lc <- log(2) - log(sigma) - log(alpha + 1 / alpha) - 0.9189385332046727
    # FS convention: alpha == gamma, so alpha > 1 skews right.
    scale <- 1 / alpha
    if (z < 0) scale <- alpha
    lp <- lc - 0.5 * (z * scale)^2
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rFSSN_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 alpha = double(0)) {
    returnType(double(0))
    a2 <- alpha * alpha
    W <- abs(rnorm(1, 0, 1))
    if (runif(1, 0, 1) < a2 / (1 + a2)) {
      z <- W * alpha        # positive side, sd alpha
    } else {
      z <- -W / alpha       # negative side, sd 1/alpha
    }
    return(mu + sigma * z)
  })

## --- FOSSEP ------------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dFOSSEP_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 alpha = double(0), theta = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    z <- (x - mu) / sigma
    az <- abs(z)
    if (z < 0) {
      base <- -0.5 * (alpha * az)^theta
    } else {
      base <- -0.5 * (az / alpha)^theta
    }
    lp <- base - log(sigma) + log(alpha) - log1p(alpha^2) -
      (1 / theta) * log(2) - lgamma(1 + 1 / theta)
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rFOSSEP_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 alpha = double(0), theta = double(0)) {
    returnType(double(0))
    a2 <- alpha * alpha
    W <- rgamma(1, shape = 1 / theta, rate = 1)
    mag <- (2 * W)^(1 / theta)
    if (runif(1, 0, 1) < a2 / (1 + a2)) {
      z <- alpha * mag
    } else {
      z <- -mag / alpha
    }
    return(mu + sigma * z)
  })

## --- FSST --------------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dFSST_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 alpha = double(0), nu = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    z <- (x - mu) / sigma
    # FS convention: alpha == gamma, so alpha > 1 skews right.
    tscale <- 1 / alpha
    if (z < 0) tscale <- alpha
    t <- z * tscale
    tlp <- lgamma((nu + 1) / 2) - lgamma(nu / 2) -
      0.5 * (log(nu) + 1.1447298858494002) -
      ((nu + 1) / 2) * log1p(t * t / nu)
    lp <- log(2) - log(sigma) - log(alpha + 1 / alpha) + tlp
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rFSST_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 alpha = double(0), nu = double(0)) {
    returnType(double(0))
    a2 <- alpha * alpha
    Zt <- rnorm(1, 0, 1)
    Wt <- rgamma(1, shape = nu / 2, rate = 0.5)
    Tt <- abs(Zt * sqrt(nu / Wt))
    if (runif(1, 0, 1) < a2 / (1 + a2)) {
      z <- Tt * alpha
    } else {
      z <- -Tt / alpha
    }
    return(mu + sigma * z)
  })

## --- JFST --------------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dJFST_k <- nimble::nimbleFunction(
  run = function(x = double(0), mu = double(0), sigma = double(0),
                 alpha = double(0), theta = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    z <- (x - mu) / sigma
    rz <- z / sqrt(alpha + theta + z * z)
    lbeta_at <- lgamma(alpha) + lgamma(theta) - lgamma(alpha + theta)
    lp <- (alpha + 0.5) * log1p(rz) + (theta + 0.5) * log1p(-rz) -
      (alpha + theta - 1) * log(2) - 0.5 * log(alpha + theta) -
      lbeta_at - log(sigma)
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rJFST_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0),
                 alpha = double(0), theta = double(0)) {
    returnType(double(0))
    B <- rbeta(1, alpha, theta)
    rz <- 2 * B - 1
    z <- rz * sqrt((alpha + theta) / (1 - rz * rz))
    return(mu + sigma * z)
  })

## --- SkewMvN -----------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dSkewMvN_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), Sigma = double(2),
                 gam = double(1), log = integer(0, default = 0)) {
    returnType(double(0))
    m <- length(x)
    U <- chol(Sigma)
    Ui <- inverse(U)
    d <- x - mu
    lp <- 0
    for (j in 1:m) {
      s <- 0
      for (i in 1:m) s <- s + d[i] * Ui[i, j]
      sc <- 1 / gam[j]
      if (s < 0) sc <- gam[j]
      lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) - 0.9189385332046727 -
        0.5 * (s * sc)^2 - log(U[j, j])
    }
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rSkewMvN_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), Sigma = double(2),
                 gam = double(1)) {
    returnType(double(1))
    m <- length(mu)
    U <- chol(Sigma)
    eps <- numeric(m)
    for (j in 1:m) {
      a2 <- gam[j] * gam[j]
      W <- abs(rnorm(1, 0, 1))
      if (runif(1, 0, 1) < a2 / (1 + a2)) eps[j] <- W * gam[j]
      else eps[j] <- -W / gam[j]
    }
    out <- numeric(m)
    for (c in 1:m) {
      s <- 0
      for (r in 1:m) s <- s + U[r, c] * eps[r]
      out[c] <- mu[c] + s
    }
    return(out)
  })

## --- SkewMvNO ----------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dSkewMvNO_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), Sigma = double(2),
                 gam = double(1), theta = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    m <- length(x)
    U <- chol(Sigma)
    Ui <- inverse(U)
    dv <- x - mu
    w <- numeric(m)
    for (j in 1:m) {
      s <- 0
      for (i in 1:m) s <- s + dv[i] * Ui[i, j]
      w[j] <- s
    }
    s1 <- sin(theta)
    c1 <- cos(theta)
    vw <- s1 * w[1] + c1 * w[2]
    eps <- numeric(2)
    eps[1] <- w[1] - 2 * s1 * vw
    eps[2] <- w[2] - 2 * c1 * vw
    lp <- 0
    for (j in 1:m) {
      sc <- 1 / gam[j]
      if (eps[j] < 0) sc <- gam[j]
      lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) - 0.9189385332046727 -
        0.5 * (eps[j] * sc)^2 - log(U[j, j])
    }
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rSkewMvNO_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), Sigma = double(2),
                 gam = double(1), theta = double(0)) {
    returnType(double(1))
    m <- length(mu)
    U <- chol(Sigma)
    eps <- numeric(m)
    for (j in 1:m) {
      a2 <- gam[j] * gam[j]
      W <- abs(rnorm(1, 0, 1))
      if (runif(1, 0, 1) < a2 / (1 + a2)) eps[j] <- W * gam[j]
      else eps[j] <- -W / gam[j]
    }
    # eta = A' eps + mu = U' O' eps + mu ; O symmetric so O' = O
    s1 <- sin(theta)
    c1 <- cos(theta)
    ve <- s1 * eps[1] + c1 * eps[2]
    oe <- numeric(2)
    oe[1] <- eps[1] - 2 * s1 * ve
    oe[2] <- eps[2] - 2 * c1 * ve
    out <- numeric(m)
    for (c in 1:m) {
      s <- 0
      for (r in 1:m) s <- s + U[r, c] * oe[r]
      out[c] <- mu[c] + s
    }
    return(out)
  })

## --- SkewMvNOG ---------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dSkewMvNOG_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), Sigma = double(2),
                 gam = double(1), theta = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    m <- length(x)
    U <- chol(Sigma)
    Ui <- inverse(U)
    dv <- x - mu
    w <- numeric(m)
    for (j in 1:m) {
      s <- 0
      for (i in 1:m) s <- s + dv[i] * Ui[i, j]
      w[j] <- s
    }
    eps <- numeric(m)
    for (i in 1:m) eps[i] <- w[i]
    for (j in 2:m) {
      start <- 1 + (j - 2) * (j - 1) / 2
      v <- numeric(j)
      cp <- 1
      v[1] <- sin(theta[start])
      if (j > 2) {
        for (i in 2:(j - 1)) {
          cp <- cp * cos(theta[start + i - 2])
          v[i] <- cp * sin(theta[start + i - 1])
        }
        cp <- cp * cos(theta[start + j - 2])
      } else {
        cp <- cos(theta[start])
      }
      v[j] <- cp
      off <- m - j
      ve <- 0
      for (i in 1:j) ve <- ve + v[i] * eps[off + i]
      vv <- 0
      for (i in 1:j) vv <- vv + v[i] * v[i]
      for (i in 1:j) eps[off + i] <- eps[off + i] - 2 * v[i] * ve / vv
    }
    lp <- 0
    for (j in 1:m) {
      sc <- 1 / gam[j]
      if (eps[j] < 0) sc <- gam[j]
      lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) - 0.9189385332046727 -
        0.5 * (eps[j] * sc)^2 - log(U[j, j])
    }
    if (log) return(lp) else return(exp(lp))
  })

## --- SkewMvIT ----------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dSkewMvIT_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), Sigma = double(2),
                 gam = double(1), nu = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    m <- length(x)
    U <- chol(Sigma)
    Ui <- inverse(U)
    d <- x - mu
    lp <- 0
    for (j in 1:m) {
      s <- 0
      for (i in 1:m) s <- s + d[i] * Ui[i, j]
      sc <- 1 / gam[j]
      if (s < 0) sc <- gam[j]
      t <- s * sc
      tlp <- lgamma((nu[j] + 1) / 2) - lgamma(nu[j] / 2) -
        0.5 * (log(nu[j]) + 1.1447298858494002) -
        ((nu[j] + 1) / 2) * log1p(t * t / nu[j])
      lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) + tlp - log(U[j, j])
    }
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rSkewMvIT_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), Sigma = double(2),
                 gam = double(1), nu = double(1)) {
    returnType(double(1))
    m <- length(mu)
    U <- chol(Sigma)
    eps <- numeric(m)
    for (j in 1:m) {
      a2 <- gam[j] * gam[j]
      Zt <- rnorm(1, 0, 1)
      Wt <- rgamma(1, shape = nu[j] / 2, rate = 0.5)
      Tt <- abs(Zt * sqrt(nu[j] / Wt))
      if (runif(1, 0, 1) < a2 / (1 + a2)) eps[j] <- Tt * gam[j]
      else eps[j] <- -Tt / gam[j]
    }
    out <- numeric(m)
    for (c in 1:m) {
      s <- 0
      for (r in 1:m) s <- s + U[r, c] * eps[r]
      out[c] <- mu[c] + s
    }
    return(out)
  })

## --- SkewMvITO ---------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dSkewMvITO_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), Sigma = double(2),
                 gam = double(1), nu = double(1), theta = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    m <- length(x)
    U <- chol(Sigma)
    Ui <- inverse(U)
    dv <- x - mu
    w <- numeric(m)
    for (j in 1:m) {
      s <- 0
      for (i in 1:m) s <- s + dv[i] * Ui[i, j]
      w[j] <- s
    }
    s1 <- sin(theta)
    c1 <- cos(theta)
    vw <- s1 * w[1] + c1 * w[2]
    eps <- numeric(2)
    eps[1] <- w[1] - 2 * s1 * vw
    eps[2] <- w[2] - 2 * c1 * vw
    lp <- 0
    for (j in 1:m) {
      sc <- 1 / gam[j]
      if (eps[j] < 0) sc <- gam[j]
      t <- eps[j] * sc
      tlp <- lgamma((nu[j] + 1) / 2) - lgamma(nu[j] / 2) -
        0.5 * (log(nu[j]) + 1.1447298858494002) -
        ((nu[j] + 1) / 2) * log1p(t * t / nu[j])
      lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) + tlp - log(U[j, j])
    }
    if (log) return(lp) else return(exp(lp))
  })
#' @rdname NeoNormalKernels
#' @export
rSkewMvITO_k <- nimble::nimbleFunction(
  run = function(n = integer(0), mu = double(1), Sigma = double(2),
                 gam = double(1), nu = double(1), theta = double(0)) {
    returnType(double(1))
    m <- length(mu)
    U <- chol(Sigma)
    eps <- numeric(m)
    for (j in 1:m) {
      a2 <- gam[j] * gam[j]
      Zt <- rnorm(1, 0, 1)
      Wt <- rgamma(1, shape = nu[j] / 2, rate = 0.5)
      Tt <- abs(Zt * sqrt(nu[j] / Wt))
      if (runif(1, 0, 1) < a2 / (1 + a2)) eps[j] <- Tt * gam[j]
      else eps[j] <- -Tt / gam[j]
    }
    s1 <- sin(theta)
    c1 <- cos(theta)
    ve <- s1 * eps[1] + c1 * eps[2]
    oe <- numeric(2)
    oe[1] <- eps[1] - 2 * s1 * ve
    oe[2] <- eps[2] - 2 * c1 * ve
    out <- numeric(m)
    for (c in 1:m) {
      s <- 0
      for (r in 1:m) s <- s + U[r, c] * oe[r]
      out[c] <- mu[c] + s
    }
    return(out)
  })

## --- SkewMvITOG --------------------------------------------------------
#' @rdname NeoNormalKernels
#' @export
dSkewMvITOG_k <- nimble::nimbleFunction(
  run = function(x = double(1), mu = double(1), Sigma = double(2),
                 gam = double(1), nu = double(1), theta = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    m <- length(x)
    U <- chol(Sigma)
    Ui <- inverse(U)
    dv <- x - mu
    w <- numeric(m)
    for (j in 1:m) {
      s <- 0
      for (i in 1:m) s <- s + dv[i] * Ui[i, j]
      w[j] <- s
    }
    eps <- numeric(m)
    for (i in 1:m) eps[i] <- w[i]
    for (j in 2:m) {
      start <- 1 + (j - 2) * (j - 1) / 2
      v <- numeric(j)
      cp <- 1
      v[1] <- sin(theta[start])
      if (j > 2) {
        for (i in 2:(j - 1)) {
          cp <- cp * cos(theta[start + i - 2])
          v[i] <- cp * sin(theta[start + i - 1])
        }
        cp <- cp * cos(theta[start + j - 2])
      } else {
        cp <- cos(theta[start])
      }
      v[j] <- cp
      off <- m - j
      ve <- 0
      for (i in 1:j) ve <- ve + v[i] * eps[off + i]
      vv <- 0
      for (i in 1:j) vv <- vv + v[i] * v[i]
      for (i in 1:j) eps[off + i] <- eps[off + i] - 2 * v[i] * ve / vv
    }
    lp <- 0
    for (j in 1:m) {
      sc <- 1 / gam[j]
      if (eps[j] < 0) sc <- gam[j]
      t <- eps[j] * sc
      tlp <- lgamma((nu[j] + 1) / 2) - lgamma(nu[j] / 2) -
        0.5 * (log(nu[j]) + 1.1447298858494002) -
        ((nu[j] + 1) / 2) * log1p(t * t / nu[j])
      lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) + tlp - log(U[j, j])
    }
    if (log) return(lp) else return(exp(lp))
  })

## --- Potts (spatial Markov random field prior) -----------------------------

#' Potts prior kernels for NIMBLE
#'
#' The Potts prior on a labelling, \eqn{p(z \mid \beta) \propto \exp(\beta
#' \sum_{i \sim j} I(z_i = z_j))} (Potts, 1952; Besag, 1974). The normalising
#' constant is intractable but depends only on \eqn{\beta}, so with \eqn{\beta}
#' fixed it cancels and the unnormalised density gives exact MCMC for the
#' labels; when \eqn{\beta} is estimated, a pseudo-likelihood sampler is used
#' instead (see \code{\link{MRFEngine}}).
#'
#' The paired \code{r} function is a stub: labels are always supplied as
#' initial values and updated by the Gibbs sweep, so exact Potts simulation is
#' never needed for inference. It exists because
#' \code{registerDistributions} requires one.
#'
#' Exported because NIMBLE's code generator resolves distribution names on the
#' search path; not intended to be called directly.
#'
#' @param x Vector of labels.
#' @param n Number of draws (required by NIMBLE).
#' @param beta Interaction strength.
#' @param e1,e2 Edge endpoints of the neighbourhood graph.
#' @param log Return the log density.
#' @return The (log) unnormalised Potts density, or a stub draw.
#' @references Besag, J. (1974). Spatial interaction and the statistical
#'   analysis of lattice systems. \emph{JRSS B}, 36(2), 192--236.
#'   \doi{10.1111/j.2517-6161.1974.tb00999.x}
#' @keywords internal
#' @name PottsKernels
NULL

#' @rdname PottsKernels
#' @export
dPottsNimix <- nimble::nimbleFunction(
  run = function(x = double(1), beta = double(0),
                 e1 = double(1), e2 = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    s <- 0
    nE <- length(e1)
    for (m in 1:nE) if (x[e1[m]] == x[e2[m]]) s <- s + 1
    lp <- beta * s
    if (log) return(lp) else return(exp(lp))
  })

#' @rdname PottsKernels
#' @export
rPottsNimix <- nimble::nimbleFunction(
  run = function(n = integer(0), beta = double(0),
                 e1 = double(1), e2 = double(1)) {
    returnType(double(1))
    out <- numeric(length = 1)
    return(out)
  })
