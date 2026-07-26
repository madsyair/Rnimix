#' @include dist-msnburr.R
#' @include dist-skewistudent-mv-og.R
#' @include dist-skewnormal-mv-og.R
#' @include dist-skewistudent-mv-o.R
#' @include dist-skewnormal-mv-o.R
#' @include dist-skewistudent-mv.R
#' @include dist-skewnormal-mv.R
#' @include dist-jfst.R
#' @include dist-fsst.R
#' @include dist-fossep.R
#' @include dist-fssn.R
#' @include dist-lep.R
#' @include dist-gmsnburr.R
#' @include dist-poisson-binomial.R
#' @include dist-mv-reg.R
#' @include dist-heavytail-reg.R
#' @include dist-glm-reg.R
#' @include dist-normal-gamma.R
#' @include dist-student-t.R
#' @include dist-student-t-mv.R
#' @include dist-normal-gamma-mv.R
NULL

## ---------------------------------------------------------------------------
## registerDistribution.R
##
## A tiny registry so advanced users can register their own DistributionSpec
## subclasses and refer to them by name in nimixClust(distribution = ...).
## Built-in distributions are registered on package load (see .onLoad).
## ---------------------------------------------------------------------------

.distRegistry <- new.env(parent = emptyenv())

#' Register a component distribution
#'
#' Adds a \code{\linkS4class{DistributionSpec}} to the registry under its
#' \code{name} slot so it can be selected by name. New built-in distributions
#' (Student-t, Poisson/Binomial) are planned for v0.4.0.
#'
#' @param spec A \code{\linkS4class{DistributionSpec}} instance.
#' @param overwrite Logical; overwrite an existing entry of the same name?
#' @return Invisibly, the registered name.
#' @examples
#' registerDistribution(NormalUvSpec(), overwrite = TRUE)
#' listDistributions()
#' @export
registerDistribution <- function(spec, overwrite = FALSE) {
  if (!methods::is(spec, "DistributionSpec"))
    stop("spec must inherit from DistributionSpec.", call. = FALSE)
  nm <- spec@name
  if (exists(nm, envir = .distRegistry, inherits = FALSE) && !overwrite)
    stop("Distribution '", nm, "' already registered; use overwrite = TRUE.",
         call. = FALSE)
  assign(nm, spec, envir = .distRegistry)
  invisible(nm)
}

#' Retrieve a registered distribution by name
#' @param name Character scalar.
#' @return A \code{\linkS4class{DistributionSpec}}.
#' @export
getDistribution <- function(name) {
  if (!exists(name, envir = .distRegistry, inherits = FALSE))
    stop("Unknown distribution '", name, "'. Registered: ",
         paste(listDistributions(), collapse = ", "), call. = FALSE)
  get(name, envir = .distRegistry, inherits = FALSE)
}

#' List registered distribution names
#' @return Character vector of registered names.
#' @export
listDistributions <- function() sort(ls(envir = .distRegistry))

#' @keywords internal
.nimixDefineMSNBurr <- function() {
  # Kernels are build-time objects now, so probing globalenv() would always
  # miss; the package-level flag is the reliable guard.
  if (isTRUE(.nimixState$msnburrDefined)) return(invisible())
  .nimixState$msnburrDefined <- TRUE
  # The scalar densities themselves are defined at build time in
  # kernels-neonormal.R and exported; only their registration happens here.
  # dMSNBurr_k / rMSNBurr_k: defined at build time in kernels-neonormal.R
  # (exported, so NIMBLE's code generator can resolve them without a global
  # assignment).

  # dMSNBurr2a_k: build time, kernels-neonormal.R
  # dGMSNBurr_k: build time, kernels-neonormal.R
  # rGMSNBurr_k: build time, kernels-neonormal.R
  # --- Batch B: SEP (symmetric exponential power) ---
  # dSEP_k: build time, kernels-neonormal.R
  # rSEP_k: build time, kernels-neonormal.R
  # --- Batch B: LEP (exponential power, alternative parameterisation) ---
  # dLEP_k: build time, kernels-neonormal.R
  # rLEP_k: build time, kernels-neonormal.R
  # --- Batch B: FSSN (Fernandez-Steel skew Normal) ---
  # dFSSN_k: build time, kernels-neonormal.R
  # rFSSN_k: build time, kernels-neonormal.R
  # --- Batch B: FOSSEP (Fernandez-Steel skew exponential power) ---
  # dFOSSEP_k: build time, kernels-neonormal.R
  # rFOSSEP_k: build time, kernels-neonormal.R
  # --- Batch B: FSST (Fernandez-Steel skew Student-t; t-kernel inlined) ---
  # dFSST_k: build time, kernels-neonormal.R
  # rFSST_k: build time, kernels-neonormal.R
  # --- Batch B: JFST (Jones-Faddy skew-t; branch-free rz) ---
  # dJFST_k: build time, kernels-neonormal.R
  # rJFST_k: build time, kernels-neonormal.R
  # --- Batch C: Ferreira-Steel skew multivariate Normal (A = chol(Sigma)) ---
  # dSkewMvN_k: build time, kernels-neonormal.R
  # rSkewMvN_k: build time, kernels-neonormal.R
  # --- Batch C: FS skew mv Normal with estimated orthogonal factor O (m = 2) ---
  # A = O U, U = chol(Sigma), O = I - 2 v v' the Householder reflection with
  # v = (sin theta, cos theta) (FS 2007 Appendix A). Note |O| = -1 always, so
  # O = I is NOT in the FS restricted set O_2; theta = 0 gives O = diag(1, -1),
  # which equals the O = I family with gamma_2 replaced by 1/gamma_2.
  # dSkewMvNO_k: build time, kernels-neonormal.R
  # rSkewMvNO_k: build time, kernels-neonormal.R
  # --- Batch C: FS skew mv Normal with estimated O, general m ---
  # O = O_{th^m} ... O_{th^2} applied to w = (U')^{-1}(x - mu); the blocks are
  # applied j = 2, ..., m, which reproduces the matrix product (verified against
  # the R reference and against the m = 2 kernel).
  # dSkewMvNOG_k: build time, kernels-neonormal.R
  # --- Batch C: FS skew mv independent-Student, estimated O, general m ---
  # dSkewMvITOG_k: build time, kernels-neonormal.R
  # --- Batch C: FS skew mv independent-Student with estimated O (m = 2) ---
  # dSkewMvITO_k: build time, kernels-neonormal.R
  # rSkewMvITO_k: build time, kernels-neonormal.R
  # --- Batch C: FS skew multivariate independent-Student (t-kernel inlined) ---
  # dSkewMvIT_k: build time, kernels-neonormal.R
  # rSkewMvIT_k: build time, kernels-neonormal.R
  # rMSNBurr2a_k: build time, kernels-neonormal.R
  invisible()
}

.onLoad <- function(libname, pkgname) {
  # The scalar densities are defined at build time in kernels-neonormal.R and
  # exported. A nimbleFunction that is merely a namespace object fails NIMBLE's
  # C++ code generation, which resolves names on the search path -- exporting
  # them is what makes the build-time definition work, with no global
  # assignment. Branch-free
  # softplus for stability. Iriawan (2000); Choir (2020).
  # Register built-ins. nimixClust() resolves "normal" to the univariate or
  # multivariate spec by data shape (see .selectClusterSpec); the "normal"
  # alias below is the univariate default for direct getDistribution() calls.
  assign("normal-uv", NormalUvSpec(), envir = .distRegistry)
  assign("normal-mv", NormalMvSpec(), envir = .distRegistry)
  assign("normal-reg", NormalRegSpec(), envir = .distRegistry)
  assign("student-t", StudentTUvSpec(), envir = .distRegistry)
  assign("normal-gamma", NormalGammaUvSpec(), envir = .distRegistry)
  assign("student-t-mv", StudentTMvSpec(), envir = .distRegistry)
  assign("normal-gamma-mv", NormalGammaMvSpec(), envir = .distRegistry)
  assign("poisson", PoissonSpec(), envir = .distRegistry)
  assign("binomial", BinomialSpec(), envir = .distRegistry)
  assign("poisson-reg", PoissonRegSpec(), envir = .distRegistry)
  assign("binomial-reg", BinomialRegSpec(), envir = .distRegistry)
  assign("msnburr-reg", MSNBurrRegSpec(), envir = .distRegistry)
  assign("sep-reg", SEPRegSpec(), envir = .distRegistry)
  assign("msnburr2a-reg", MSNBurr2aRegSpec(), envir = .distRegistry)
  assign("fssn-reg", FSSNRegSpec(), envir = .distRegistry)
  assign("gmsnburr-reg", GMSNBurrRegSpec(), envir = .distRegistry)
  assign("lep-reg", LEPRegSpec(), envir = .distRegistry)
  assign("fsst-reg", FSSTRegSpec(), envir = .distRegistry)
  assign("fossep-reg", FOSSEPRegSpec(), envir = .distRegistry)
  assign("jfst-reg", JFSTRegSpec(), envir = .distRegistry)
  assign("student-t-reg", StudentTRegSpec(), envir = .distRegistry)
  assign("normal-gamma-reg", NormalGammaRegSpec(), envir = .distRegistry)
  assign("normal-mv-reg", NormalMvRegSpec(), envir = .distRegistry)
  assign("student-t-mv-reg", StudentTMvRegSpec(), envir = .distRegistry)
  assign("normal-gamma-mv-reg", NormalGammaMvRegSpec(), envir = .distRegistry)
  assign("msnburr", MSNBurrUvSpec(), envir = .distRegistry)
  assign("msnburr2a", MSNBurr2aUvSpec(), envir = .distRegistry)
  assign("gmsnburr", GMSNBurrUvSpec(), envir = .distRegistry)
  assign("sep", SEPUvSpec(), envir = .distRegistry)
  assign("lep", LEPUvSpec(), envir = .distRegistry)
  assign("fssn", FSSNUvSpec(), envir = .distRegistry)
  assign("fossep", FOSSEPUvSpec(), envir = .distRegistry)
  assign("fsst", FSSTUvSpec(), envir = .distRegistry)
  assign("jfst", JFSTUvSpec(), envir = .distRegistry)
  assign("skewnormal-mv", SkewNormalMvSpec(), envir = .distRegistry)
  assign("skewistudent-mv", SkewIStudentMvSpec(), envir = .distRegistry)
  assign("skewnormal-mv-o", SkewNormalMvOSpec(), envir = .distRegistry)
  assign("skewistudent-mv-o", SkewIStudentMvOSpec(), envir = .distRegistry)
  assign("skewnormal-mv-og", SkewNormalMvOGenSpec(), envir = .distRegistry)
  assign("skewistudent-mv-og", SkewIStudentMvOGenSpec(), envir = .distRegistry)
  assign("normal", NormalUvSpec(), envir = .distRegistry)
  # Register the user-defined multivariate-t density with NIMBLE so the
  # StudentTMvSpec kernel resolves at model-build time.
  suppressMessages(suppressWarnings(try(
    nimble::registerDistributions(list(
      dmvt_nimix = list(
        BUGSdist = "dmvt_nimix(mu, cov, df)",
        types = c("value = double(1)", "mu = double(1)",
                  "cov = double(2)", "df = double(0)")))),
    silent = TRUE)))
  # The unnormalised Potts prior for the MRF engine is built and registered
  # lazily in globalenv by .nimixDefinePotts()/.nimixEnsureMSNBurr(); doing it
  # here (namespace frame) makes NIMBLE fail to find rPottsNimix at code-gen.
  invisible(NULL)
}

# Lazily ensure the MSNBurr densities exist in the global environment AND are
# registered with NIMBLE. Registration is deferred to first use (not .onLoad)
# because a registration performed while the objects are being (re)built during
# package load can bind a distribution name to a namespace-frame object that
# fails C++ code generation; binding the name once, here, to the global-frame
# objects avoids that. Idempotent and cheap.
.nimixEnsureMSNBurr <- function() {
  .nimixDefineMSNBurr()
  .nimixDefinePotts()
  if (isTRUE(.nimixState$msnburrRegistered)) return(invisible())
  # Registration used to be evaluated in globalenv() so that it could see the
  # kernels; they are exported now, so the search path resolves them and a
  # plain call suffices.
  suppressMessages(suppressWarnings(try(nimble::registerDistributions(list(
    dMSNBurr_k = list(
      BUGSdist = "dMSNBurr_k(mu, sigma, alpha)",
      types = c("value = double(0)", "mu = double(0)",
                "sigma = double(0)", "alpha = double(0)"),
      discrete = FALSE),
    dMSNBurr2a_k = list(
      BUGSdist = "dMSNBurr2a_k(mu, sigma, alpha)",
      types = c("value = double(0)", "mu = double(0)",
                "sigma = double(0)", "alpha = double(0)"),
      discrete = FALSE),
    dGMSNBurr_k = list(
      BUGSdist = "dGMSNBurr_k(mu, sigma, alpha, theta)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)", "theta = double(0)"),
      discrete = FALSE),
    dPottsNimix = list(
      BUGSdist = "dPottsNimix(beta, e1, e2)",
      types = c("value = double(1)", "beta = double(0)",
                "e1 = double(1)", "e2 = double(1)"),
      discrete = TRUE, mixedSizes = TRUE),
    dSEP_k = list(
      BUGSdist = "dSEP_k(mu, sigma, nu)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "nu = double(0)"), discrete = FALSE),
    dLEP_k = list(
      BUGSdist = "dLEP_k(mu, sigma, nu)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "nu = double(0)"), discrete = FALSE),
    dFSSN_k = list(
      BUGSdist = "dFSSN_k(mu, sigma, alpha)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)"), discrete = FALSE),
    dFOSSEP_k = list(
      BUGSdist = "dFOSSEP_k(mu, sigma, alpha, theta)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)", "theta = double(0)"), discrete = FALSE),
    dFSST_k = list(
      BUGSdist = "dFSST_k(mu, sigma, alpha, nu)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)", "nu = double(0)"), discrete = FALSE),
    dJFST_k = list(
      BUGSdist = "dJFST_k(mu, sigma, alpha, theta)",
      types = c("value = double(0)", "mu = double(0)", "sigma = double(0)",
                "alpha = double(0)", "theta = double(0)"),
      discrete = FALSE),
    dSkewMvN_k = list(
      BUGSdist = "dSkewMvN_k(mu, Sigma, gam)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)"),
      discrete = FALSE),
    dSkewMvIT_k = list(
      BUGSdist = "dSkewMvIT_k(mu, Sigma, gam, nu)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "nu = double(1)"),
      discrete = FALSE),
    dSkewMvNO_k = list(
      BUGSdist = "dSkewMvNO_k(mu, Sigma, gam, theta)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "theta = double(0)"),
      discrete = FALSE),
    dSkewMvITO_k = list(
      BUGSdist = "dSkewMvITO_k(mu, Sigma, gam, nu, theta)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "nu = double(1)", "theta = double(0)"),
      discrete = FALSE),
    dSkewMvNOG_k = list(
      BUGSdist = "dSkewMvNOG_k(mu, Sigma, gam, theta)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "theta = double(1)"),
      discrete = FALSE),
    dSkewMvITOG_k = list(
      BUGSdist = "dSkewMvITOG_k(mu, Sigma, gam, nu, theta)",
      types = c("value = double(1)", "mu = double(1)", "Sigma = double(2)",
                "gam = double(1)", "nu = double(1)", "theta = double(1)"),
      discrete = FALSE))), silent = TRUE)))
  .nimixState$msnburrRegistered <- TRUE
  invisible()
}

# Build the Potts prior's d/r functions in the GLOBAL environment. Registering
# them from a namespace frame makes NIMBLE fail to find rPottsNimix during code
# generation for the latent label node (it is invoked to simulate z), the same
# class of failure that affects the scalar neo-normal densities. Building here,
# in globalenv, resolves it.
# The Potts d/r pair is defined at build time in kernels-neonormal.R (exported,
# so NIMBLE's code generator can resolve it). Kept as a no-op so existing call
# sites stay valid.
.nimixDefinePotts <- function() invisible()

.nimixState <- new.env(parent = emptyenv())
