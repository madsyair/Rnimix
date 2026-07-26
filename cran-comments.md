# cran-comments

## Test environments

* local Ubuntu 24.04, R 4.3.3
* (add win-builder / R-hub / macOS results before submitting)

## R CMD check results

0 errors | 0 warnings | 1 note (new submission)

### Note: new submission

This is a first submission, so `checking CRAN incoming feasibility` reports it
as such.

An earlier version of this package produced a note about
`integer(0, default = 0)` and about assignments to the global environment,
both arising from NIMBLE's requirement that user-defined distributions be
resolvable by name. Those kernels are now defined at build time and exported,
so no global assignment remains and the note no longer occurs.

## Notes seen only in restricted environments

If suggested packages are unavailable, an additional note lists them
(`rmarkdown`, `bayesplot`, `loo`, `cluster`, `fpc`, `sf`). All are used
conditionally, guarded by `requireNamespace()`, and all corresponding tests
and vignette chunks are skipped when the package is absent.

## Tests

`tests/testthat/` holds the automatic suite (distribution mathematics, S4
contracts, prior scaling, validation, diagnostics) and completes in well under
a minute.

MCMC-heavy tests — parameter recovery, engine integration, prediction and
workflow checks — live in `tests/manual/` and are **not** run by `R CMD check`,
because each fits real models and the suite takes tens of minutes. They are run
deliberately before releases via `Rscript tests/manual/run-manual-tests.R`.
