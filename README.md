# PriorPostPlots

A visual check that your data - not just your priors - are driving your Bayesian model's
conclusions. Give `PriorPosteriorPlot()` a fitted Stan model and the Stan code it was fitted from,
and it plots each parameter's prior against its posterior on the same axis.

**Note:** This package is still in development.

## Installation

Before installing `PriorPostPlots`, ensure you have the `rstan` package installed:

```r
install.packages("rstan")
```

The development version of `PriorPostPlots` can be installed from GitHub. Passing
`build_vignettes = TRUE` also builds the tutorial (see [Getting started](#getting-started) below)
so it's available locally via `vignette("tutorial", package = "PriorPostPlots")`:

```r
# install remotes if you don't have it
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# install from GitHub, including the tutorial vignette
remotes::install_github("shanearichards/PriorPostPlots", build_vignettes = TRUE)
```

## Troubleshooting

**RStudio's Help pane shows "Internal Server Error", or a topic is missing, after installing.**
This isn't specific to `PriorPostPlots` — it's a general RStudio/R quirk: RStudio caches a
help-topics index and runs a local help server, and reinstalling *any* package over a version
already loaded in that same live R session can leave that cache stale. It's most likely to show up
if you `remotes::install_github()` over a version of the package you'd already `library()`-loaded
earlier in the session. (As of v0.0.3, every `.Rd` file passes `tools::checkRd()` and `R CMD
check`'s documentation checks cleanly, so a genuinely malformed help page shipped with the package
is unlikely to be the cause.)

The fix is the same either way — restart R, then remove and reinstall forcing a rebuild:

```r
# restart R session (recommended, e.g. Session > Restart R in RStudio)
remove.packages("PriorPostPlots")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("shanearichards/PriorPostPlots", dependencies = TRUE, force = TRUE,
                         build = TRUE, build_vignettes = TRUE)
```

As a quick temporary workaround you can view any raw Rd directly from the repository and render it
locally without reinstalling, for example:

```r
tmp_rd <- tempfile(fileext = ".Rd")
download.file("https://raw.githubusercontent.com/shanearichards/PriorPostPlots/master/man/PriorPosteriorPlot.Rd", tmp_rd)
out <- tempfile(fileext = ".html")
tools::Rd2HTML(tmp_rd, out)
browseURL(out)
```

Developer note: Before pushing releases, regenerate documentation and rebuild the package so man/ and the compiled help DB are up-to-date:

```r
# from the package source directory
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("roxygen2", quietly = TRUE)) install.packages("roxygen2")

devtools::document()
devtools::check()
```

# Overview

`PriorPostPlots` has one function you need: **`PriorPosteriorPlot()`**. Give it a fitted Stan model
(`rstan` or `cmdstanr`) and the Stan code it was fitted from, and it reads the priors and hard
bounds straight out of the Stan code, reads the posterior draws straight out of the fit, and
produces one panel per parameter comparing the two - no intermediate data frame to build by hand.

Supported priors (Stan's own names/argument order): `normal`, `lognormal`, `exponential`, `gamma`,
`uniform`, `double_exponential` (Laplace), `beta`, `cauchy`, `student_t`, and `lkj_corr_cholesky`
(on a `cholesky_factor_corr[K]` parameter, plotted via its pairwise correlations). Hard bounds
(`<lower=, upper=>`) are detected automatically, including half-normal/half-Cauchy truncation.

# Getting started

```r
library(PriorPostPlots)

vignette("tutorial", package = "PriorPostPlots")  # step-by-step walkthrough

# full reference, including every supported prior's exact argument order
?PriorPosteriorPlot
```

# Author

Shane Richards

University of Tasmania
