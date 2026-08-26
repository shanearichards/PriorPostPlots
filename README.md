# PriorPostPlots

An R package for comparing prior and posterior distributions of parameters from Stan models.

**Note:** This package is still in development.

## Installation

Before installing `PriorPostPlots`, ensure you have the `rstan` package installed:

```r
install.packages("rstan")
```

The development version of `PriorPostPlots` can be installed from GitHub:

```r
# install remotes if you don't have it
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# normal install from GitHub
remotes::install_github("shanearichards/PriorPostPlots")
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
remotes::install_github("shanearichards/PriorPostPlots", dependencies = TRUE, force = TRUE, build = TRUE)
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

`PriorPostPlots` provides functions for visualising prior and posterior distributions of selected parameters from Stan models.
The main functions are:

- `Create_df_priors()` — creates a data frame specifying prior distributions.
- `LongParameters()` — extracts selected parameters from Stan output.
- `PriorCurves()` — generates prior distribution curves.
- `PriorPosteriorPlotStan()` — produces plots comparing prior and posterior distributions.

# Example

```r
library(PriorPostPlots)

# See the function documentation for examples
?PriorPosteriorPlotStan
```

# Author

Shane Richards

University of Tasmania
