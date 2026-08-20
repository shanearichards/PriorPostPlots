# PriorPostPlots

An R package for comparing prior and posterior distributions of parameters from Stan models.

## Installation

The development version of `PriorPostPlots` can be installed from GitHub:

```r
install.packages("remotes")
remotes::install_github("shanearichards/PriorPostPlots")
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
