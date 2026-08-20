#' @importFrom rlang .data
NULL

# Avoid R CMD check notes for variables used inside tidyverse/ggplot2
# non-standard evaluation.
utils::globalVariables(c("par", "x", "y", "f", "val", "density"))

#' Probability density for the Laplace distribution
#'
#' @description
#' Probability density for the Laplace distribution.
#'
#' @param x vector of quantities.
#' @param mu vector of means.
#' @param b vector of inverse decay rates.
#'
#' @return Vector of densities.
#' @details
#' Laplace probability density function has two arguments, the mean (mu) and the decay rate (b), and is described by
#'  exp(-abs(x - mu) / b) / (2 * b).
#'
#' @export
dlaplace <- function(x, mu = 0, b = 1) {
  return(1 / (2 * b) * exp(-abs(x - mu) / b))
}

#' Cumulative probability for the Laplace distribution
#'
#' @description
#' Cumulative probability for the Laplace distribution.
#'
#' @param x vector of quantities.
#' @param mu vector of means.
#' @param b vector of inverse decay rates.
#'
#' @return Vector of probabilities.
#' @details
#' Laplace probability density function has two arguments, the mean (mu) and the decay rate (b), and is described by
#'  exp(-abs(x - mu) / b) / (2 * b).
#'
#' @export
plaplace <- function(x, mu = 0, b = 1) {
  ifelse(x < mu,
    0.5 * exp((x - mu) / b),
    1 - 0.5 * exp(-(x - mu) / b))
}

#' Generates random draws from the Laplace distribution
#'
#' @description
#' Generates random draws from the Laplace distribution.
#'
#' @param n number of observations.
#' @param mu mean.
#' @param b decay rate.
#'
#' @return Vector of observations.
#' @details
#' Laplace probability density function has two arguments, the mean (mu) and the decay rate (b), and is described by
#'  exp(-abs(x - mu) / b) / (2 * b).
#'
#' @export
rlaplace <- function(n, mu = 0, b = 1) {
  u <- stats::runif(n, min = -0.5, max = 0.5)
  return(mu - b * sign(u) * log(1 - 2 * abs(u)))
}

# convert list of posterior estimates into long format data frame
LongParameters <- function(l_params, pars_use) {

  # Flatten any matrix/array-valued entries into individual scalar
  # columns, matching Create_df_priors()'s naming convention:
  #   - vector[N] parameters come back from rstan::extract() as an
  #     iterations x N matrix -> "name[1]", "name[2]", ...
  #   - matrix[R, C] parameters come back as an iterations x R x C
  #     array -> "name[1,1]", "name[2,1]", ..., "name[R,1]", "name[1,2]",
  #     ... (column-major, matching Stan's own parameter-name ordering).
  flat_list <- list()
  for (nm in names(l_params)) {
    val <- l_params[[nm]]
    if (is.matrix(val)) {
      for (j in seq_len(ncol(val))) {
        flat_list[[paste0(nm, "[", j, "]")]] <- val[, j]
      }
    } else if (is.array(val) && length(dim(val)) == 3) {
      d  <- dim(val)
      nr <- d[2]
      nc <- d[3]
      for (j in seq_len(nc)) {
        for (i in seq_len(nr)) {
          flat_list[[paste0(nm, "[", i, ",", j, "]")]] <- val[, i, j]
        }
      }
    } else {
      flat_list[[nm]] <- val
    }
  }

  df_params <- flat_list |>
    tibble::as_tibble() |>
    dplyr::select(tidyselect::all_of(pars_use)) |>
    tidyr::pivot_longer(names_to = "par", values_to = "val", 1:length(pars_use))

  df_params$par <- factor(df_params$par, levels = pars_use)

  df_params <- df_params |>
    dplyr::arrange(.data$par, .data$val)

  return(df_params)
}

# Create theoretical curves data frame
PriorCurves <- function(df_param) {

  df_curves <- purrr::pmap_dfr(
    list(par = df_param$par,
    	dist     = df_param$dist,
      arg1   = df_param$arg1,
    	arg2   = df_param$arg2,
    	v_min  = df_param$v_min,
    	v_max  = df_param$v_max,
  	  v_lwr  = pmax(df_param$v_min, df_param$v_lwr, na.rm=TRUE),
  	  v_upr  = pmin(df_param$v_max, df_param$v_upr, na.rm=TRUE)
    ),

    function(par, dist, arg1, arg2, v_min, v_max, v_lwr, v_upr) {
      x_vals <- seq(v_min, v_max, length.out = 100)

      y_vals <- switch(dist,
        normal      = stats::dnorm(x_vals, mean = arg1, sd = arg2),
        exponential = stats::dexp(x_vals, rate = arg1),
        gamma       = stats::dgamma(x_vals, shape = arg1, rate = arg2),
    	  log_normal  = stats::dlnorm(x_vals, meanlog = arg1, sdlog = arg2),
      	laplace     = dlaplace(x_vals, mu = arg1, b = arg2),
      	uniform     = stats::dunif(x_vals, min = arg1, max = arg2),
      	beta        = stats::dbeta(x_vals, shape1 = arg1, shape2 = arg2),
        rep(NA, length(x_vals))
      )

      d <- switch(dist,
        normal = stats::pnorm(v_upr, mean = arg1, sd = arg2) -
          stats::pnorm(v_lwr, mean = arg1, sd = arg2),
        exponential = stats::pexp(v_upr, rate = arg1) -
          stats::pexp(v_lwr, rate = arg1),
        gamma = stats::pgamma(v_upr, shape = arg1, rate = arg2) -
          stats::pgamma(v_lwr, shape = arg1, rate = arg2),
    	  log_normal = stats::plnorm(v_upr, meanlog = arg1, sdlog = arg2) -
          stats::plnorm(v_lwr, meanlog = arg1, sdlog = arg2),
    	  laplace = plaplace(v_upr, mu = arg1, b = arg2) -
          plaplace(v_lwr, mu = arg1, b = arg2),
      	uniform = stats::punif(v_upr, min = arg1, max = arg2) -
          stats::punif(v_lwr, min = arg1, max = arg2),
        beta = stats::pbeta(v_upr, shape1 = arg1, shape2 = arg2) -
          stats::pbeta(v_lwr, shape1 = arg1, shape2 = arg2),
        rep(NA, length(x_vals))
      )

      tibble::tibble(par = par, x = x_vals, y = y_vals, d = d) |>
        dplyr::mutate(f = .data$y / .data$d) |>
        dplyr::select(.data$par, .data$x, y = .data$f)
    }
  )

  return(df_curves)
}

#' Generates prior-posterior plots
#'
#' @description
#' Generates a ggplot comparing prior and posterior distributions by combining sampling output from a Stan model with a data frame specifying the prior distributions for selected model parameters.
#'
#' @param stan_fit An object returned by rstan::sampling.
#' @param pars A vector of parameter names associated with the stan model for plotting. Entries may be exact \code{df_priors$par} values (e.g. a scalar parameter, or an already bracket-indexed element like \code{"etaR[1]"}), or the *base* name of a vector/matrix parameter (e.g. \code{"etaR"}, \code{"beta_genus"}), which expands to all of that parameter's elements as they appear in \code{df_priors} (e.g. \code{"etaR[1]"}, \code{"etaR[2]"}, \code{"etaR[3]"}) -- so elements don't need to be listed individually.
#' @param df_priors A data frame containing information that describes the prior for each parameter. The data frame must have the columns par, dist, arg1, arg2, v_min, v_max, v_lwr, v_upr. par = parameter name, dist = prior distribution (normal, log_normal, exponential, gamma, uniform, laplace, beta), arg1 = first distribution parameter, arg2 = second distribution parameter (NA if not needed), v_min and v_max are the bounds of the plotted prior, v_lwr and v_upr are the stan-imposed parameter bounds (NA if none are set).
#' @param ncol Number of columns provided to facet_wrap. Square arrangement is produced when no value is provided.
#' @param nbins Number of bins used to display histograms (default is 25).
#'
#' @return A ggplot object.
#' @details
#' Currently, the distributions allowed are uniform, normal, log_normal, exponential,
#'  gamma, laplace and beta. See corresponding d-functions for arguments required.
#'  Laplace has two arguments, the mean (mu) and the decay rate (b), and is described by
#'  exp(-abs(x - mu) / b) / (2 * b).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pars <- c("A", "B", "C", "E") # parameters to plot
#' df_priors <- data.frame(
#'   par   = c("A",  "B",  "C",  "D", "E", "F", "G"), # parameter names
#'   v_min = c(  0, 0.01, 0.01, 0.01,  -2,   0, 0.1), # display range max
#'   v_max = c(  5,   20,    9,   50,   1,  20, 0.8), # display range min
#'   v_lwr = c(  0, 0.01, 0.01, 0.01,  NA,   0, 0.1), # stan lower bound, NA?
#'   v_upr = c(  5,   NA,   25,   60,  NA,  20, 0.8), # stan upper bound, NA?
#'   dist  = c("normal", "exponential", "gamma",
#'             "log_normal", "laplace", "uniform",
#'             "beta"),                               # distribution
#'   arg1  = c(  0,  0.5,    2,  1.5,   0,   0,   2), # distribution argument 1
#'   arg2  = c(  3,   NA,  0.5,  1.0, 0.5,  10,   3)  # distribution argument 2
#' )
#' PriorPosteriorPlotStan(stan_fit, pars, df_priors) # produce the prior-post plot
#' }
PriorPosteriorPlotStan <- function(stan_fit, pars, df_priors, ncol = NA, nbins = 25) {
  # check validy of stan_fit
	if (!inherits(stan_fit, "stanfit")) {
    stop("Error: 'stan_fit' argument must be an object of class 'stanfit'.")
  }

	# check validity of ncol
  if (!is.na(ncol) && (!is.numeric(ncol) || ncol %% 1 != 0 || ncol <= 0)) {
    stop("Error: 'ncol' must be NA or a positive integer.")
  }

  if (!(is.data.frame(df_priors) || tibble::is_tibble(df_priors))) {
    stop("Error: 'df_priors' argument must be either a data frame or a tibble.")
  }

  # check required column names are present in df_priors
	required_names <- c("par", "v_min", "v_max", "v_lwr", "v_upr", "dist",
    "arg1", "arg2")

  if (!all(required_names %in% colnames(df_priors))) {
    stop("Error: df_priors argument does not contain all required column names.")
  }

  pars_use <- unique(pars) # remove any duplicate parameter names

  # Allow a 'pars' entry to be either an exact match to a df_priors$par
  # value (e.g. a scalar parameter "sigma", or an already bracket-indexed
  # element like "etaR[1]"), OR a vector/matrix *base* name (e.g. "etaR",
  # "beta_genus") that stands in for all of its expanded elements
  # ("etaR[1]", "etaR[2]", ..., "beta_genus[1,1]", "beta_genus[2,1]", ...)
  # as they appear in df_priors$par. This means callers don't have to
  # enumerate every vector/matrix element by hand.
  all_par_names <- as.character(df_priors$par)
  pars_use <- unlist(lapply(pars_use, function(p) {
    if (p %in% all_par_names) return(p)
    matches <- all_par_names[grepl(paste0("^", p, "\\["), all_par_names, perl = TRUE)]
    if (length(matches) == 0) {
      warning(
        "Parameter '", p, "' was not found in 'df_priors' (checked for ",
        "an exact match and for '", p, "[...]' elements); it will be dropped."
      )
      return(character(0))
    }
    matches
  }))
  pars_use <- unique(pars_use)

  # set v_lwr and v_upr to v_min and v_max if NA
  df_priors_use <- df_priors |>
    dplyr::filter(.data$par %in% pars_use)

  if (!(nrow(df_priors_use) > 0)) { # at least 1 parameter is present
    stop("Error: no valid model parameters provided in pars.")
  }

  df_priors_use$par <- factor(df_priors_use$par, levels = pars_use)
	df_priors_use$v_lwr <- pmax(df_priors_use$v_min, df_priors_use$v_lwr, na.rm=TRUE)
	df_priors_use$v_upr <- pmin(df_priors_use$v_max, df_priors_use$v_upr, na.rm=TRUE)

  # extract posterior samples and place in long format
  # rstan::extract() only accepts base parameter names (e.g. "beta"),
  # not bracket-indexed names (e.g. "beta[1]") -- strip any "[...]"
  # suffix before requesting from rstan; LongParameters() re-expands
  # any resulting matrix-valued entries back into "name[i]" columns.
  base_pars <- unique(sub("\\[.*\\]$", "", pars_use))
  l_posteriors <- rstan::extract(object = stan_fit, pars = base_pars)
  df_posteriors_plot <- LongParameters(l_posteriors, pars_use)

  if (!(nrow(df_posteriors_plot) > 0)) { # at least 1 parameter is present
    stop("Error: missing prior information")
  }

  # generate data for prior distribution curves in long format
  df_priors_plot <- PriorCurves(df_priors_use)

  # set the number of facet columns
  ncol <- ifelse(is.na(ncol), ceiling(0.5*length(pars_use)), ncol)

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = df_priors_plot,
      ggplot2::aes(x = .data$x, ymin = 0, ymax = .data$y), fill = "salmon") +
    ggplot2::geom_histogram(data = df_posteriors_plot,
      ggplot2::aes(x = .data$val, y = ggplot2::after_stat(.data$density)), bins = nbins,
      fill = "grey85", color = "black") +
    ggplot2::facet_wrap(ggplot2::vars(par), scales = "free", ncol = ncol) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white", color = "black")) +
    ggplot2::labs(x = "Parameter value", y = "Density")
}

# =============================================================================
# Create_df_priors.R
#
# Parses Stan model code to automatically build the df_priors data frame
# required by PriorPosteriorPlotStan(). Extracts scalar `real` parameters
# and their hard bounds from the `parameters` block, and their sampling
# distribution (if any) from the `model` block.
# =============================================================================

# ---- quantile function for the Laplace distribution ------------------------

#' Quantile function for the Laplace distribution
#'
#' @param p vector of probabilities.
#' @param mu vector of means.
#' @param b vector of decay rates.
#' @return Vector of quantiles.
#' @keywords internal
qlaplace <- function(p, mu = 0, b = 1) {
  ifelse(p < 0.5,
    mu + b * log(2 * p),
    mu - b * log(2 * (1 - p)))
}


# ---- low-level helpers (not exported) ---------------------------------------

#' Strip // and /* */ comments from Stan code
#' @keywords internal
strip_comments <- function(code) {
  code <- gsub("//.*", "", code, perl = TRUE)
  code <- gsub("(?s)/\\*.*?\\*/", "", code, perl = TRUE)
  code
}

#' Extract a named block (e.g. "parameters", "model") from Stan code,
#' using brace matching so nested braces (loops, if-statements) don't
#' truncate the block early.
#' @keywords internal
extract_block <- function(code, block_name) {
  pattern <- paste0("\\b", block_name, "\\b\\s*\\{")
  m <- regexpr(pattern, code, perl = TRUE)
  if (m[1] == -1) return(NA_character_)

  brace_start <- m[1] + attr(m, "match.length") - 1  # index of the opening "{"
  chars <- strsplit(code, "")[[1]]

  depth <- 0
  end <- NA_integer_
  for (i in brace_start:length(chars)) {
    if (chars[i] == "{") depth <- depth + 1
    if (chars[i] == "}") depth <- depth - 1
    if (depth == 0) {
      end <- i
      break
    }
  }

  if (is.na(end)) {
    stop("Unbalanced braces while extracting the '", block_name, "' block.")
  }

  substr(code, brace_start + 1, end - 1)
}

#' Try to resolve a Stan expression (e.g. "0.5" or "sigma_prior") to a
#' number, optionally looking it up in a supplied data list.
#' @keywords internal
resolve_numeric <- function(x, data_list = NULL, par_name = NULL) {
  x <- trimws(x)
  val <- suppressWarnings(as.numeric(x))
  if (!is.na(val)) return(val)

  if (!is.null(data_list) && x %in% names(data_list)) {
    val <- suppressWarnings(as.numeric(data_list[[x]]))
    if (!is.na(val)) return(val)
  }

  warning(
    "Could not resolve value '", x, "'",
    if (!is.null(par_name)) paste0(" for parameter '", par_name, "'") else "",
    " to a number (it is not a literal and was not found in 'data_list'). ",
    "Setting to NA."
  )
  NA_real_
}

#' Parse the parameters block: identify scalar `real`, `vector[N]` and
#' `matrix[R, C]` parameters, and any hard lower/upper bounds declared on
#' them.
#' @keywords internal
parse_parameters_block <- function(block, data_list = NULL) {
  statements <- strsplit(block, ";")[[1]]
  statements <- trimws(gsub("\\s+", " ", statements))
  statements <- statements[statements != ""]

  unsupported_types <- c("array", "int", "simplex",
    "ordered", "positive_ordered", "row_vector", "cholesky_factor_corr",
    "cholesky_factor_cov", "corr_matrix", "cov_matrix", "unit_vector")

  results <- list()

  add_result <- function(pname, base, v_lwr, v_upr) {
    results[[length(results) + 1]] <<- data.frame(
      par = pname, base = base, v_lwr = v_lwr, v_upr = v_upr,
      stringsAsFactors = FALSE
    )
  }

  parse_bounds <- function(bounds_str, data_list) {
    v_lwr <- NA_real_
    v_upr <- NA_real_
    if (!is.na(bounds_str) && nzchar(bounds_str)) {
      lwr_m <- regmatches(bounds_str, regexpr("lower\\s*=\\s*[^,>]+", bounds_str))
      upr_m <- regmatches(bounds_str, regexpr("upper\\s*=\\s*[^,>]+", bounds_str))
      if (length(lwr_m) > 0) {
        v_lwr <- resolve_numeric(trimws(sub("lower\\s*=\\s*", "", lwr_m)), data_list)
      }
      if (length(upr_m) > 0) {
        v_upr <- resolve_numeric(trimws(sub("upper\\s*=\\s*", "", upr_m)), data_list)
      }
    }
    list(v_lwr = v_lwr, v_upr = v_upr)
  }

  for (stmt in statements) {

    # vector declaration, e.g. "vector<lower=0>[4] beta" or "vector[K] beta, gamma"
    m_vec <- regmatches(stmt, regexec(
      "^vector(<([^>]*)>)?\\s*\\[\\s*([^\\]]+)\\s*\\]\\s+(.+)$", stmt,
      perl = TRUE))[[1]]

    if (length(m_vec) > 0) {
      bounds <- parse_bounds(m_vec[3], data_list)
      size_expr <- trimws(m_vec[4])
      names_str <- sub("=.*$", "", m_vec[5])  # drop any "= expr" definition

      n_elem <- resolve_numeric(size_expr, data_list)
      if (is.na(n_elem)) {
        warning(
          "Could not resolve size '", size_expr, "' for vector parameter(s) '",
          trimws(names_str), "'; pass its value via 'data_list' if it's a ",
          "data-block constant. Skipping."
        )
        next
      }
      n_elem <- as.integer(round(n_elem))

      for (pname in trimws(strsplit(names_str, ",")[[1]])) {
        if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", pname)) next
        for (i in seq_len(n_elem)) {
          add_result(paste0(pname, "[", i, "]"), pname, bounds$v_lwr, bounds$v_upr)
        }
      }
      next
    }

    # matrix declaration, e.g. "matrix<lower=0>[R, C] beta_genus" or
    # "matrix[n_genus_minus1, 4] beta_genus, other_mat"
    m_mat <- regmatches(stmt, regexec(
      "^matrix(<([^>]*)>)?\\s*\\[\\s*([^,\\]]+)\\s*,\\s*([^\\]]+)\\s*\\]\\s+(.+)$",
      stmt, perl = TRUE))[[1]]

    if (length(m_mat) > 0) {
      bounds <- parse_bounds(m_mat[3], data_list)
      nrow_expr <- trimws(m_mat[4])
      ncol_expr <- trimws(m_mat[5])
      names_str <- sub("=.*$", "", m_mat[6])  # drop any "= expr" definition

      n_row <- resolve_numeric(nrow_expr, data_list)
      n_col <- resolve_numeric(ncol_expr, data_list)
      if (is.na(n_row) || is.na(n_col)) {
        warning(
          "Could not resolve dimensions '[", nrow_expr, ", ", ncol_expr,
          "]' for matrix parameter(s) '", trimws(names_str), "'; pass ",
          "the missing value(s) via 'data_list' if they are data-block ",
          "constants. Skipping."
        )
        next
      }
      n_row <- as.integer(round(n_row))
      n_col <- as.integer(round(n_col))

      # Elements are named "name[i,j]" in column-major order (j outer,
      # i inner), matching Stan's own parameter-naming convention, and
      # matching how LongParameters() flattens the corresponding
      # iterations x R x C array returned by rstan::extract().
      for (pname in trimws(strsplit(names_str, ",")[[1]])) {
        if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", pname)) next
        for (j in seq_len(n_col)) {
          for (i in seq_len(n_row)) {
            add_result(paste0(pname, "[", i, ",", j, "]"), pname,
              bounds$v_lwr, bounds$v_upr)
          }
        }
      }
      next
    }

    # scalar real declaration, e.g. "real<lower=0,upper=1> a, b"
    m_real <- regmatches(stmt, regexec("^real(<([^>]*)>)?\\s+(.+)$", stmt))[[1]]

    if (length(m_real) > 0) {
      bounds <- parse_bounds(m_real[3], data_list)
      names_str <- sub("=.*$", "", m_real[4])  # drop any "= expr" definition

      for (pname in trimws(strsplit(names_str, ",")[[1]])) {
        if (grepl("\\[", pname)) {
          warning(
            "Parameter '", pname, "' appears to be array/indexed-valued; ",
            "skipping (only scalar 'real' and 'vector' parameters are supported)."
          )
          next
        }
        if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", pname)) next
        add_result(pname, pname, bounds$v_lwr, bounds$v_upr)
      }
      next
    }

    if (any(startsWith(stmt, unsupported_types))) {
      warning(
        "Parameter declaration '", stmt, "' has an unsupported type; ",
        "only scalar 'real', 'vector' and 'matrix' parameters are ",
        "currently handled by Create_df_priors(). Skipping."
      )
    }
  }

  if (length(results) == 0) {
    stop("No supported ('real', 'vector' or 'matrix') parameters were found in the block.")
  }

  df <- do.call(rbind, results)
  rownames(df) <- NULL
  df
}

#' Parse the model block for `par ~ dist(args);` sampling statements
#' associated with the supplied parameter names. Also recognises the
#' `to_vector(par) ~ dist(args);` idiom commonly used to put a vectorized
#' prior on a whole `matrix` parameter -- `to_vector(...)` is stripped and
#' `par` is treated exactly as if it had been written bare on the LHS
#' (i.e. the prior is applied to every expanded element of `par`).
#'
#' Note: this does not (yet) detect priors specified via target +=
#' *_lpdf(...)/*_lpmf(...) syntax, nor does it adjust for truncation
#' (`T[lower, upper]`) -- these fall back to "no prior" (dist = NA).
#' @keywords internal
parse_priors_block <- function(block, par_names) {
  pattern <- paste0(
    "(?:to_vector\\(\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\)|",
    "([A-Za-z_][A-Za-z0-9_]*))",
    "\\s*~\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\(([^)]*)\\)"
  )
  m <- gregexpr(pattern, block, perl = TRUE)
  matches <- regmatches(block, m)[[1]]

  priors <- list()

  for (mt in matches) {
    parts <- regmatches(mt, regexec(pattern, mt, perl = TRUE))[[1]]
    pname     <- if (nzchar(parts[2])) parts[2] else parts[3]
    dist_stan <- parts[4]
    args_str  <- parts[5]

    if (!(pname %in% par_names)) next

    args <- trimws(strsplit(args_str, ",")[[1]])

    if (pname %in% names(priors)) {
      warning(
        "Multiple sampling statements found for parameter '", pname,
        "'; using the first one encountered."
      )
      next
    }

    priors[[pname]] <- list(dist_stan = dist_stan, args = args)
  }

  priors
}

#' Map of supported Stan distribution names to the internal `dist` labels
#' used by PriorPosteriorPlotStan(), and how many arguments each needs.
#' @keywords internal
.dist_map <- list(
  normal             = list(dist = "normal",      nargs = 2),
  lognormal          = list(dist = "log_normal",  nargs = 2),
  exponential        = list(dist = "exponential", nargs = 1),
  gamma               = list(dist = "gamma",       nargs = 2),
  uniform            = list(dist = "uniform",     nargs = 2),
  double_exponential = list(dist = "laplace",     nargs = 2),
  beta               = list(dist = "beta",        nargs = 2)
)

#' Evaluate the quantile function for one of the seven supported
#' distributions, using PriorPosteriorPlotStan()'s internal naming
#' convention (dist, arg1, arg2).
#' @keywords internal
quantile_fn <- function(dist, p, arg1, arg2) {
  switch(dist,
    normal      = stats::qnorm(p, mean = arg1, sd = arg2),
    log_normal  = stats::qlnorm(p, meanlog = arg1, sdlog = arg2),
    exponential = stats::qexp(p, rate = arg1),
    gamma       = stats::qgamma(p, shape = arg1, rate = arg2),
    uniform     = stats::qunif(p, min = arg1, max = arg2),
    laplace     = qlaplace(p, mu = arg1, b = arg2),
    beta        = stats::qbeta(p, shape1 = arg1, shape2 = arg2),
    NA_real_
  )
}


# ---- main function -----------------------------------------------------------

#' Build a df_priors data frame directly from Stan model code
#'
#' @description
#' Parses Stan code to identify scalar `real` model parameters, any hard
#' bounds declared for them in the `parameters` block, and their sampling
#' distribution (if any) declared in the `model` block via `par ~ dist(...)`
#' statements. The result is a data frame in the format required by the
#' `df_priors` argument of \code{PriorPosteriorPlotStan()}.
#'
#' @param stan_code A file path to a .stan file, a character string/vector
#'   containing Stan code, or a compiled \code{stanmodel} object (e.g. the
#'   object created by \code{rstan::stan_model()}, or by an R Markdown Stan
#'   code chunk with an \code{output.var} argument -- in the latter case,
#'   pass the resulting output object directly).
#' @param data_list Optional named list (e.g. the same list passed to
#'   \code{rstan::sampling(data = ...)}). Used to resolve prior arguments
#'   or parameter bounds that reference named constants (e.g.
#'   \code{normal(mu_prior, sigma_prior)}) rather than numeric literals.
#'
#' @return A tibble with columns \code{par, v_min, v_max, v_lwr, v_upr,
#'   dist, arg1, arg2}, suitable for use as the \code{df_priors} argument
#'   to \code{PriorPosteriorPlotStan()}.
#'
#' @details
#' Currently supported prior distributions (matching Stan's names on the
#' left): \code{normal}, \code{lognormal} (-> \code{log_normal}),
#' \code{exponential}, \code{gamma}, \code{uniform},
#' \code{double_exponential} (-> \code{laplace}), \code{beta}.
#'
#' Scalar \code{real}, \code{vector[N]} and \code{matrix[R, C]} parameters
#' are currently supported (not \code{array}, \code{int}, \code{simplex},
#' etc.) -- unsupported types are skipped with a warning. \code{vector[N]}
#' and \code{matrix[R, C]} parameters (declared in either the
#' \code{parameters} or \code{transformed parameters} block) are expanded
#' into one row per element: \code{vector[N]} elements are named
#' \code{"name[1]"}, \code{"name[2]"}, ...; \code{matrix[R, C]} elements
#' are named \code{"name[1,1]"}, \code{"name[2,1]"}, ..., \code{"name[R,1]"},
#' \code{"name[1,2]"}, ... (column-major order, matching Stan's own
#' parameter-naming convention and what \code{rstan::extract()} returns).
#' If \code{N}, \code{R} or \code{C} reference a data-block constant
#' rather than a literal (e.g. \code{vector[K] beta} or
#' \code{matrix[R, C] beta_genus}), pass its value via \code{data_list}.
#'
#' A vectorized sampling statement (e.g. \code{beta ~ normal(0, 5);}) is
#' applied to every element of that vector. Per-element statements written
#' inside a loop (e.g. \code{beta[i] ~ normal(0, 5);}) are \strong{not}
#' detected -- such parameters fall back to having no prior detected
#' (or an implicit uniform prior, if both bounds are hard-coded).
#'
#' If a parameter has no `~` sampling statement in the model block (e.g.
#' its prior is specified via `target += ..._lpdf(...)`, or it genuinely
#' has no prior), \code{dist}, \code{arg1} and \code{arg2} are set to NA.
#'
#' If a parameter has no hard bound on one or both sides, the
#' corresponding \code{v_min}/\code{v_max} display range is instead set to
#' the 2.5% / 97.5% quantile of its prior distribution (i.e. the range
#' encompassing the central 95% of the distribution). If no prior is
#' available either, \code{v_min}/\code{v_max} are left as NA and should
#' be filled in manually before plotting.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' stan_code <- "
#' data {
#'   real mu_prior;
#' }
#' parameters {
#'   real<lower=0, upper=5> A;
#'   real<lower=0.01> B;
#'   real C;
#' }
#' model {
#'   A ~ normal(0, 3);
#'   B ~ exponential(0.5);
#'   C ~ normal(mu_prior, 1);
#' }
#' "
#' df_priors <- Create_df_priors(stan_code, data_list = list(mu_prior = 0))
#'
#' # equivalently, from a compiled stanmodel object, e.g. produced by
#' # a ```{stan output.var = "wounds_model_2"} code chunk:
#' df_priors2 <- Create_df_priors(wounds_model_2, data_list = list(mu_prior = 0))
#' }
Create_df_priors <- function(stan_code, data_list = NULL) {

  # ---- read input ------------------------------------------------------
  # stan_code may be: (1) a compiled stanmodel S4 object (e.g. from
  # rstan::stan_model(), or from a ```{stan output.var = "x"} chunk),
  # (2) a path to a .stan file, or (3) a character string/vector of Stan
  # code itself.
  if (methods::is(stan_code, "stanmodel")) {
    code <- paste(methods::slot(stan_code, "model_code"), collapse = "\n")
  } else if (is.character(stan_code) && length(stan_code) == 1 &&
      !grepl("\\{", stan_code)) {
    # a single string with no "{" can't be real Stan code (every Stan
    # program has at least a data/parameters/model block), so treat it
    # as an intended file path -- and fail loudly if it doesn't exist,
    # rather than silently trying to parse the path string itself as code.
    if (!file.exists(stan_code)) {
      stop(
        "'stan_code' looks like a file path ('", stan_code, "') but no ",
        "such file exists relative to the current working directory ('",
        getwd(), "'). Check the path/filename, or use setwd()/here::here() ",
        "to point at the right location."
      )
    }
    code <- paste(readLines(stan_code, warn = FALSE), collapse = "\n")
  } else if (is.character(stan_code)) {
    code <- paste(stan_code, collapse = "\n")
  } else {
    stop(
      "'stan_code' must be a file path, a character string of Stan code, ",
      "or a compiled 'stanmodel' object (got class: ",
      paste(class(stan_code), collapse = "/"), ")."
    )
  }

  code <- strip_comments(code)

  # ---- extract blocks ----------------------------------------------------
  parameters_block             <- extract_block(code, "parameters")
  transformed_parameters_block <- extract_block(code, "transformed parameters")
  model_block                  <- extract_block(code, "model")

  if (is.na(parameters_block)) {
    stop("Could not find a 'parameters' block in the supplied Stan code.")
  }
  if (is.na(model_block)) {
    warning(
      "Could not find a 'model' block; all parameters will be set to ",
      "have no prior (dist = NA)."
    )
    model_block <- ""
  }

  # ---- parameters + hard bounds -------------------------------------------
  df_params <- parse_parameters_block(parameters_block, data_list)

  if (!is.na(transformed_parameters_block)) {
    df_tp <- tryCatch(
      parse_parameters_block(transformed_parameters_block, data_list),
      error = function(e) NULL  # e.g. transformed parameters block has no
    )                           # supported declarations, only statements
    if (!is.null(df_tp)) {
      dup <- intersect(df_tp$par, df_params$par)
      if (length(dup) > 0) {
        warning(
          "Parameter(s) '", paste(dup, collapse = ", "), "' declared in ",
          "both 'parameters' and 'transformed parameters'; keeping the ",
          "'parameters' block version."
        )
        df_tp <- df_tp[!(df_tp$par %in% dup), ]
      }
      df_params <- rbind(df_params, df_tp)
    }
  }

  # ---- priors ----------------------------------------------------------
  # priors are matched by 'base' name (e.g. "beta"), since a Stan
  # vectorized sampling statement like `beta ~ normal(0, 5);` applies the
  # same prior to every expanded element ("beta[1]", "beta[2]", ...).
  priors <- parse_priors_block(model_block, unique(df_params$base))

  # ---- assemble one row per parameter (or per vector element) -------------
  rows <- lapply(seq_len(nrow(df_params)), function(i) {

    pname <- df_params$par[i]
    base  <- df_params$base[i]
    v_lwr <- df_params$v_lwr[i]
    v_upr <- df_params$v_upr[i]

    dist <- NA_character_
    arg1 <- NA_real_
    arg2 <- NA_real_

    if (!is.null(priors[[base]])) {
      dist_stan <- priors[[base]]$dist_stan
      args      <- priors[[base]]$args

      if (dist_stan %in% names(.dist_map)) {
        mapping <- .dist_map[[dist_stan]]
        dist <- mapping$dist

        if (length(args) >= 1 && nzchar(args[1])) {
          arg1 <- resolve_numeric(args[1], data_list, base)
        }
        if (mapping$nargs == 2 && length(args) >= 2 && nzchar(args[2])) {
          arg2 <- resolve_numeric(args[2], data_list, base)
        }
      } else {
        warning(
          "Unsupported prior distribution '", dist_stan, "' for parameter '",
          base, "'; setting its prior to NA. Supported: ",
          paste(names(.dist_map), collapse = ", ")
        )
      }
    }

    # ---- implicit uniform prior ------------------------------------------
    # If both hard bounds are present but no explicit prior was found,
    # Stan implicitly treats the parameter as uniform over that range --
    # represent that explicitly so PriorPosteriorPlotStan() draws a flat
    # prior instead of leaving the panel without a prior ribbon.
    if (is.na(dist) && !is.na(v_lwr) && !is.na(v_upr)) {
      dist <- "uniform"
      arg1 <- v_lwr
      arg2 <- v_upr
    }

    # ---- display range (v_min / v_max) ---------------------------------
    v_min <- v_lwr
    v_max <- v_upr

    if (is.na(v_min) && !is.na(dist) && !is.na(arg1)) {
      v_min <- quantile_fn(dist, 0.025, arg1, arg2)
    }
    if (is.na(v_max) && !is.na(dist) && !is.na(arg1)) {
      v_max <- quantile_fn(dist, 0.975, arg1, arg2)
    }

    data.frame(
      par   = pname,
      v_min = v_min,
      v_max = v_max,
      v_lwr = v_lwr,
      v_upr = v_upr,
      dist  = dist,
      arg1  = arg1,
      arg2  = arg2,
      stringsAsFactors = FALSE
    )
  })

  df_priors <- do.call(rbind, rows)
  rownames(df_priors) <- NULL

  tibble::as_tibble(df_priors)
}


#' Generate a prior-posterior plot directly from Stan code and a stanfit object
#'
#' @description
#' Convenience wrapper around \code{Create_df_priors()} and
#' \code{PriorPosteriorPlotStan()}: it parses \code{df_priors} directly from
#' the supplied Stan code and immediately produces the prior-posterior plot,
#' so there is no need to call \code{Create_df_priors()} separately first.
#'
#' @param stan_fit An object returned by rstan::sampling. This supplies the
#'   posterior samples; it is distinct from \code{stan_code}, which supplies
#'   the parameter bounds and priors and does not itself contain samples.
#' @param stan_code A file path to a .stan file, a character string/vector
#'   containing Stan code, or a compiled \code{stanmodel} object. See
#'   \code{Create_df_priors()} for details on what is parsed from this.
#' @param pars A vector of parameter names to plot. Entries may be exact
#'   \code{df_priors$par} values (e.g. a scalar parameter, or an already
#'   bracket-indexed element like \code{"etaR[1]"}), or the *base* name of
#'   a vector/matrix parameter (e.g. \code{"etaR"}, \code{"beta_genus"}),
#'   which expands to all of that parameter's elements automatically. If
#'   omitted (\code{NULL},
#'   the default), every parameter identified by \code{Create_df_priors()}
#'   from \code{stan_code} is plotted -- note this means only the scalar
#'   \code{real} parameters it recognises (see \code{Create_df_priors()}
#'   for the current limitations on supported parameter types).
#' @param data_list Optional named list (e.g. the same list passed to
#'   \code{rstan::sampling(data = ...)}), used to resolve prior arguments or
#'   parameter bounds that reference named constants rather than numeric
#'   literals. See \code{Create_df_priors()}.
#' @param ncol Number of columns provided to facet_wrap. Square arrangement
#'   is produced when no value is provided.
#' @param nbins Number of bins used to display histograms (default is 25).
#'
#' @return A ggplot object.
#'
#' @details
#' This function performs no plotting logic of its own -- it simply calls
#' \code{Create_df_priors(stan_code, data_list)} to build \code{df_priors},
#' then passes it straight to \code{PriorPosteriorPlotStan()}. Use
#' \code{Create_df_priors()} and \code{PriorPosteriorPlotStan()} directly
#' instead if you need to inspect or edit \code{df_priors} (e.g. to
#' override an auto-derived bound, or add a prior the parser couldn't
#' detect) before plotting.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pars <- c("ymin", "r_F", "r_M", "aStar_F", "aStar_M")
#' PriorPosteriorPlot(stan_fit, "wounds_2.stan", pars)
#'
#' # omit pars to plot every parameter Create_df_priors() could identify
#' PriorPosteriorPlot(stan_fit, "wounds_2.stan")
#' }
PriorPosteriorPlot <- function(stan_fit, stan_code, pars = NULL, data_list = NULL,
                                ncol = NA, nbins = 25) {
  df_priors <- Create_df_priors(stan_code, data_list = data_list)

  if (is.null(pars)) {
    pars <- as.character(df_priors$par)
  }

  PriorPosteriorPlotStan(stan_fit, pars, df_priors, ncol = ncol, nbins = nbins)
}
