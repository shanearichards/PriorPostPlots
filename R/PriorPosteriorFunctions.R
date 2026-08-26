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

#' Marginal distribution of a single off-diagonal element of a
#' correlation matrix under an LKJ prior
#'
#' @description
#' Under Stan's \code{lkj_corr_cholesky(eta)} prior (equivalently
#' \code{lkj_corr(eta)} on the correlation matrix \code{R} itself) in
#' \code{K} dimensions, any single off-diagonal entry \code{R[i,j]}
#' marginally follows a shifted, scaled Beta distribution on \code{(-1,
#' 1)}: \code{(R[i,j] + 1) / 2 ~ Beta(alpha, alpha)} with
#' \code{alpha = eta + (K - 2) / 2} (Lewandowski, Kurowicka & Joe, 2009).
#' \code{dlkjcorr1()}/\code{plkjcorr1()}/\code{qlkjcorr1()}/\code{rlkjcorr1()}
#' are that marginal's density, CDF, quantile function, and random-draw
#' generator, expressed in terms of Stan's own \code{eta} and \code{K}
#' rather than the Beta distribution's \code{alpha} directly.
#'
#' @param x,q vector of correlation values in \code{(-1, 1)}.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param eta the \code{lkj_corr_cholesky()}/\code{lkj_corr()} shape
#'   parameter as declared in the Stan model (\code{eta = 1} is uniform
#'   over all valid correlation matrices; larger \code{eta} concentrates
#'   more mass near 0, i.e. favours weaker correlation).
#' @param K the dimension of the correlation matrix (\code{K >= 2}).
#'
#' @return Vector of densities (\code{dlkjcorr1}), probabilities
#'   (\code{plkjcorr1}), quantiles (\code{qlkjcorr1}), or random draws
#'   (\code{rlkjcorr1}).
#' @name lkjcorr1
NULL

#' @rdname lkjcorr1
#' @export
dlkjcorr1 <- function(x, eta, K) {
  alpha <- eta + (K - 2) / 2
  stats::dbeta((x + 1) / 2, alpha, alpha) / 2
}

#' @rdname lkjcorr1
#' @export
plkjcorr1 <- function(q, eta, K) {
  alpha <- eta + (K - 2) / 2
  stats::pbeta((q + 1) / 2, alpha, alpha)
}

#' @rdname lkjcorr1
#' @export
qlkjcorr1 <- function(p, eta, K) {
  alpha <- eta + (K - 2) / 2
  2 * stats::qbeta(p, alpha, alpha) - 1
}

#' @rdname lkjcorr1
#' @export
rlkjcorr1 <- function(n, eta, K) {
  alpha <- eta + (K - 2) / 2
  2 * stats::rbeta(n, alpha, alpha) - 1
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
        lkj_corr    = dlkjcorr1(x_vals, eta = arg1, K = arg2),
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
        lkj_corr = plkjcorr1(v_upr, eta = arg1, K = arg2) -
          plkjcorr1(v_lwr, eta = arg1, K = arg2),
        rep(NA, length(x_vals))
      )

      tibble::tibble(par = par, x = x_vals, y = y_vals, d = d) |>
        dplyr::mutate(f = .data$y / .data$d) |>
        dplyr::select(.data$par, .data$x, y = .data$f)
    }
  )

  return(df_curves)
}

#' Extract posterior draws from a fitted Stan model into the same
#' shape `rstan::extract()` produces - a named list with one entry per
#' base parameter: a plain vector for a scalar parameter, an
#' `iterations x J` matrix for a `vector[J]` parameter, or an
#' `iterations x R x C` array for a `matrix[R, C]` parameter.
#'
#' Two input types are handled: an `rstan::sampling()`/`rstan::stan()`
#' fit (class `"stanfit"`) is passed straight to `rstan::extract()`,
#' unchanged from previous versions of this package. Anything else with
#' a `$draws(variables = ..., format = "matrix")` method - which
#' includes every `cmdstanr` fit class (`CmdStanMCMC`, and more broadly
#' `CmdStanFit`) without `cmdstanr` needing to be an explicit dependency
#' of this package - has its `draws_matrix` output (column names
#' `"name"`/`"name[i]"`/`"name[i,j]"`, one column per scalar element)
#' reshaped into the same vector/matrix/array-per-parameter list
#' structure by hand, since `cmdstanr` has no equivalent of
#' `rstan::extract()`'s automatic per-parameter reshaping.
#' @keywords internal
extract_posterior_list <- function(stan_fit, base_pars) {
  if (inherits(stan_fit, "stanfit")) {
    return(rstan::extract(object = stan_fit, pars = base_pars))
  }

  has_draws_method <- tryCatch(is.function(stan_fit$draws), error = function(e) FALSE)
  if (!has_draws_method) {
    stop(
      "'stan_fit' must be an object of class 'stanfit' (from rstan::sampling()/",
      "rstan::stan()), or a cmdstanr fit object (or similar) with a '$draws()' ",
      "method (e.g. from cmdstanr::cmdstan_model()$sample()) -- got an object ",
      "of class: ", paste(class(stan_fit), collapse = "/"), "."
    )
  }

  draws_mat <- tryCatch(
    stan_fit$draws(variables = base_pars, format = "matrix"),
    error = function(e) {
      stop(
        "Failed to extract draws for parameter(s) '",
        paste(base_pars, collapse = "', '"), "' from 'stan_fit' via its ",
        "$draws() method: ", conditionMessage(e)
      )
    }
  )
  # Strip any draws_matrix/draws S3 classing immediately -- indexing such
  # an object is not always equivalent to indexing a plain matrix (e.g. it
  # can carry extra required attributes/metadata), and all that's needed
  # from here on is a plain numeric matrix with parameter names as columns.
  draws_mat <- matrix(as.numeric(draws_mat), nrow = nrow(draws_mat),
                       dimnames = list(NULL, colnames(draws_mat)))

  col_names <- colnames(draws_mat)
  n_iter <- nrow(draws_mat)
  out <- list()

  for (p in base_pars) {
    cols <- col_names[col_names == p | grepl(paste0("^", p, "\\["), col_names, perl = TRUE)]
    if (length(cols) == 0) next  # not present in this fit; silently omitted,

    if (length(cols) == 1 && cols == p) {
      out[[p]] <- as.numeric(draws_mat[, p])
      next
    }

    # "name[i]" (vector) or "name[i,j]" (matrix) columns - recover each
    # element's index/indices from its column name, matching how
    # parse_parameters_block() (Create_df_priors.R) names them.
    idx_strs <- sub(paste0("^", p, "\\[(.*)\\]$"), "\\1", cols)

    if (any(grepl(",", idx_strs, fixed = TRUE))) {
      idx <- do.call(rbind, lapply(strsplit(idx_strs, ",", fixed = TRUE), as.integer))
      arr <- array(NA_real_, dim = c(n_iter, max(idx[, 1]), max(idx[, 2])))
      for (k in seq_along(cols)) arr[, idx[k, 1], idx[k, 2]] <- draws_mat[, cols[k]]
      out[[p]] <- arr
    } else {
      idx <- as.integer(idx_strs)
      m <- matrix(NA_real_, nrow = n_iter, ncol = max(idx))
      for (k in seq_along(cols)) m[, idx[k]] <- draws_mat[, cols[k]]
      out[[p]] <- m
    }
  }

  out
}

#' Generates prior-posterior plots
#'
#' @description
#' Generates a ggplot comparing prior and posterior distributions by combining sampling output from a Stan model with a data frame specifying the prior distributions for selected model parameters.
#'
#' @param stan_fit A fitted Stan model: either an object returned by \code{rstan::sampling()}/\code{rstan::stan()} (class \code{"stanfit"}), or a \code{cmdstanr} fit object (e.g. from \code{cmdstanr::cmdstan_model()$sample()}, class \code{"CmdStanMCMC"}/\code{"CmdStanFit"}) or anything else exposing the same \code{$draws(variables = ..., format = "matrix")} method. See \code{extract_posterior_list()} for exactly how each is handled.
#' @param pars A vector of parameter names associated with the stan model for plotting. Entries may be exact \code{df_priors$par} values (e.g. a scalar parameter, or an already bracket-indexed element like \code{"etaR[1]"}), or the *base* name of a vector/matrix parameter (e.g. \code{"etaR"}, \code{"beta_genus"}), which expands to all of that parameter's elements as they appear in \code{df_priors} (e.g. \code{"etaR[1]"}, \code{"etaR[2]"}, \code{"etaR[3]"}) -- so elements don't need to be listed individually.
#' @param df_priors A data frame containing information that describes the prior for each parameter. The data frame must have the columns par, dist, arg1, arg2, v_min, v_max, v_lwr, v_upr. par = parameter name, dist = prior distribution (normal, log_normal, exponential, gamma, uniform, laplace, beta, lkj_corr), arg1 = first distribution parameter, arg2 = second distribution parameter (NA if not needed; for lkj_corr, arg1 is eta and arg2 is K, the correlation matrix's dimension - see ?lkjcorr1), v_min and v_max are the bounds of the plotted prior, v_lwr and v_upr are the stan-imposed parameter bounds (NA if none are set).
#' @param ncol Number of columns provided to facet_wrap. Square arrangement is produced when no value is provided.
#' @param nbins Number of bins used to display histograms (default is 25).
#'
#' @return A ggplot object.
#' @details
#' Currently, the distributions allowed are uniform, normal, log_normal, exponential,
#'  gamma, laplace, beta, and lkj_corr (the marginal distribution of a single pairwise
#'  correlation under an lkj_corr_cholesky() prior - see ?lkjcorr1 and
#'  Create_df_priors()'s own documentation for how these rows get built). See
#'  corresponding d-functions for arguments required.
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
#' PriorPosteriorPlotStan(stan_fit, pars, df_priors) # rstan fit
#'
#' # equally, a cmdstanr fit works directly, no conversion needed:
#' # fit <- cmdstanr::cmdstan_model("model.stan")$sample(data = stan_data)
#' # PriorPosteriorPlotStan(fit, pars, df_priors)
#' }
PriorPosteriorPlotStan <- function(stan_fit, pars, df_priors, ncol = NA, nbins = 25) {
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

  # Record which parameters have a GENUINE Stan hard bound before the
  # fallback below overwrites any NA with v_min/v_max (needed so
  # PriorCurves() can renormalize an unbounded density's display window as
  # if it were "bounded" there, for its own truncated-density integral) -
  # that fallback makes "real bound" vs. "just the display window"
  # indistinguishable from v_lwr/v_upr alone afterwards, but the later
  # widening step below needs exactly that distinction.
  df_priors_use$true_v_lwr <- df_priors_use$v_lwr
  df_priors_use$true_v_upr <- df_priors_use$v_upr

	df_priors_use$v_lwr <- pmax(df_priors_use$v_min, df_priors_use$v_lwr, na.rm=TRUE)
	df_priors_use$v_upr <- pmin(df_priors_use$v_max, df_priors_use$v_upr, na.rm=TRUE)

  # extract posterior samples and place in long format
  # rstan::extract()/extract_posterior_list() only accept base parameter
  # names (e.g. "beta"), not bracket-indexed names (e.g. "beta[1]") --
  # strip any "[...]" suffix before requesting the draws; LongParameters()
  # re-expands any resulting matrix-valued entries back into "name[i]"
  # columns.
  base_pars <- unique(sub("\\[.*\\]$", "", pars_use))
  l_posteriors <- extract_posterior_list(stan_fit, base_pars)
  df_posteriors_plot <- LongParameters(l_posteriors, pars_use)

  if (!(nrow(df_posteriors_plot) > 0)) { # at least 1 parameter is present
    stop("Error: missing prior information")
  }

  # ---- widen the display range to also cover the actual posterior draws --
  # v_min/v_max are derived purely from the PRIOR (a hard Stan bound if one
  # exists, else the prior's own 2.5%/97.5% quantile - see
  # Create_df_priors()) and know nothing about where the posterior actually
  # landed. Since the whole point of this plot is to show the data pulling
  # the posterior away from the prior, the posterior routinely extends
  # beyond that illustrative window - and PriorCurves() draws the ribbon
  # only across [v_min, v_max], so without this adjustment the histogram
  # would show bars past the edge of the ribbon with no prior shading
  # there at all, looking like the posterior violated a bound it never
  # actually violated. Widen (never narrow) each parameter's display range
  # to the union of its own [v_min, v_max] and its observed posterior
  # range.
  #
  # Widening to the raw draws' own min/max is not quite enough on its own:
  # geom_histogram()'s bin edges routinely extend past the raw data (its
  # default bin placement is not anchored to the data's exact min/max), so
  # a bin can still peek out past a ribbon that only reaches the raw
  # min/max. Precompute the histogram bars ourselves - from the posterior
  # draws ALONE, with no ribbon layer present - and widen to their actual
  # bin edges instead.
  #
  # This precomputation has to happen before the ribbon is drawn, and the
  # final plot below has to render these same precomputed bars (via
  # geom_rect(), not a second live geom_histogram() call) rather than
  # letting ggplot2 recompute the histogram itself. With
  # facet_wrap(scales = "free"), stat_bin()'s automatic bin-width
  # calculation reads each panel's fully-trained x-scale, which is the
  # union of every layer sharing that panel - so a live geom_histogram()
  # sitting alongside the ribbon would have its own bin edges pulled wider
  # by the very ribbon we widened to match it, which pulls the ribbon
  # wider still on the next look, and so on (confirmed empirically against
  # a real model: re-deriving bin edges from a two-layer probe grew them
  # on every iteration instead of settling down). Precomputing the bars
  # from a posterior-only probe and drawing them as static geom_rect() data
  # removes that feedback loop entirely, since nothing about the ribbon can
  # feed back into a computation that already happened.
  #
  # A genuine Stan hard bound (v_lwr/v_upr) is a different matter: a valid
  # posterior draw can never legitimately fall outside one (Stan's own
  # constraining transform enforces it), so if the RAW DRAWS do, something
  # is actually wrong (e.g. df_priors has the wrong bound attached to this
  # parameter) - flag it rather than silently widening past it. A
  # HISTOGRAM BIN edge, on the other hand, can extend past a true bound
  # purely as a placement artifact with no real draws actually out there
  # (stat_bin() has no notion of the parameter's valid support), so the
  # widened range is re-clamped to any true bound afterwards rather than
  # visually implying the parameter's support extends past where it really
  # can.
  post_range <- df_posteriors_plot |>
    dplyr::summarise(
      post_min = min(.data$val, na.rm = TRUE),
      post_max = max(.data$val, na.rm = TRUE),
      .by = "par"
    )

  hist_probe <- ggplot2::ggplot_build(
    ggplot2::ggplot(df_posteriors_plot, ggplot2::aes(x = .data$val)) +
      ggplot2::geom_histogram(bins = nbins) +
      ggplot2::facet_wrap(ggplot2::vars(par), scales = "free")
  )
  df_bins_plot <- dplyr::left_join(
    hist_probe$data[[1]], hist_probe$layout$layout, by = "PANEL"
  )
  bin_range <- df_bins_plot |>
    dplyr::summarise(bin_min = min(.data$xmin), bin_max = max(.data$xmax), .by = "par")

  df_priors_use <- df_priors_use |>
    dplyr::left_join(post_range, by = "par") |>
    dplyr::left_join(bin_range, by = "par") |>
    dplyr::mutate(
      out_of_bounds = (!is.na(.data$true_v_lwr) & .data$post_min < .data$true_v_lwr) |
        (!is.na(.data$true_v_upr) & .data$post_max > .data$true_v_upr),
      v_min = pmin(.data$v_min, .data$bin_min, na.rm = TRUE),
      v_max = pmax(.data$v_max, .data$bin_max, na.rm = TRUE),
      v_min = dplyr::if_else(!is.na(.data$true_v_lwr), pmax(.data$v_min, .data$true_v_lwr), .data$v_min),
      v_max = dplyr::if_else(!is.na(.data$true_v_upr), pmin(.data$v_max, .data$true_v_upr), .data$v_max)
    )

  if (any(df_priors_use$out_of_bounds, na.rm = TRUE)) {
    warning(
      "Posterior draws fall outside the hard bound(s) recorded in 'df_priors' ",
      "for: '", paste(df_priors_use$par[df_priors_use$out_of_bounds], collapse = "', '"),
      "' - this should not be possible for a genuine Stan <lower=, upper=> ",
      "constraint, so double-check df_priors has the right bounds attached ",
      "to the right parameter."
    )
  }

  # The bars themselves need the same true-bound clamp as the ribbon just
  # got, and for the same reason: a bin edge dipping past a true hard bound
  # is stat_bin()'s placement artifact, not a real draw out there, so the
  # bar shouldn't visually claim support past it either - otherwise the
  # ribbon and the bar would each be clamped to a *different* edge right at
  # the boundary, leaving a sliver of bar poking out with no ribbon under
  # it (the one situation the widening above doesn't otherwise reach).
  df_bins_plot <- df_bins_plot |>
    dplyr::left_join(
      dplyr::select(df_priors_use, "par", "true_v_lwr", "true_v_upr"), by = "par"
    ) |>
    dplyr::mutate(
      xmin = dplyr::if_else(!is.na(.data$true_v_lwr), pmax(.data$xmin, .data$true_v_lwr), .data$xmin),
      xmax = dplyr::if_else(!is.na(.data$true_v_upr), pmin(.data$xmax, .data$true_v_upr), .data$xmax)
    )

  df_priors_use <- dplyr::select(
    df_priors_use,
    -"post_min", -"post_max", -"bin_min", -"bin_max", -"out_of_bounds", -"true_v_lwr", -"true_v_upr"
  )

  # generate data for prior distribution curves in long format
  df_priors_plot <- PriorCurves(df_priors_use)

  # set the number of facet columns
  ncol <- ifelse(is.na(ncol), ceiling(0.5*length(pars_use)), ncol)

  # df_bins_plot$par came off hist_probe's own facet layout, so it needs
  # the same factor levels as df_priors_plot/df_posteriors_plot to facet
  # into matching, correctly-ordered panels below.
  df_bins_plot$par <- factor(df_bins_plot$par, levels = pars_use)

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = df_priors_plot,
      ggplot2::aes(x = .data$x, ymin = 0, ymax = .data$y), fill = "salmon") +
    ggplot2::geom_rect(data = df_bins_plot,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = .data$density),
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

#' Try to resolve a Stan expression (e.g. "0.5", "sigma_prior", or
#' "log(0.5)") to a number: first as a bare numeric literal, then as a
#' name looked up in `data_list`, and finally - for anything else, e.g. a
#' function call or arithmetic expression using syntax common to both
#' Stan and R (`log(0.5)`, `1/2`, `sqrt(2)`) - by evaluating it as an R
#' expression, with any `data_list` entries available as named values
#' inside that expression too (e.g. `normal(0, 2 * sigma_prior)`).
#' @keywords internal
resolve_numeric <- function(x, data_list = NULL, par_name = NULL) {
  x <- trimws(x)
  val <- suppressWarnings(as.numeric(x))
  if (!is.na(val)) return(val)

  if (!is.null(data_list) && x %in% names(data_list)) {
    val <- suppressWarnings(as.numeric(data_list[[x]]))
    if (!is.na(val)) return(val)
  }

  val <- tryCatch({
    env <- list2env(if (is.null(data_list)) list() else data_list, parent = baseenv())
    result <- eval(parse(text = x), envir = env)
    if (is.numeric(result) && length(result) == 1) result else NA_real_
  }, error = function(e) NA_real_)
  if (!is.na(val)) return(val)

  warning(
    "Could not resolve value '", x, "'",
    if (!is.null(par_name)) paste0(" for parameter '", par_name, "'") else "",
    " to a number (it is not a literal, a 'data_list' entry, or a ",
    "resolvable numeric expression). Setting to NA."
  )
  NA_real_
}

#' Parse the parameters block: identify scalar `real`, `vector[N]` and
#' `matrix[R, C]` parameters, and any hard lower/upper bounds declared on
#' them.
#' @keywords internal
parse_parameters_block <- function(block, data_list = NULL, pars_filter = NULL, skip_names = NULL) {
  statements <- strsplit(block, ";")[[1]]
  statements <- trimws(gsub("\\s+", " ", statements))
  statements <- statements[statements != ""]

  # cholesky_factor_corr is deliberately NOT in this list: it's recognised
  # and silently skipped here (no row, no warning) because it's handled by
  # a separate pass, parse_lkj_priors(), which turns its
  # lkj_corr_cholesky(eta) prior into one row per pairwise correlation
  # instead of one row for the raw Cholesky factor itself (see that
  # function's documentation for why).
  unsupported_types <- c("array", "int", "simplex",
    "ordered", "positive_ordered", "row_vector",
    "cholesky_factor_cov", "corr_matrix", "cov_matrix", "unit_vector")

  # When the caller only wants a handful of parameters (e.g.
  # PriorPosteriorPlotStan()'s `pars`), a statement whose declared name(s)
  # aren't in that set is skipped entirely, before any size-resolution or
  # type-support check runs - so a model's other parameters (random-effect
  # z-scores, correlation Cholesky factors, etc.) generate no warnings at
  # all rather than warnings about parameters nobody asked to plot.
  #
  # `skip_names` is a second, unconditional skip list (irrespective of
  # `pars_filter`): a `corr_matrix[K]` declared inside `transformed
  # parameters` (rather than `generated quantities`) would otherwise be
  # flagged here as an "unsupported type" even when it's one
  # parse_lkj_priors() has already successfully matched to an
  # lkj_corr_cholesky() prior and IS being plotted, just via that separate
  # pass rather than this one.
  wanted <- function(candidate_names) {
    (is.null(pars_filter) || any(candidate_names %in% pars_filter)) &&
      !any(candidate_names %in% skip_names)
  }

  # Generic name extractor for declaration shapes this function doesn't
  # otherwise parse (the `unsupported_types` catch-all covers many - array,
  # simplex, ordered, corr_matrix, etc. - each with a different
  # bracket/constraint layout). Stan always closes any <...>/[...] type
  # annotation before listing variable names, so the name list is
  # whatever's left after the last "]" (if present), else the last ">"
  # (if present), else the leading type keyword.
  extract_declared_names <- function(stmt) {
    after <- if (grepl("\\]", stmt)) {
      sub("^.*\\]\\s*", "", stmt)
    } else if (grepl(">", stmt)) {
      sub("^.*>\\s*", "", stmt)
    } else {
      sub("^[A-Za-z_][A-Za-z0-9_]*\\s+", "", stmt)
    }
    trimws(strsplit(sub("=.*$", "", after), ",")[[1]])
  }

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
      if (!wanted(trimws(strsplit(names_str, ",")[[1]]))) next

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
      if (!wanted(trimws(strsplit(names_str, ",")[[1]]))) next

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
      if (!wanted(trimws(strsplit(names_str, ",")[[1]]))) next

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

    if (any(startsWith(stmt, unsupported_types)) && wanted(extract_declared_names(stmt))) {
      warning(
        "Parameter declaration '", stmt, "' has an unsupported type; ",
        "only scalar 'real', 'vector' and 'matrix' parameters are ",
        "currently handled by Create_df_priors(). Skipping."
      )
    }
  }

  if (length(results) == 0) {
    # Not necessarily an error - a block can legitimately contain only a
    # cholesky_factor_corr parameter (handled separately by
    # parse_lkj_priors(), not here). Create_df_priors() checks the
    # combined result of both passes for real emptiness.
    return(data.frame(
      par = character(0), base = character(0),
      v_lwr = numeric(0), v_upr = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  df <- do.call(rbind, results)
  rownames(df) <- NULL
  df
}

#' Find the position of the `)` matching the `(` at `open_pos` in `text`,
#' by scanning forward and tracking parenthesis depth. Handles arbitrary
#' nesting (e.g. the `)` closing `normal(log(0.5), 0.75)`'s outer call).
#' Returns `NA` if the parentheses from `open_pos` onward are unbalanced.
#' @keywords internal
find_matching_paren <- function(text, open_pos) {
  chars <- strsplit(substr(text, open_pos, nchar(text)), "")[[1]]
  depth <- 0L
  for (i in seq_along(chars)) {
    if (chars[i] == "(") {
      depth <- depth + 1L
    } else if (chars[i] == ")") {
      depth <- depth - 1L
      if (depth == 0L) return(open_pos + i - 1L)
    }
  }
  NA_integer_
}

#' Split a Stan sampling-statement argument list on top-level commas only
#' (i.e. commas not nested inside a further `(...)`), so an argument that
#' is itself a function call with its own comma-separated arguments (e.g.
#' a 2-argument nested call) is kept intact rather than being split at
#' that inner comma.
#' @keywords internal
split_args_respecting_parens <- function(args_str) {
  chars <- strsplit(args_str, "")[[1]]
  depth <- 0L
  piece_start <- 1L
  pieces <- character(0)
  for (i in seq_along(chars)) {
    if (chars[i] == "(") {
      depth <- depth + 1L
    } else if (chars[i] == ")") {
      depth <- depth - 1L
    } else if (chars[i] == "," && depth == 0L) {
      pieces <- c(pieces, substr(args_str, piece_start, i - 1L))
      piece_start <- i + 1L
    }
  }
  pieces <- c(pieces, substr(args_str, piece_start, nchar(args_str)))
  trimws(pieces)
}

#' Parse the model block for `par ~ dist(args);` sampling statements
#' associated with the supplied parameter names. Recognises three LHS
#' forms: a bare vectorized statement (`par ~ dist(args);`, applied to
#' every expanded element of `par`); the `to_vector(par) ~ dist(args);`
#' idiom commonly used to put a vectorized prior on a whole `matrix`
#' parameter (handled the same way as a bare vectorized statement); and a
#' per-element statement with a literal index (`par[1] ~ dist(args);` for
#' a `vector` element, `par[1,2] ~ dist(args);` for a `matrix` element) -
#' if both a per-element and a vectorized statement exist for the same
#' base parameter, the per-element one takes precedence for that specific
#' element (resolved in `Create_df_priors()`, which looks up each
#' expanded element's own indexed name before falling back to its base
#' name).
#'
#' The argument list itself is extracted with `find_matching_paren()`
#' rather than a single-level `"\\\\(([^)]*)\\\\)"` regex capture, and then
#' split into individual arguments with `split_args_respecting_parens()`
#' rather than a plain `strsplit(..., ",")` - both are parenthesis-depth
#' aware, so an argument that is itself a function call (e.g. the
#' `log(0.5)` in `normal(log(0.5), 0.75)`) doesn't truncate/corrupt the
#' rest of the argument list (a plain, depth-unaware capture stops at the
#' FIRST `)`, which closes the nested call rather than the sampling
#' statement itself).
#'
#' Note: this does not (yet) detect priors specified via target +=
#' *_lpdf(...)/*_lpmf(...) syntax, priors written with a `for`-loop
#' variable as the index (`for (i in 1:N) par[i] ~ normal(0, sigma);`  -
#' only a *literal* index like `par[1]` is recognised, not a loop
#' variable), nor does it adjust for truncation (`T[lower, upper]`) -
#' these fall back to "no prior" (dist = NA).
#' @keywords internal
parse_priors_block <- function(block, par_names) {
  # Matches the LHS + distribution name + opening "(" of a sampling
  # statement only, stopping right after that "(" - see the argument-
  # extraction comment above for why the argument list itself is found
  # separately, not as part of this regex.
  pattern <- paste0(
    "(?:to_vector\\(\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\)|",
    "([A-Za-z_][A-Za-z0-9_]*)(\\[\\s*[0-9]+\\s*(?:,\\s*[0-9]+\\s*)?\\])?)",
    "\\s*~\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\("
  )
  m <- gregexpr(pattern, block, perl = TRUE)
  match_starts <- m[[1]]
  if (match_starts[1] == -1) return(list())
  match_lengths <- attr(m[[1]], "match.length")

  priors <- list()

  for (k in seq_along(match_starts)) {
    mt <- substr(block, match_starts[k], match_starts[k] + match_lengths[k] - 1)
    parts <- regmatches(mt, regexec(pattern, mt, perl = TRUE))[[1]]

    if (nzchar(parts[2])) {
      # to_vector(par) ~ dist(...) - applies to every element of `par`
      pname <- parts[2]
      base  <- parts[2]
    } else {
      base  <- parts[3]
      index <- gsub("\\s+", "", parts[4])  # "" (bare/vectorized), "[i]" or "[i,j]"
      pname <- paste0(base, index)
    }
    dist_stan <- parts[5]

    if (!(base %in% par_names)) next

    # The match above ends right after the distribution's opening "(";
    # walk forward from there tracking parenthesis depth to find its
    # TRUE matching close, so e.g. "log(0.5)" as an argument survives
    # intact instead of truncating the match at its own closing paren.
    open_pos  <- match_starts[k] + match_lengths[k] - 1
    close_pos <- find_matching_paren(block, open_pos)
    if (is.na(close_pos)) next  # unbalanced parens - malformed Stan code
    args_str <- substr(block, open_pos + 1, close_pos - 1)
    args <- split_args_respecting_parens(args_str)

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

#' Scan the `transformed parameters`/`generated quantities` blocks for the
#' standard Stan idiom that materialises a correlation matrix from a
#' `cholesky_factor_corr` parameter: `<name> = multiply_lower_tri_self_
#' transpose(<L name>);` (with an optional leading type/size declaration).
#' Returns a named list keyed by `<L name>`, valued `<name>` - used both to
#' build `parse_lkj_priors()`'s prior rows and, by `Create_df_priors()`, to
#' keep `parse_parameters_block()` from separately flagging `<name>` as an
#' "unsupported type" (`corr_matrix`) when it's declared inside
#' `transformed parameters` rather than `generated quantities`.
#' @keywords internal
find_lkj_corr_map <- function(transformed_parameters_block, generated_quantities_block) {
  gq_code <- paste(
    c(if (!is.na(transformed_parameters_block)) transformed_parameters_block,
      if (!is.na(generated_quantities_block)) generated_quantities_block),
    collapse = "\n"
  )
  gq_statements <- trimws(gsub("\\s+", " ", strsplit(gq_code, ";")[[1]]))
  gq_statements <- gq_statements[gq_statements != ""]

  corr_pattern <- paste0(
    "(?:^[A-Za-z_][A-Za-z0-9_]*(?:<[^>]*>)?\\s*\\[[^\\]]*\\]\\s+)?",
    "([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*multiply_lower_tri_self_transpose",
    "\\s*\\(\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\)"
  )

  corr_map <- list()
  for (stmt in gq_statements) {
    m <- regmatches(stmt, regexec(corr_pattern, stmt, perl = TRUE))[[1]]
    if (length(m) == 0) next
    corr_map[[m[3]]] <- m[2]
  }
  corr_map
}

#' Turn `cholesky_factor_corr[K] L; ... L ~ lkj_corr_cholesky(eta);` into
#' one prior row per pairwise correlation.
#'
#' @description
#' Unlike every other supported parameter type, the raw elements of a
#' `cholesky_factor_corr[K]` parameter are not individually meaningful:
#' `L[1,1]` is a hard-coded constant 1, and (for `K > 2`) the remaining
#' entries are components of a spherical/stick-breaking parameterisation,
#' not correlations themselves. The quantities worth plotting are the
#' off-diagonal entries of the correlation matrix `R = L %*% t(L)`, each
#' individually bounded in `(-1, 1)` and following a shifted-Beta marginal
#' distribution under `L ~ lkj_corr_cholesky(eta)` - see `?lkjcorr1`.
#'
#' `R` itself has to come from somewhere in the *posterior*, though - a
#' `cholesky_factor_corr[K]` Stan parameter's own draws are `L`, not `R`.
#' This looks for the standard Stan idiom for materialising `R` from `L`:
#' `corr_matrix[K] Corr_L = multiply_lower_tri_self_transpose(L);` (an
#' optional leading type/size declaration, in case `Corr_L` was already
#' declared elsewhere and only assigned here), in the `transformed
#' parameters` or `generated quantities` block. If a matching statement is
#' found, one prior row is emitted per pairwise correlation, named to
#' match that generated quantity's own draws (e.g. `"Corr_L[1,2]"`) so it
#' flows through `PriorPosteriorPlotStan()`'s existing matrix-element
#' machinery unchanged - `Corr_L` is just an ordinary `matrix[K,K]`-shaped
#' generated quantity as far as the rest of the package is concerned. If
#' no matching statement is found, a warning suggests one to add (skipped
#' entirely, not warned about, for any `L` not in `pars_filter`).
#'
#' Only searches the `transformed parameters`/`generated quantities`
#' blocks, and only recognises that one idiom - a correlation matrix
#' declared in the `parameters` block instead of computed from `L`, or
#' materialised via some other expression (e.g. a manual `L * L'`), is not
#' detected.
#' @keywords internal
parse_lkj_priors <- function(parameters_block, model_block,
                              transformed_parameters_block,
                              generated_quantities_block,
                              data_list, pars_filter) {

  # ---- find cholesky_factor_corr[K] name(s) ------------------------------
  statements <- strsplit(parameters_block, ";")[[1]]
  statements <- trimws(gsub("\\s+", " ", statements))
  statements <- statements[statements != ""]

  l_names <- character(0)
  l_sizes <- list()

  for (stmt in statements) {
    m <- regmatches(stmt, regexec(
      "^cholesky_factor_corr\\s*\\[\\s*([^\\]]+)\\s*\\]\\s+(.+)$", stmt,
      perl = TRUE))[[1]]
    if (length(m) == 0) next

    size_expr <- trimws(m[2])
    names_str <- sub("=.*$", "", m[3])
    k_val <- resolve_numeric(size_expr, data_list)
    if (is.na(k_val)) {
      warning(
        "Could not resolve size '", size_expr, "' for cholesky_factor_corr ",
        "parameter(s) '", trimws(names_str), "'; pass its value via ",
        "'data_list' if it's a data-block constant. Skipping its ",
        "pairwise-correlation prior(s)."
      )
      next
    }
    k_val <- as.integer(round(k_val))

    for (nm in trimws(strsplit(names_str, ",")[[1]])) {
      if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", nm)) next
      l_names <- c(l_names, nm)
      l_sizes[[nm]] <- k_val
    }
  }

  if (length(l_names) == 0) return(NULL)

  # ---- find each L's lkj_corr_cholesky(eta) sampling statement -----------
  lkj_priors <- parse_priors_block(model_block, l_names)

  # ---- auto-detect each L's companion correlation-matrix statement -------
  corr_map <- find_lkj_corr_map(transformed_parameters_block, generated_quantities_block)

  # ---- assemble one row per pairwise correlation -------------------------
  rows <- list()

  for (l_name in l_names) {
    prior_entry <- lkj_priors[[l_name]]
    if (is.null(prior_entry)) next  # no sampling statement found; no prior

    if (prior_entry$dist_stan != "lkj_corr_cholesky") {
      warning(
        "Unsupported prior distribution '", prior_entry$dist_stan,
        "' for cholesky_factor_corr parameter '", l_name, "'; only ",
        "lkj_corr_cholesky(eta) is currently supported. Skipping."
      )
      next
    }

    eta   <- resolve_numeric(prior_entry$args[1], data_list, l_name)
    k_val <- l_sizes[[l_name]]
    corr_name <- corr_map[[l_name]]

    if (is.null(corr_name)) {
      if (is.null(pars_filter) || l_name %in% pars_filter) {
        suggested <- if (grepl("^L_", l_name)) sub("^L_", "Corr_", l_name) else paste0("Corr_", l_name)
        warning(
          "'", l_name, "' has an lkj_corr_cholesky(", eta, ") prior, but no ",
          "generated quantity of the form '<name> = ",
          "multiply_lower_tri_self_transpose(", l_name, ")' was found in the ",
          "Stan code, so its pairwise correlation(s) can't be matched to ",
          "posterior draws and will be skipped. Add a line such as:\n    ",
          "corr_matrix[", k_val, "] ", suggested,
          " = multiply_lower_tri_self_transpose(", l_name, ");\n",
          "to your model's 'transformed parameters' or 'generated quantities' ",
          "block to enable this check."
        )
      }
      next
    }

    if (!(is.null(pars_filter) || corr_name %in% pars_filter || l_name %in% pars_filter)) next

    for (i in seq_len(k_val - 1)) {
      for (j in seq(i + 1, k_val)) {
        rows[[length(rows) + 1]] <- data.frame(
          par   = paste0(corr_name, "[", i, ",", j, "]"),
          v_min = -1, v_max = 1, v_lwr = -1, v_upr = 1,
          dist  = "lkj_corr", arg1 = eta, arg2 = k_val,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
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
    lkj_corr    = qlkjcorr1(p, eta = arg1, K = arg2),
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
#' @param pars Optional vector of parameter names. If supplied, only
#'   declarations matching these (exact matches, or the *base* name of a
#'   vector/matrix parameter, e.g. \code{"beta"} for \code{"beta[1]"}) are
#'   parsed - every other parameter in the model is skipped silently,
#'   generating no warning even if its type or prior isn't one
#'   \code{Create_df_priors()} can handle. Leave as \code{NULL} (the
#'   default) to parse and report on every parameter in the model, e.g.
#'   when building a \code{df_priors} you intend to reuse across several
#'   different \code{pars} subsets later.
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
#' applied to every element of that vector. A per-element statement with a
#' \emph{literal} index (e.g. \code{beta[1] ~ normal(0.4, 1); beta[2] ~
#' normal(-0.9, 1);}) is also detected and takes precedence, for that
#' element, over any vectorized statement on the same parameter. A
#' statement indexed by a \code{for}-loop variable (e.g. \code{for (i in
#' 1:N) beta[i] ~ normal(0, sigma);}) is \strong{not} detected -- such
#' parameters fall back to having no prior detected (or an implicit
#' uniform prior, if both bounds are hard-coded). Sampling-statement
#' arguments may themselves contain nested function calls or further
#' parentheses (e.g. \code{normal(log(0.5), 0.75)}); these are resolved
#' correctly (\code{log(0.5)} is evaluated as an R expression once
#' extracted intact) rather than truncated at the first closing
#' parenthesis.
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
#' A \code{cholesky_factor_corr[K]} parameter \code{L} with an
#' \code{L ~ lkj_corr_cholesky(eta);} prior is handled differently from
#' every other type: its own raw elements aren't individually meaningful
#' (see \code{?lkjcorr1}), so instead of a row for \code{L} itself, one row
#' is emitted per pairwise correlation - \code{K*(K-1)/2} of them - each
#' with \code{dist = "lkj_corr"}, \code{v_lwr}/\code{v_upr} = -1/1, and
#' named to match a companion correlation-matrix generated quantity
#' (\code{corr_matrix[K] Corr_L = multiply_lower_tri_self_transpose(L);}),
#' auto-detected from the \code{transformed parameters}/
#' \code{generated quantities} block. If \code{L} has an
#' \code{lkj_corr_cholesky()} prior but no such generated quantity is
#' found, a warning suggests the exact line to add; nothing is plotted for
#' it in the meantime, since Stan itself never samples anything on the
#' interpretable correlation scale in that case.
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
Create_df_priors <- function(stan_code, data_list = NULL, pars = NULL) {

  # A caller who only wants a handful of parameters (the common case: a
  # `pars` subset headed straight for PriorPosteriorPlotStan()/
  # PriorPosteriorPlot()) shouldn't have to see warnings about every OTHER
  # parameter in the model that Create_df_priors() can't handle (random-
  # effect z-scores, correlation Cholesky factors, etc.) - strip any
  # bracket index first, since `pars` may name either a base vector/matrix
  # parameter ("beta") or one specific element ("beta[1]"), but the
  # declaration-level filtering below only ever sees the base name.
  pars_filter <- if (!is.null(pars)) unique(sub("\\[.*\\]$", "", pars)) else NULL

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
  generated_quantities_block   <- extract_block(code, "generated quantities")

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
  # Computed up front (not just inside parse_lkj_priors(), later) so a
  # corr_matrix[K] generated quantity living in 'transformed parameters'
  # rather than 'generated quantities' doesn't ALSO get flagged as an
  # "unsupported type" by the ordinary parse_parameters_block() calls below
  # - it's already accounted for, just via parse_lkj_priors()'s own pass.
  lkj_corr_names <- unname(unlist(
    find_lkj_corr_map(transformed_parameters_block, generated_quantities_block)
  ))

  df_params <- parse_parameters_block(parameters_block, data_list, pars_filter, lkj_corr_names)

  if (!is.na(transformed_parameters_block)) {
    df_tp <- tryCatch(
      parse_parameters_block(transformed_parameters_block, data_list, pars_filter, lkj_corr_names),
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
  # priors are returned keyed by whatever LHS parse_priors_block() found:
  # a 'base' name (e.g. "beta") for a vectorized statement
  # (`beta ~ normal(0, 5);`, applying the same prior to every expanded
  # element "beta[1]", "beta[2]", ...), or a specific expanded element's
  # own name (e.g. "beta[1]") for a per-element statement
  # (`beta[1] ~ normal(...);`). Each row below looks up its own exact
  # name first, falling back to the base name.
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

    prior_entry <- priors[[pname]]
    if (is.null(prior_entry)) prior_entry <- priors[[base]]

    if (!is.null(prior_entry)) {
      dist_stan <- prior_entry$dist_stan
      args      <- prior_entry$args

      if (dist_stan %in% names(.dist_map)) {
        mapping <- .dist_map[[dist_stan]]
        dist <- mapping$dist

        if (length(args) >= 1 && nzchar(args[1])) {
          arg1 <- resolve_numeric(args[1], data_list, pname)
        }
        if (mapping$nargs == 2 && length(args) >= 2 && nzchar(args[2])) {
          arg2 <- resolve_numeric(args[2], data_list, pname)
        }
      } else {
        warning(
          "Unsupported prior distribution '", dist_stan, "' for parameter '",
          pname, "'; setting its prior to NA. Supported: ",
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

  # ---- lkj_corr_cholesky() parameters: one row per pairwise correlation --
  # A completely separate pass from the assembly loop above - see
  # parse_lkj_priors()'s own documentation for why a cholesky_factor_corr
  # parameter can't be handled the same way as every other type.
  df_lkj <- parse_lkj_priors(
    parameters_block, model_block, transformed_parameters_block,
    generated_quantities_block, data_list, pars_filter
  )
  if (!is.null(df_lkj)) {
    df_priors <- rbind(df_priors, df_lkj)
  }

  if (is.null(df_priors) || nrow(df_priors) == 0) {
    stop(
      "No supported parameters ('real', 'vector', 'matrix', or ",
      "'cholesky_factor_corr' with a matching correlation-matrix generated ",
      "quantity) were found in the supplied Stan code."
    )
  }

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
#' @param stan_fit A fitted Stan model: either an object returned by
#'   \code{rstan::sampling()}/\code{rstan::stan()}, or a \code{cmdstanr} fit
#'   object (or anything else exposing a \code{$draws()} method - see
#'   \code{extract_posterior_list()}). This supplies the posterior samples;
#'   it is distinct from \code{stan_code}, which supplies the parameter
#'   bounds and priors and does not itself contain samples.
#' @param stan_code A file path to a .stan file, a character string/vector
#'   containing Stan code, or a compiled \code{stanmodel} object. See
#'   \code{Create_df_priors()} for details on what is parsed from this.
#' @param pars A vector of parameter names to plot. Entries may be exact
#'   \code{df_priors$par} values (e.g. a scalar parameter, or an already
#'   bracket-indexed element like \code{"etaR[1]"}), or the *base* name of
#'   a vector/matrix parameter (e.g. \code{"etaR"}, \code{"beta_genus"}),
#'   which expands to all of that parameter's elements automatically. This
#'   is also forwarded into \code{Create_df_priors()}, so any OTHER
#'   parameter in \code{stan_code} - one this call was never going to plot -
#'   generates no warning even if it has an unsupported type or prior; see
#'   \code{Create_df_priors()}'s own \code{pars} argument. If omitted
#'   (\code{NULL}, the default), every parameter identified by
#'   \code{Create_df_priors()} from \code{stan_code} is plotted (and
#'   warned about, if unsupported) -- note this means only the parameter
#'   types \code{Create_df_priors()} recognises (see its documentation for
#'   the current limitations).
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
  # `pars` is forwarded into Create_df_priors() itself (not just applied
  # as a filter afterwards) so that any OTHER parameter in the model - one
  # nobody asked to plot - generates no warning at all, rather than a
  # warning about a parameter this call was never going to show.
  df_priors <- Create_df_priors(stan_code, data_list = data_list, pars = pars)

  if (is.null(pars)) {
    pars <- as.character(df_priors$par)
  }

  PriorPosteriorPlotStan(stan_fit, pars, df_priors, ncol = ncol, nbins = nbins)
}
