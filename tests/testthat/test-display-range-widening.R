# Regression test for a real bug found 2026-08-26: v_min/v_max are derived
# purely from the prior (a hard Stan bound, or else the prior's own
# 2.5%/97.5% quantile), with no knowledge of where the posterior actually
# lands. PriorCurves() draws the ribbon only across [v_min, v_max], so
# whenever the data pulls the posterior outside that illustrative window -
# which is exactly the interesting case this whole plot exists to show -
# histogram bars appeared past the edge of the ribbon with no prior shading
# there, looking like the posterior had violated a bound it never actually
# violated. Confirmed against a real model: 7 of 13 checked parameters had
# posterior draws outside the old window, in one case by over 50% of the
# window's own width (sigma_block).

test_that("the prior ribbon's display range is widened to cover posterior draws that fall outside the prior's own window", {
  n_draws <- 500
  # deliberately shift the "posterior" well outside the prior's own
  # 2.5%/97.5% quantile window - normal(0, 0.5) has ~[-1, 1] as that window
  set.seed(1)
  mock_draws <- matrix(rnorm(n_draws, mean = 5, sd = 0.2), nrow = n_draws, dimnames = list(NULL, "a"))
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws)

  df_priors <- data.frame(
    par = "a", v_min = qnorm(0.025, 0, 0.5), v_max = qnorm(0.975, 0, 0.5),
    v_lwr = NA_real_, v_upr = NA_real_, dist = "normal", arg1 = 0, arg2 = 0.5
  )

  # the OLD window would not contain the shifted "posterior" at all
  expect_true(df_priors$v_max < min(mock_draws))

  p <- PriorPosteriorPlotStan(mock_fit, "a", df_priors)
  ribbon_layer <- p$layers[[1]]
  ribbon_x_range <- range(ribbon_layer$data$x)

  expect_true(ribbon_x_range[1] <= min(mock_draws))
  expect_true(ribbon_x_range[2] >= max(mock_draws))
})

test_that("a genuine hard-bound violation (which should never happen for a valid fit) is warned about, not silently widened past", {
  n_draws <- 50
  # simulate an impossible scenario: draws below a declared hard lower bound
  mock_draws <- matrix(c(-1, rep(1, n_draws - 1)), nrow = n_draws, dimnames = list(NULL, "a"))
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws)

  df_priors <- data.frame(
    par = "a", v_min = 0, v_max = 5, v_lwr = 0, v_upr = NA_real_,
    dist = "exponential", arg1 = 1, arg2 = NA_real_
  )

  expect_warning(
    PriorPosteriorPlotStan(mock_fit, "a", df_priors),
    "outside the hard bound"
  )
})

test_that("display-range widening does not trigger a spurious hard-bound warning for the normal, expected case", {
  n_draws <- 50
  mock_draws <- matrix(rexp(n_draws), nrow = n_draws, dimnames = list(NULL, "a"))
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws)

  df_priors <- data.frame(
    par = "a", v_min = 0, v_max = 3, v_lwr = 0, v_upr = NA_real_,
    dist = "exponential", arg1 = 1, arg2 = NA_real_
  )

  expect_no_warning(PriorPosteriorPlotStan(mock_fit, "a", df_priors))
})
