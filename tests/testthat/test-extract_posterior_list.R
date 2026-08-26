# Regression tests for accepting cmdstanr fits (or anything else exposing a
# $draws(variables, format) method) alongside rstan::stanfit objects, added
# 2026-08-26. cmdstanr fits are mocked here (a plain list with a $draws
# function) rather than fitted for real, so these tests don't need CmdStan
# installed and stay fast; the rstan dispatch path is mocked via
# local_mocked_bindings() for the same reason (no need to compile a real
# Stan model just to exercise the dispatch logic).

test_that("a cmdstanr-like fit ($draws method) is reshaped into rstan::extract()'s list shape", {
  n_draws <- 100
  mock_draws <- matrix(
    c(rnorm(n_draws), rnorm(n_draws), abs(rnorm(n_draws))),
    nrow = n_draws,
    dimnames = list(NULL, c("a[1]", "a[2]", "sigma"))
  )
  class(mock_draws) <- c("draws_matrix", "draws", "matrix", "array")  # mimic posterior's class

  mock_fit <- list(draws = function(variables, format = "matrix") {
    keep <- grepl(paste0("^(", paste(variables, collapse = "|"), ")(\\[|$)"), colnames(mock_draws))
    mock_draws[, keep, drop = FALSE]
  })

  out <- extract_posterior_list(mock_fit, c("a", "sigma"))

  expect_true(is.matrix(out$a))
  expect_false(inherits(out$a, "draws_matrix"))  # class stripped, see function docs
  expect_equal(dim(out$a), c(n_draws, 2))
  expect_equal(out$a[, 1], as.numeric(mock_draws[, "a[1]"]))
  expect_equal(out$a[, 2], as.numeric(mock_draws[, "a[2]"]))
  expect_equal(out$sigma, as.numeric(mock_draws[, "sigma"]))
})

test_that("a matrix parameter's '[i,j]' columns are reshaped into an iterations x R x C array", {
  n_draws <- 50
  cols <- c("Sigma[1,1]", "Sigma[2,1]", "Sigma[1,2]", "Sigma[2,2]")
  mock_draws <- matrix(rnorm(n_draws * 4), nrow = n_draws, dimnames = list(NULL, cols))
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws)

  out <- extract_posterior_list(mock_fit, "Sigma")

  expect_equal(dim(out$Sigma), c(n_draws, 2, 2))
  expect_equal(out$Sigma[, 1, 1], mock_draws[, "Sigma[1,1]"])
  expect_equal(out$Sigma[, 2, 1], mock_draws[, "Sigma[2,1]"])
  expect_equal(out$Sigma[, 1, 2], mock_draws[, "Sigma[1,2]"])
})

test_that("a genuine 'stanfit'-classed object is still dispatched to rstan::extract() (no regression)", {
  fake_stanfit <- structure(list(), class = "stanfit")
  local_mocked_bindings(extract = function(object, pars) list(mocked = TRUE), .package = "rstan")

  expect_equal(extract_posterior_list(fake_stanfit, "a"), list(mocked = TRUE))
})

test_that("an object with neither class gives an informative error", {
  expect_error(
    extract_posterior_list(list(x = 1), "a"),
    "must be an object of class 'stanfit'.*cmdstanr fit object"
  )
})

test_that("PriorPosteriorPlotStan() accepts a cmdstanr-like fit end-to-end", {
  n_draws <- 200
  mock_draws <- matrix(rnorm(n_draws), nrow = n_draws, dimnames = list(NULL, "a"))
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws)

  df_priors <- data.frame(
    par = "a", v_min = -3, v_max = 3, v_lwr = NA, v_upr = NA,
    dist = "normal", arg1 = 0, arg2 = 1
  )

  p <- PriorPosteriorPlotStan(mock_fit, "a", df_priors)
  expect_s3_class(p, "ggplot")
})
