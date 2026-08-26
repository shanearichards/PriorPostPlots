# Tests for inferring an unresolved vector/matrix size from a fitted
# model's own posterior draws, as a fallback when the size can't be
# resolved from a literal or from `data_list` - added when
# infer_size_from_fit()/the `size_resolver` argument were introduced.

mock_fit_matrix <- function(draws) {
  list(draws = function(variables, format = "matrix") {
    cols <- unlist(lapply(variables, function(v) {
      grep(paste0("^", v, "(\\[|$)"), colnames(draws), value = TRUE)
    }))
    draws[, cols, drop = FALSE]
  })
}

test_that("infer_size_from_fit() reads a vector's length from the fit's columns", {
  draws <- matrix(rnorm(300), ncol = 3,
    dimnames = list(NULL, paste0("beta[", 1:3, "]")))
  fit <- mock_fit_matrix(draws)

  expect_equal(infer_size_from_fit(fit, "beta", "vector"), 3)
})

test_that("infer_size_from_fit() reads a matrix's dimensions from the fit's columns", {
  draws <- matrix(rnorm(400), ncol = 4,
    dimnames = list(NULL, c("Sigma[1,1]", "Sigma[2,1]", "Sigma[1,2]", "Sigma[2,2]")))
  fit <- mock_fit_matrix(draws)

  expect_equal(infer_size_from_fit(fit, "Sigma", "matrix"), c(2L, 2L))
})

test_that("infer_size_from_fit() returns NA when the parameter isn't in the fit", {
  draws <- matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "mu"))
  fit <- mock_fit_matrix(draws)

  expect_true(is.na(infer_size_from_fit(fit, "beta", "vector")))
  expect_true(all(is.na(infer_size_from_fit(fit, "Sigma", "matrix"))))
})

test_that("infer_size_from_fit() returns NA (not an error) when extraction itself fails", {
  fit <- list(draws = function(variables, format = "matrix") stop("no such variable"))

  expect_true(is.na(infer_size_from_fit(fit, "beta", "vector")))
})

test_that("Create_df_priors() infers a vector size from stan_fit when data_list can't resolve it", {
  stan_code <- "
    data {
      int<lower=1> K;
    }
    parameters {
      vector[K] beta;
    }
    model {
      beta ~ normal(0, 1);
    }
  "
  draws <- matrix(rnorm(3000), ncol = 3,
    dimnames = list(NULL, paste0("beta[", 1:3, "]")))
  fit <- mock_fit_matrix(draws)

  df <- Create_df_priors(stan_code, stan_fit = fit)
  df <- df[order(df$par), ]

  expect_equal(df$par, c("beta[1]", "beta[2]", "beta[3]"))
})

test_that("Create_df_priors() infers a matrix's dimensions from stan_fit when data_list can't resolve them", {
  stan_code <- "
    data {
      int<lower=1> R;
      int<lower=1> C;
    }
    parameters {
      matrix[R, C] Sigma;
    }
    model {
      to_vector(Sigma) ~ normal(0, 1);
    }
  "
  draws <- matrix(rnorm(400), ncol = 4,
    dimnames = list(NULL, c("Sigma[1,1]", "Sigma[2,1]", "Sigma[1,2]", "Sigma[2,2]")))
  fit <- mock_fit_matrix(draws)

  df <- Create_df_priors(stan_code, stan_fit = fit)
  df <- df[order(df$par), ]

  expect_equal(df$par, c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,1]", "Sigma[2,2]"))
})

test_that("data_list still takes priority over stan_fit when both resolve a size", {
  stan_code <- "
    data {
      int<lower=1> K;
    }
    parameters {
      vector[K] beta;
    }
    model {
      beta ~ normal(0, 1);
    }
  "
  # fit only exposes 2 elements, but data_list says K = 3 - data_list wins
  draws <- matrix(rnorm(2000), ncol = 2,
    dimnames = list(NULL, paste0("beta[", 1:2, "]")))
  fit <- mock_fit_matrix(draws)

  df <- Create_df_priors(stan_code, data_list = list(K = 3), stan_fit = fit)

  expect_equal(nrow(df), 3)
})

test_that("an unresolved size still warns and skips when stan_fit doesn't have it either", {
  stan_code <- "
    data {
      int<lower=1> K;
    }
    parameters {
      vector[K] beta;
      real mu;
    }
    model {
      beta ~ normal(0, 1);
      mu ~ normal(0, 1);
    }
  "
  draws <- matrix(rnorm(1000), ncol = 1, dimnames = list(NULL, "mu"))
  fit <- mock_fit_matrix(draws)

  expect_warning(
    df <- Create_df_priors(stan_code, stan_fit = fit),
    "posterior draws so its size can be inferred"
  )
  expect_equal(df$par, "mu")
})

test_that("omitting stan_fit preserves the original warn-and-skip behaviour", {
  stan_code <- "
    data {
      int<lower=1> K;
    }
    parameters {
      vector[K] beta;
      real mu;
    }
    model {
      beta ~ normal(0, 1);
      mu ~ normal(0, 1);
    }
  "
  expect_warning(
    df <- Create_df_priors(stan_code),
    "posterior draws so its size can be inferred"
  )
  expect_equal(df$par, "mu")
})

test_that("PriorPosteriorPlot() succeeds without data_list when stan_fit resolves the size", {
  stan_code <- "
    data {
      int<lower=1> K;
    }
    parameters {
      vector[K] beta;
    }
    model {
      beta ~ normal(0, 1);
    }
  "
  draws <- matrix(rnorm(3000, 0, 1), ncol = 3,
    dimnames = list(NULL, paste0("beta[", 1:3, "]")))
  fit <- mock_fit_matrix(draws)

  p <- PriorPosteriorPlot(fit, stan_code)
  built <- ggplot2::ggplot_build(p)

  expect_equal(length(unique(built$layout$layout$PANEL)), 3)
})
