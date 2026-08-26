# Tests for lkj_corr_cholesky() support added 2026-08-26: a
# cholesky_factor_corr[K] parameter's raw elements aren't individually
# meaningful, so instead of a row for the Cholesky factor itself,
# Create_df_priors() emits one row per pairwise correlation - matched to a
# companion correlation-matrix generated quantity, auto-detected via the
# standard `multiply_lower_tri_self_transpose(L)` idiom.

test_that("dlkjcorr1/plkjcorr1/qlkjcorr1/rlkjcorr1 implement the shifted-Beta marginal correctly", {
  eta <- 3; K <- 4
  alpha <- eta + (K - 2) / 2

  # by construction: (r+1)/2 ~ Beta(alpha, alpha)
  r <- seq(-0.9, 0.9, by = 0.1)
  expect_equal(dlkjcorr1(r, eta, K), stats::dbeta((r + 1) / 2, alpha, alpha) / 2)
  expect_equal(plkjcorr1(r, eta, K), stats::pbeta((r + 1) / 2, alpha, alpha))
  expect_equal(qlkjcorr1(plkjcorr1(r, eta, K), eta, K), r)

  # symmetric around 0, and integrates to 1 over (-1, 1)
  expect_equal(plkjcorr1(0, eta, K), 0.5)
  expect_equal(stats::integrate(dlkjcorr1, -1, 1, eta = eta, K = K)$value, 1, tolerance = 1e-6)

  set.seed(1)
  draws <- rlkjcorr1(20000, eta, K)
  expect_true(all(draws > -1 & draws < 1))
  expect_equal(mean(draws), 0, tolerance = 0.02)
})

test_that("K=2 (the common two-correlated-random-effects case) produces exactly one row", {
  stan_code <- "
    parameters {
      cholesky_factor_corr[2] L_block;
    }
    transformed parameters {
      matrix[2, 2] junk = L_block * L_block';
    }
    generated quantities {
      corr_matrix[2] Corr_block = multiply_lower_tri_self_transpose(L_block);
    }
    model {
      L_block ~ lkj_corr_cholesky(4);
    }
  "
  df <- Create_df_priors(stan_code)
  df <- df[!is.na(df$dist) & df$dist == "lkj_corr", ]

  expect_equal(nrow(df), 1)
  expect_equal(df$par, "Corr_block[1,2]")
  expect_equal(df$v_min, -1); expect_equal(df$v_max, 1)
  expect_equal(df$v_lwr, -1); expect_equal(df$v_upr, 1)
  expect_equal(df$arg1, 4)  # eta
  expect_equal(df$arg2, 2)  # K
})

test_that("K=3 produces K*(K-1)/2 = 3 rows, one per pair, all sharing the same eta/K", {
  stan_code <- "
    parameters {
      cholesky_factor_corr[3] L;
    }
    transformed parameters {
      corr_matrix[3] Corr_L = multiply_lower_tri_self_transpose(L);
    }
    model {
      L ~ lkj_corr_cholesky(2);
    }
  "
  df <- Create_df_priors(stan_code)
  df <- df[order(df$par), ]

  expect_equal(df$par, c("Corr_L[1,2]", "Corr_L[1,3]", "Corr_L[2,3]"))
  expect_true(all(df$dist == "lkj_corr"))
  expect_true(all(df$arg1 == 2))
  expect_true(all(df$arg2 == 3))
})

test_that("no matching correlation-matrix generated quantity: warns with a copy-pasteable suggested line, and adds no rows", {
  stan_code <- "
    parameters {
      real<lower=0.01> epsilon;
      cholesky_factor_corr[2] L_block;
    }
    model {
      epsilon ~ normal(0, 0.15);
      L_block ~ lkj_corr_cholesky(4);
    }
  "
  expect_warning(
    df <- Create_df_priors(stan_code),
    "corr_matrix\\[2\\] Corr_block = multiply_lower_tri_self_transpose\\(L_block\\)"
  )
  expect_false("lkj_corr" %in% df$dist)
  expect_equal(df$par, "epsilon")
})

test_that("an L not requested via `pars` produces no warning even when its correlation matrix is missing", {
  stan_code <- "
    parameters {
      real<lower=0.01> epsilon;
      cholesky_factor_corr[2] L_block;
    }
    model {
      epsilon ~ normal(0, 0.15);
      L_block ~ lkj_corr_cholesky(4);
    }
  "
  df <- expect_no_warning(Create_df_priors(stan_code, pars = "epsilon"))
  expect_equal(df$par, "epsilon")
})

test_that("a correlation matrix declared inside 'transformed parameters' doesn't also trigger a redundant 'unsupported type' warning", {
  stan_code <- "
    parameters {
      cholesky_factor_corr[2] L_block;
    }
    transformed parameters {
      corr_matrix[2] Corr_block = multiply_lower_tri_self_transpose(L_block);
    }
    model {
      L_block ~ lkj_corr_cholesky(4);
    }
  "
  df <- expect_no_warning(Create_df_priors(stan_code))
  expect_equal(df$par, "Corr_block[1,2]")
})

test_that("PriorPosteriorPlotStan() renders an lkj_corr row using the true [-1, 1] hard bound, no widening warning", {
  n_draws <- 500
  set.seed(1)
  mock_draws <- matrix(rlkjcorr1(n_draws, eta = 4, K = 2), nrow = n_draws,
    dimnames = list(NULL, "Corr_block[1,2]"))
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws)

  df_priors <- data.frame(
    par = "Corr_block[1,2]", v_min = -1, v_max = 1, v_lwr = -1, v_upr = 1,
    dist = "lkj_corr", arg1 = 4, arg2 = 2
  )

  p <- expect_no_warning(PriorPosteriorPlotStan(mock_fit, "Corr_block[1,2]", df_priors))
  built <- ggplot2::ggplot_build(p)
  expect_equal(range(built$data[[1]]$x), c(-1, 1))
})

test_that("an lkj_corr_cholesky() sampling statement using an unsupported alternative distribution warns clearly", {
  stan_code <- "
    parameters {
      real<lower=0.01> epsilon;
      cholesky_factor_corr[2] L_block;
    }
    transformed parameters {
      corr_matrix[2] Corr_block = multiply_lower_tri_self_transpose(L_block);
    }
    model {
      epsilon ~ normal(0, 0.15);
      L_block ~ normal(0, 1);
    }
  "
  expect_warning(
    df <- Create_df_priors(stan_code),
    "Unsupported prior distribution 'normal' for cholesky_factor_corr parameter 'L_block'"
  )
  expect_false("lkj_corr" %in% df$dist)
})
