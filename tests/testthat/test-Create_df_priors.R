# Regression tests for the generality fixes made 2026-08-26:
#   1. Per-element indexed sampling statements (`par[1] ~ dist(...);`,
#      `par[1,2] ~ dist(...);`) are now detected, for both vector and
#      matrix parameters - previously only a bare/vectorized statement on
#      the whole parameter was recognised.
#   2. Sampling-statement arguments containing nested parentheses (e.g. a
#      function call like `log(0.5)`, or a compound expression like
#      `log(0.5) + 1`) are now extracted and resolved correctly -
#      previously the argument-capturing regex stopped at the FIRST
#      closing paren, silently truncating/losing later arguments.
# These were found while applying the package to real Stan models with
# exactly these patterns (see the package's GitHub issue tracker/commit
# history for the originating report).

test_that("per-element indexed vector priors are detected", {
  stan_code <- "
    parameters {
      vector[2] a;
    }
    model {
      a[1] ~ normal(0.41, 1);
      a[2] ~ normal(-0.94, 1);
    }
  "
  df <- suppressWarnings(Create_df_priors(stan_code))
  df <- df[order(df$par), ]

  expect_equal(df$par, c("a[1]", "a[2]"))
  expect_equal(df$dist, c("normal", "normal"))
  expect_equal(df$arg1, c(0.41, -0.94))
  expect_equal(df$arg2, c(1, 1))
})

test_that("per-element indexed matrix priors are detected, including a mix of literal and expression arguments", {
  stan_code <- "
    parameters {
      matrix[2, 2] Sigma;
    }
    model {
      Sigma[1,1] ~ normal(0, 1);
      Sigma[2,1] ~ normal(0, sqrt(2));
      Sigma[1,2] ~ normal(0, 1);
      Sigma[2,2] ~ normal(0, 1);
    }
  "
  df <- suppressWarnings(Create_df_priors(stan_code))
  df <- df[order(df$par), ]

  expect_equal(df$par, c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,1]", "Sigma[2,2]"))
  expect_equal(df$arg2[df$par == "Sigma[2,1]"], sqrt(2))
  expect_equal(df$arg2[df$par == "Sigma[1,1]"], 1)
})

test_that("a vectorized statement still applies to every expanded element (no regression)", {
  stan_code <- "
    parameters {
      vector[2] sigma_block;
    }
    model {
      sigma_block ~ normal(0, 0.5);
    }
  "
  df <- suppressWarnings(Create_df_priors(stan_code))
  df <- df[order(df$par), ]

  expect_equal(df$par, c("sigma_block[1]", "sigma_block[2]"))
  expect_true(all(df$dist == "normal"))
  expect_true(all(df$arg1 == 0))
  expect_true(all(df$arg2 == 0.5))
})

test_that("a per-element statement takes precedence over a vectorized statement on the same base parameter", {
  stan_code <- "
    parameters {
      vector[2] a;
    }
    model {
      a ~ normal(0, 1);
      a[2] ~ normal(5, 1);
    }
  "
  df <- suppressWarnings(Create_df_priors(stan_code))
  df <- df[order(df$par), ]

  expect_equal(df$arg1[df$par == "a[1]"], 0)   # falls back to the vectorized statement
  expect_equal(df$arg1[df$par == "a[2]"], 5)   # per-element statement wins
})

test_that("a nested function call as a prior argument is resolved, not truncated", {
  stan_code <- "
    parameters {
      real log_alpha;
    }
    model {
      log_alpha ~ normal(log(0.5), 0.75);
    }
  "
  df <- suppressWarnings(Create_df_priors(stan_code))

  expect_equal(df$dist, "normal")
  expect_equal(df$arg1, log(0.5))
  expect_equal(df$arg2, 0.75)  # previously lost entirely (regex stopped at log(0.5)'s own ")")
})

test_that("a compound expression (nested call + arithmetic) as a prior argument resolves correctly", {
  stan_code <- "
    parameters {
      real c;
    }
    model {
      c ~ normal(log(0.5) + 1, 0.75);
    }
  "
  df <- suppressWarnings(Create_df_priors(stan_code))
  expect_equal(df$arg1, log(0.5) + 1)
})

test_that("resolve_numeric() still resolves plain literals and data_list constants (no regression)", {
  expect_equal(resolve_numeric("0.75"), 0.75)
  expect_equal(resolve_numeric("K", data_list = list(K = 4)), 4)
})

test_that("resolve_numeric() resolves an expression combining a data_list constant with arithmetic", {
  expect_equal(resolve_numeric("1 / K", data_list = list(K = 2)), 0.5)
})

test_that("find_matching_paren() handles nesting and reports NA when unbalanced", {
  text <- "normal(log(0.5), 0.75)"
  open_pos <- regexpr("(", text, fixed = TRUE)[1]
  expect_equal(find_matching_paren(text, open_pos), nchar(text))
  expect_true(is.na(find_matching_paren("normal(log(0.5), 0.75", 7)))
})

test_that("split_args_respecting_parens() doesn't split inside a nested call", {
  expect_equal(
    split_args_respecting_parens("log(0.5), 0.75"),
    c("log(0.5)", "0.75")
  )
  expect_equal(
    split_args_respecting_parens("some_func(a, b), 0.75"),
    c("some_func(a, b)", "0.75")
  )
})
