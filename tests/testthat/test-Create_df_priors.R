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

# Regression test for a real usability issue found 2026-08-26: calling
# Create_df_priors() (directly, or via PriorPosteriorPlot()) on a real
# hierarchical model with a `pars` subset still warned about every OTHER
# parameter in the model - unsupported types (e.g. cholesky_factor_corr),
# unsupported prior families (e.g. std_normal) on random-effect z-scores -
# none of which the caller ever asked to plot. The `pars` argument now lets
# Create_df_priors() skip anything not requested, silently.
stan_code_mixed <- "
  parameters {
    real<lower=0.01> epsilon;
    cholesky_factor_corr[2] L_block;
    vector[2] z_block;
  }
  model {
    epsilon ~ normal(0, 0.15);
    z_block ~ std_normal();
  }
"

test_that("without a `pars` filter, every parameter is parsed and warned about as before (no regression)", {
  expect_warning(
    expect_warning(
      df <- Create_df_priors(stan_code_mixed),
      "unsupported type"
    ),
    "Unsupported prior distribution"
  )
  expect_true(all(c("epsilon", "z_block[1]", "z_block[2]") %in% df$par))
})

test_that("with a `pars` filter, unrequested parameters generate no warning and are excluded", {
  df <- expect_no_warning(Create_df_priors(stan_code_mixed, pars = "epsilon"))
  expect_equal(df$par, "epsilon")
})

test_that("a `pars` filter naming one specific vector element still parses that whole vector's declaration, and still warns about that vector's own (genuinely unsupported) prior", {
  # z_block's declaration is processed (its bounds/size are needed to
  # produce the z_block[1] row actually asked for), so its std_normal
  # prior - actually attached to the parameter that WAS requested, unlike
  # cholesky_factor_corr[2] L_block or the other filtered-out declarations
  # above - is still correctly flagged as unsupported, not silenced.
  expect_warning(
    df <- Create_df_priors(stan_code_mixed, pars = "z_block[1]"),
    "Unsupported prior distribution"
  )
  expect_equal(sort(df$par), c("z_block[1]", "z_block[2]"))
})

test_that("PriorPosteriorPlot() forwards `pars` into Create_df_priors() so unrequested parameters stay silent", {
  n_draws <- 50
  # kept comfortably inside epsilon's real=<lower=0.01> bound and its
  # normal(0, 0.15) prior's own display window, so the only thing this
  # test can catch is stray warnings about z_block/L_block leaking through
  mock_draws <- matrix(seq(0.05, 0.2, length.out = n_draws), nrow = n_draws,
    dimnames = list(NULL, "epsilon"))
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws)

  expect_no_warning(PriorPosteriorPlot(mock_fit, stan_code_mixed, pars = "epsilon"))
})
