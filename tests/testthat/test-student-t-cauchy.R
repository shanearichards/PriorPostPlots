# Tests for student_t and cauchy prior support added 2026-08-26.
#
# student_t(nu, mu, sigma) is Stan's three-argument form; base R's dt()/pt()/
# qt() only implement the standard (location 0, scale 1) case, so dstudent_t()
# etc. add the location-scale transform. This is also the first
# three-argument family the package has supported, so it's the first time
# df_priors carries a genuinely-used `arg3` column - these tests also cover
# that the addition is backward compatible with existing two-argument-only
# df_priors data frames that predate arg3's existence.
#
# cauchy(mu, sigma) is a plain two-argument family (no schema changes
# needed), included here because a <lower=0> bound on it renders as
# half-cauchy via the same generic bound-truncation machinery normal +
# <lower=0> already uses for half-normal.

test_that("dstudent_t/pstudent_t/qstudent_t/rstudent_t implement Stan's location-scale student_t(nu, mu, sigma)", {
  nu <- 4; mu <- 2; sigma <- 3
  x <- seq(-8, 12, by = 0.5)

  expect_equal(dstudent_t(x, nu, mu, sigma), stats::dt((x - mu) / sigma, nu) / sigma)
  expect_equal(pstudent_t(x, nu, mu, sigma), stats::pt((x - mu) / sigma, nu))
  expect_equal(qstudent_t(pstudent_t(x, nu, mu, sigma), nu, mu, sigma), x)

  # default location/scale (0, 1) reduces to the standard student_t
  expect_equal(dstudent_t(x, nu), stats::dt(x, nu))

  set.seed(1)
  draws <- rstudent_t(20000, nu, mu, sigma)
  expect_equal(mean(draws), mu, tolerance = 0.1)
})

test_that("Create_df_priors() parses a student_t() sampling statement into a 3-argument row", {
  stan_code <- "
    parameters {
      real beta;
    }
    model {
      beta ~ student_t(4, 0, 2.5);
    }
  "
  df <- Create_df_priors(stan_code)
  expect_equal(df$dist, "student_t")
  expect_equal(df$arg1, 4)    # nu
  expect_equal(df$arg2, 0)    # mu
  expect_equal(df$arg3, 2.5)  # sigma
  # display range falls back to the 2.5%/97.5% quantile of student_t(4, 0, 2.5)
  expect_equal(df$v_min, qstudent_t(0.025, 4, 0, 2.5))
  expect_equal(df$v_max, qstudent_t(0.975, 4, 0, 2.5))
})

test_that("Create_df_priors() parses a cauchy() sampling statement, and <lower=0> renders half-cauchy", {
  stan_code <- "
    parameters {
      real<lower=0> sigma;
    }
    model {
      sigma ~ cauchy(0, 1);
    }
  "
  df <- Create_df_priors(stan_code)
  expect_equal(df$dist, "cauchy")
  expect_equal(df$arg1, 0)
  expect_equal(df$arg2, 1)
  expect_true(is.na(df$arg3))
  expect_equal(df$v_lwr, 0)  # the true hard bound - half-cauchy, not full cauchy
})

test_that("PriorPosteriorPlotStan() renders student_t and cauchy rows without error, ribbon covering the histogram", {
  n_draws <- 500
  set.seed(1)
  mock_draws <- cbind(
    beta  = rstudent_t(n_draws, 4, 0, 2.5) * 0.3,
    sigma = abs(stats::rcauchy(n_draws, 0, 1)) * 0.3 + 0.05
  )
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws[, variables, drop = FALSE])

  df_priors <- data.frame(
    par = c("beta", "sigma"),
    v_min = c(NA, 0), v_max = c(NA, NA), v_lwr = c(NA, 0), v_upr = c(NA, NA),
    dist = c("student_t", "cauchy"), arg1 = c(4, 0), arg2 = c(0, 1), arg3 = c(2.5, NA)
  )
  df_priors$v_min[1] <- qstudent_t(0.025, 4, 0, 2.5)
  df_priors$v_max[1] <- qstudent_t(0.975, 4, 0, 2.5)
  df_priors$v_max[2] <- stats::qcauchy(0.975, 0, 1)

  p <- expect_no_warning(PriorPosteriorPlotStan(mock_fit, c("beta", "sigma"), df_priors, nbins = 15))
  built <- ggplot2::ggplot_build(p)
  layout <- built$layout$layout

  hist_range <- dplyr::left_join(built$data[[2]], layout, by = "PANEL") |>
    dplyr::summarise(hmin = min(.data$xmin), hmax = max(.data$xmax), .by = "par")
  ribbon_range <- dplyr::left_join(built$data[[1]], layout, by = "PANEL") |>
    dplyr::summarise(rmin = min(.data$x), rmax = max(.data$x), .by = "par")
  check <- dplyr::left_join(hist_range, ribbon_range, by = "par")

  expect_true(all(check$hmin >= check$rmin - 1e-9))
  expect_true(all(check$hmax <= check$rmax + 1e-9))
})

test_that("a df_priors data frame with no arg3 column (pre-existing two-argument-only usage) still works", {
  df_old <- data.frame(
    par = "a", v_min = -3, v_max = 3, v_lwr = NA_real_, v_upr = NA_real_,
    dist = "normal", arg1 = 0, arg2 = 1
  )
  expect_false("arg3" %in% colnames(df_old))

  mock_draws <- matrix(rnorm(200), ncol = 1, dimnames = list(NULL, "a"))
  mock_fit <- list(draws = function(variables, format = "matrix") mock_draws)

  expect_no_error(PriorPosteriorPlotStan(mock_fit, "a", df_old))
})
