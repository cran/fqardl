## Regression tests for the 1.0.3 correctness release.
## These exist so that the defects fixed in 1.0.3 cannot silently return.

test_that("the design's left-hand side is Delta y, not y", {
  data(macro_data, package = "fqardl")
  y <- macro_data$gdp
  X <- as.matrix(macro_data[, c("inflation", "interest_rate")])
  fr <- generate_fourier_terms(length(y), 1)
  d <- fqardl:::build_ardl_design(y, X, fr, p = 1, q = 1)
  n <- length(y); m <- d$max_lag
  expect_equal(d$y, y[(m + 1):n] - y[m:(n - 1)])
  expect_false(isTRUE(all.equal(d$y, y[(m + 1):n])))
})

test_that("the error-correction coefficient is negative on a cointegrated DGP", {
  set.seed(101); n <- 250
  x <- cumsum(rnorm(n)); y <- numeric(n); y[1] <- 0.8 * x[1]
  for (t in 2:n)
    y[t] <- y[t-1] - 0.35*(y[t-1] - 2 - 0.8*x[t-1]) + 0.8*(x[t]-x[t-1]) + rnorm(1, 0, .5)
  set.seed(7)
  f <- fqardl(y ~ x1, data = data.frame(y = y, x1 = x), tau = c(.25, .5, .75),
              max_k = 1, max_p = 2, max_q = 2, verbose = FALSE)
  ect <- vapply(f$qardl_results, function(r) unname(r$coefficients["y_lag1"]), numeric(1))
  expect_true(all(ect < 0))
  expect_true(all(f$bounds_test$summary$t_stat < 0))
})

test_that("the estimator recovers a known long-run coefficient", {
  # A consistency check, so it is run over replications rather than one draw:
  # a single seed is a fragile test of an asymptotic property.
  gen <- function(seed, n = 300) {
    set.seed(seed)
    x <- cumsum(rnorm(n)); y <- numeric(n); y[1] <- 0.8 * x[1]
    for (t in 2:n)
      y[t] <- y[t-1] - 0.35*(y[t-1] - 2 - 0.8*x[t-1]) + 0.8*(x[t]-x[t-1]) + rnorm(1, 0, .5)
    data.frame(y = y, x1 = x)
  }
  lr <- ect <- rep(NA_real_, 15)
  for (i in 1:15) {
    set.seed(7)
    f <- try(fqardl(y ~ x1, data = gen(300 + i), tau = 0.5,
                    max_k = 1, max_p = 2, max_q = 2, verbose = FALSE), silent = TRUE)
    if (inherits(f, "try-error")) next
    lr[i]  <- unname(f$long_run[1, "x1"])
    ect[i] <- unname(f$qardl_results[[1]]$coefficients["y_lag1"])
  }
  expect_equal(median(lr,  na.rm = TRUE),  0.80, tolerance = 0.08)
  expect_equal(median(ect, na.rm = TRUE), -0.35, tolerance = 0.10)
  expect_true(all(lr > 0, na.rm = TRUE))   # sign was wrong in 100% of cases in 1.0.2
})

test_that("the bounds F is a Wald statistic, not mean(t^2)", {
  data(macro_data, package = "fqardl")
  set.seed(123)
  f <- fqardl(gdp ~ inflation + interest_rate, data = macro_data,
              tau = 0.5, max_k = 1, verbose = FALSE)
  r <- f$qardl_results[[1]]
  lev <- c("y_lag1", "x1_lag1", "x2_lag1")
  pos <- match(lev, names(r$coefficients))
  R <- matrix(0, length(pos), length(r$coefficients))
  R[cbind(seq_along(pos), pos)] <- 1
  Rb <- R %*% r$coefficients
  Fman <- as.numeric(t(Rb) %*% solve(R %*% r$vcov %*% t(R)) %*% Rb) / length(pos)
  expect_equal(f$bounds_test$F_stat, Fman)

  # and it is NOT the old construction
  Fold <- mean(r$t_statistics[lev]^2)
  expect_false(isTRUE(all.equal(f$bounds_test$F_stat, Fold)))
})

test_that("the bounds test is reported at every quantile", {
  data(macro_data, package = "fqardl")
  taus <- c(0.25, 0.5, 0.75)
  set.seed(123)
  f <- fqardl(gdp ~ inflation, data = macro_data, tau = taus,
              max_k = 1, verbose = FALSE)
  expect_equal(nrow(f$bounds_test$summary), length(taus))
  expect_equal(f$bounds_test$summary$tau, taus)
})

test_that("the inconclusive region is reachable", {
  verdicts <- vapply(c(3.0, 5.0, 4.2), function(s)
    fqardl:::bounds_verdict(s, 3.79, 4.85), character(1))
  expect_match(verdicts[1], "No cointegration")
  expect_match(verdicts[2], "Cointegration")
  expect_match(verdicts[3], "Inconclusive")
})

test_that("case other than 3 raises an error rather than reporting a wrong table", {
  data(macro_data, package = "fqardl")
  expect_error(
    fqardl(gdp ~ inflation, data = macro_data, tau = 0.5, max_k = 1,
           case = 5, verbose = FALSE),
    "Only case = 3")
  expect_error(fqardl:::get_pss_critical_values(2, case = 4), "Only case = 3")
})

test_that("coefficients, standard errors and t statistics are named consistently", {
  data(macro_data, package = "fqardl")
  fr <- generate_fourier_terms(nrow(macro_data), 1)
  r <- estimate_qardl(macro_data$gdp,
                      as.matrix(macro_data[, c("inflation", "interest_rate")]),
                      fr, p = 1, q = 1, tau = 0.5)
  expect_identical(names(r$coefficients), names(r$std_errors))
  expect_identical(names(r$coefficients), names(r$t_statistics))
  expect_false(is.na(r$std_errors["y_lag1"]))
})

test_that("the shipped critical values match PSS Tables CI(iii) and CII(iii) at 5 percent", {
  expected <- rbind(
    c(1,  4.94, 5.73, -2.86, -3.22),
    c(2,  3.79, 4.85, -2.86, -3.53),
    c(3,  3.23, 4.35, -2.86, -3.78),
    c(4,  2.86, 4.01, -2.86, -3.99),
    c(5,  2.62, 3.79, -2.86, -4.19),
    c(6,  2.45, 3.61, -2.86, -4.38),
    c(7,  2.32, 3.50, -2.86, -4.57),
    c(8,  2.22, 3.39, -2.86, -4.72),
    c(9,  2.14, 3.30, -2.86, -4.88),
    c(10, 2.06, 3.24, -2.86, -5.03))
  for (i in seq_len(nrow(expected))) {
    cv <- fqardl:::get_pss_critical_values(expected[i, 1], case = 3)
    expect_equal(c(cv$F_lower[2], cv$F_upper[2], cv$t_lower[2], cv$t_upper[2]),
                 as.numeric(expected[i, 2:5]),
                 info = paste("k =", expected[i, 1]))
  }
})

test_that("mtnardl returns a t statistic and a negative ECT", {
  data(oil_gdp_data, package = "fqardl")
  set.seed(11)
  mt <- mtnardl(gdp ~ oil_price, data = oil_gdp_data, verbose = FALSE)
  expect_true(is.finite(mt$bounds_test$t_stat))
  expect_lt(unname(mt$model$coefficients["y_lag1"]), 0)
})

test_that("routines with known limitations warn at the point of use", {
  data(macro_data, package = "fqardl")
  set.seed(1)
  f <- fqardl(gdp ~ inflation, data = macro_data, tau = c(0.25, 0.75),
              max_k = 1, verbose = FALSE)
  expect_warning(quantile_wald_test(f$qardl_results, "y_lag1"),
                 "independence across quantiles")
})

test_that("only the verified 5 percent critical values are supplied", {
  # An adversarial audit of 1.0.3 found that the 1 and 10 percent F bounds had
  # no traceable provenance. They are now NA rather than plausible-looking
  # numbers of unknown origin, and the decision rule uses 5 percent only.
  for (k in 1:10) {
    cv <- fqardl:::get_pss_critical_values(k, case = 3)
    expect_true(is.na(cv$F_lower[1]) && is.na(cv$F_lower[3]))
    expect_true(is.na(cv$F_upper[1]) && is.na(cv$F_upper[3]))
    expect_true(is.na(cv$t_upper[1]) && is.na(cv$t_upper[3]))
    expect_false(is.na(cv$F_lower[2]))
    expect_false(is.na(cv$F_upper[2]))
    expect_identical(cv$verified, c(FALSE, TRUE, FALSE))
  }
})

test_that("the shipped table differs from the wrong one shipped in 1.0.2", {
  # 1.0.2's Case III table was correct only at k = 1; from k = 2 it diverged
  # from PSS Table CI(iii), and for k > 6 it returned one row for every k.
  wrong_102 <- list("2" = c(3.55, 4.38), "3" = c(3.10, 3.99), "5" = c(2.69, 3.61))
  for (kk in names(wrong_102)) {
    cv <- fqardl:::get_pss_critical_values(as.integer(kk), case = 3)
    expect_false(isTRUE(all.equal(c(cv$F_lower[2], cv$F_upper[2]),
                                 wrong_102[[kk]], tolerance = 1e-8)))
  }
  # and it must still vary with k beyond 6, which 1.0.2 did not
  a <- fqardl:::get_pss_critical_values(7,  case = 3)$F_upper[2]
  b <- fqardl:::get_pss_critical_values(10, case = 3)$F_upper[2]
  expect_true(a > b)
})
