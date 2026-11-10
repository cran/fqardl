## Checks of test_asymmetry() and dynamic_multipliers() against
## Shin, Yu and Greenwood-Nimmo (2014), sections 9.2.2 to 9.2.4.

nardl_fit <- function(seed, n = 150, p = 3, q = 3) {
  set.seed(seed)
  x <- cumsum(rnorm(n))
  d <- decompose_variables(data.frame(x = x), "x")
  xp <- d$x$positive; xn <- d$x$negative
  y <- numeric(n)
  for (t in 2:n)
    y[t] <- y[t - 1] - 0.3 * (y[t - 1] - 1.0 * xp[t - 1] - 0.4 * xn[t - 1]) +
      0.5 * (xp[t] - xp[t - 1]) + 0.1 * (xn[t] - xn[t - 1]) + rnorm(1, sd = 0.5)
  X <- cbind(x_pos = xp, x_neg = xn)
  estimate_nardl(y, X, generate_fourier_terms(n, 1), p, q, 3)
}

test_that("long-run symmetry: delta-method Wald on L(+) = L(-)", {
  r <- nardl_fit(1)
  b <- r$coefficients; V <- vcov(r$model)
  i_r <- which(names(b) == "y_lag1")
  i_p <- which(names(b) == "x_pos_lag1"); i_n <- which(names(b) == "x_neg_lag1")
  Lp <- -b[i_p] / b[i_r]; Ln <- -b[i_n] / b[i_r]
  g <- numeric(length(b))
  g[i_p] <- -1 / b[i_r]; g[i_n] <- 1 / b[i_r]
  g[i_r] <- (b[i_p] - b[i_n]) / b[i_r]^2
  W <- unname((Lp - Ln)^2 / drop(t(g) %*% V %*% g))
  a <- test_asymmetry(r, "x")$x
  expect_equal(a$wald_stat, W, tolerance = 1e-10)
  expect_equal(a$p_value, pchisq(W, 1, lower.tail = FALSE), tolerance = 1e-10)
  expect_equal(a$lr_positive, unname(Lp)); expect_equal(a$lr_negative, unname(Ln))
  expect_equal(a$difference, unname(Lp - Ln))
  expect_equal(a$coef_positive, unname(b[i_p]))
})

test_that("short-run symmetry: Wald equals q times the F test of the restricted OLS fit", {
  r <- nardl_fit(2, p = 2, q = 3)
  a <- test_asymmetry(r, "x")$x
  Xd <- model.matrix(r$model); colnames(Xd) <- names(r$coefficients)
  yd <- model.response(model.frame(r$model))
  full <- lm(yd ~ Xd - 1)
  dp <- c("d_x_pos", "d_x_pos_lag1", "d_x_pos_lag2")
  dn <- c("d_x_neg", "d_x_neg_lag1", "d_x_neg_lag2")
  # pairwise: pi(+)_j = pi(-)_j, j = 0, 1, 2
  Xr <- cbind(Xd[, setdiff(colnames(Xd), c(dp, dn))], Xd[, dp] + Xd[, dn])
  Fp <- anova(lm(yd ~ Xr - 1), full)$F[2]
  expect_equal(a$short_run_pairwise$df, 3L)
  expect_equal(a$short_run_pairwise$wald_stat, 3 * Fp, tolerance = 1e-8)
  # additive: sum pi(+) = sum pi(-); reparameterise pi(-)_0 = sum pi(+) - pi(-)_1 - pi(-)_2
  Xa <- Xd[, setdiff(colnames(Xd), "d_x_neg")]
  Xa[, dp] <- Xa[, dp] + Xd[, "d_x_neg"]
  Xa[, c("d_x_neg_lag1", "d_x_neg_lag2")] <- Xa[, c("d_x_neg_lag1", "d_x_neg_lag2")] - Xd[, "d_x_neg"]
  Fa <- anova(lm(yd ~ Xa - 1), full)$F[2]
  expect_equal(a$short_run_additive$df, 1L)
  expect_equal(a$short_run_additive$wald_stat, Fa, tolerance = 1e-8)
  b <- r$coefficients
  expect_equal(a$short_run_additive$sum_positive, sum(b[dp]))
  expect_equal(a$short_run_additive$sum_negative, sum(b[dn]))
})

test_that("short-run symmetry: the two forms coincide when q = 1", {
  r <- nardl_fit(3, p = 2, q = 1)
  a <- test_asymmetry(r, "x")$x
  expect_equal(a$short_run_pairwise$df, 1L)
  expect_equal(a$short_run_pairwise$wald_stat, a$short_run_additive$wald_stat)
})

test_that("dynamic multipliers follow the levels recursion of footnotes 11 and 12", {
  r <- nardl_fit(4, p = 3, q = 3)
  b <- r$coefficients
  rho <- b[["y_lag1"]]; gam <- unname(b[c("dy_lag1", "dy_lag2")])
  H <- 40
  lev <- function(theta, pi) {
    # phi(L) of the ARDL in levels: phi_1 = 1 + rho + gamma_1,
    # phi_i = gamma_i - gamma_(i-1), phi_p = -gamma_(p-1)
    phi <- c(1 + rho + gam[1], gam[2] - gam[1], -gam[2])
    # theta(L): psi_0 = pi_0, psi_1 = theta - pi_0 + pi_1,
    # psi_i = pi_i - pi_(i-1), psi_q = -pi_(q-1)
    psi <- c(pi[1], theta - pi[1] + pi[2], pi[3] - pi[2], -pi[3])
    lam <- numeric(H + 1)
    for (j in 0:H) {
      s <- if (j + 1 <= length(psi)) psi[j + 1] else 0
      for (i in seq_along(phi)) if (j - i >= 0) s <- s + phi[i] * lam[j - i + 1]
      lam[j + 1] <- s
    }
    cumsum(lam)
  }
  mp <- lev(b[["x_pos_lag1"]], unname(b[c("d_x_pos", "d_x_pos_lag1", "d_x_pos_lag2")]))
  mn <- lev(b[["x_neg_lag1"]], unname(b[c("d_x_neg", "d_x_neg_lag1", "d_x_neg_lag2")]))
  dm <- dynamic_multipliers(r, H, "x")$x
  expect_equal(dm$m_positive, mp, tolerance = 1e-10)
  expect_equal(dm$m_negative, mn, tolerance = 1e-10)
  expect_equal(dm$asymmetry, mp - mn, tolerance = 1e-10)
  expect_equal(dm$m_positive[1], b[["d_x_pos"]])
})

test_that("dynamic multipliers tend to the long-run multipliers", {
  r <- nardl_fit(5, p = 2, q = 2)
  dm <- dynamic_multipliers(r, 300, "x")$x
  expect_equal(dm$m_positive[301], attr(dm, "lr_positive"), tolerance = 1e-6)
  expect_equal(dm$m_negative[301], attr(dm, "lr_negative"), tolerance = 1e-6)
  expect_error(dynamic_multipliers(r, -1, "x"), "horizon")
})

test_that("fnardl object: multipliers named by variable, summary and plots", {
  skip_on_cran()
  set.seed(11)
  n <- 120
  dat <- data.frame(y = cumsum(rnorm(n)), x = cumsum(rnorm(n)))
  invisible(capture.output(f <- fnardl(y ~ x, data = dat, decompose = "x",
                                       max_k = 1, verbose = FALSE)))
  expect_identical(names(f$multipliers$long_run_positive), "x")
  expect_equal(unname(f$multipliers$long_run_positive),
               f$asymmetry_tests$x$lr_positive)
  out <- capture.output(summary(f))
  expect_true(any(grepl("Short-run symmetry \\(additive\\)", out)))
  expect_false(any(grepl("NA", grep("^x\\(", out, value = TRUE))))
  expect_s3_class(plot(f, type = "dynamic"), "ggplot")
  expect_s3_class(plot(f, type = "cumulative"), "ggplot")
  expect_s3_class(plot(f, type = "comparison"), "ggplot")
  expect_s3_class(plot(f, type = "asymmetry"), "ggplot")
  dm <- dynamic_multipliers(f, 10)
  expect_identical(names(dm), "x")
  expect_equal(nrow(dm$x), 11L)
})
