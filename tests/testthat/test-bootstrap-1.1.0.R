## Tests of the 1.1.0 recursive bootstrap and of the Fourier bounds caveat.

# AND rule written out by hand: cointegration only if Fov, t and Find all
# reject; Fov without t, and Fov and t without Find, are the degenerate cases.
hand_and <- function(r) {
  if (!r[["Fov"]]) "NO_COINTEGRATION"
  else if (!r[["t"]]) "DEGENERATE_1"
  else if (!r[["Find"]]) "DEGENERATE_2"
  else "COINTEGRATION"
}

sim_rw <- function(seed, n = 80, k = 1) {
  set.seed(seed)
  d <- data.frame(y = cumsum(rnorm(n)))
  for (j in seq_len(k)) d[[paste0("x", j)]] <- cumsum(rnorm(n))
  d
}

test_that("the engine copy is version 1.1.0 and unedited", {
  expect_identical(fqardl:::.abe_engine_version, "1.1.0")
  src <- testthat::test_path("..", "..", "R", "ardl_boot_engine.R")
  skip_if_not(file.exists(src), "engine source not available")
  expect_identical(unname(tools::md5sum(src)), "bc08ba0a908a6cc92fb1ae491eced973")
})

test_that("bootstrap statistics equal restricted-RSS F and t tests by lm()", {
  d <- sim_rw(2, n = 90, k = 2)
  X <- as.matrix(d[, c("x1", "x2")])
  fr <- generate_fourier_terms(nrow(d), 2)
  st <- fqardl:::.fq_stats_qardl(d$y, X, list(fourier = fr, p = 3, q = 2))
  D <- fqardl:::build_ardl_design(d$y, X, fr, 3, 2)
  u <- lm(D$y ~ D$X - 1)
  drop_cols <- function(cols) lm(D$y ~ D$X[, !colnames(D$X) %in% cols] - 1)
  Fov <- anova(drop_cols(c("y_lag1", "x1_lag1", "x2_lag1")), u)$F[2]
  Find <- anova(drop_cols(c("x1_lag1", "x2_lag1")), u)$F[2]
  tt <- summary(u)$coefficients["D$Xy_lag1", 3]
  expect_equal(unname(st), c(Fov, tt, Find), tolerance = 1e-10)
})

test_that("bootstrap_bounds_test reproduces the data and uses the AND rule", {
  d <- sim_rw(3)
  X <- as.matrix(d["x1"])
  fr <- generate_fourier_terms(nrow(d), 1)
  for (sc in c("bvz", "mcnown")) {
    b <- bootstrap_bounds_test(d$y, X, fr, p = 2, q = 2, n_boot = 9,
                               scheme = sc, seed = 1)
    expect_lt(max(b$engine$dgpcheck), 1e-8)
    expect_lt(b$engine$design_check, 1e-8)
    expect_equal(b$engine$settings$nulls, if (sc == "bvz") "separate" else "joint")
    expect_length(b$boot_F, 9)
    # hand-coded AND rule on the observed statistics and the order-statistic
    # critical values at 5 percent
    cv <- b$engine$cv_level
    r <- c(Fov = b$orig_F > cv[["Fov"]], t = b$orig_t < cv[["t"]],
           Find = b$orig_Find > cv[["Find"]])
    expect_identical(b$decision, hand_and(r))
  }
  st <- fqardl:::.fq_stats_qardl(d$y, X, list(fourier = fr, p = 2, q = 2))
  expect_equal(c(b$orig_F, b$orig_t, b$orig_Find), unname(st))
  # the AND rule on known reject vectors: cointegration only if all reject
  known <- list(c(TRUE, TRUE, TRUE), c(FALSE, TRUE, TRUE), c(TRUE, FALSE, TRUE),
                c(TRUE, TRUE, FALSE), c(FALSE, FALSE, FALSE), c(FALSE, TRUE, FALSE))
  expected <- c("COINTEGRATION", "NO_COINTEGRATION", "DEGENERATE_1",
                "DEGENERATE_2", "NO_COINTEGRATION", "NO_COINTEGRATION")
  for (i in seq_along(known)) {
    r <- stats::setNames(known[[i]], c("Fov", "t", "Find"))
    expect_identical(fqardl:::.abe_decision(r)$decision, expected[i])
    expect_identical(hand_and(r), expected[i])
  }
})

test_that("the bootstrap leaves the caller's random numbers untouched", {
  d <- sim_rw(4, n = 60)
  fr <- generate_fourier_terms(60, 1)
  set.seed(9); a <- runif(1)
  set.seed(9)
  b1 <- bootstrap_bounds_test(d$y, as.matrix(d["x1"]), fr, 1, 1, n_boot = 5, seed = 2)
  expect_identical(runif(1), a)
  b2 <- bootstrap_bounds_test(d$y, as.matrix(d["x1"]), fr, 1, 1, n_boot = 5, seed = 2)
  expect_identical(b1$boot_F, b2$boot_F)
})

test_that("re-selection runs fqardl's own selection on every bootstrap sample", {
  d <- sim_rw(5, n = 70)
  X <- as.matrix(d["x1"])
  rs <- list(max_k = 2, max_p = 2, max_q = 2, criterion = "BIC")
  sel <- fqardl:::.fq_make_select(rs)
  s <- sel(d$y, X)
  expect_identical(s$k, select_fourier_frequency(d$y, X, 2, "BIC")$optimal_k)
  lg <- select_optimal_lags(d$y, X, generate_fourier_terms(70, s$k), 2, 2, "BIC")
  expect_identical(c(s$p, s$q), c(lg$optimal_p, lg$optimal_q))
  # a frequency fixed by the user is kept
  s2 <- fqardl:::.fq_make_select(c(rs, k = 2))(d$y, X)
  expect_identical(s2$k, 2)
  b <- bootstrap_bounds_test(d$y, X, generate_fourier_terms(70, 1), 1, 1,
                             n_boot = 5, reselect = rs, seed = 3)
  expect_true(b$reselect)
  expect_true(b$engine$settings$reselect)
})

test_that("fnardl bootstrap rebuilds the partial sums from x*", {
  d <- sim_rw(6, n = 80, k = 2)
  dec <- fqardl:::.fn_decomposer(c("x1", "x2"), "x1")
  X <- as.matrix(d[, c("x1", "x2")])
  W <- dec$full(X)
  dx <- diff(d$x1)
  expect_equal(unname(W[, 1]), cumsum(c(0, pmax(dx, 0))))
  expect_equal(unname(W[, 2]), cumsum(c(0, pmin(dx, 0))))
  expect_equal(unname(W[, 3]), d$x2)
  expect_equal(dec$step(c(-0.5, 0.3)), c(0, -0.5, 0.3))
  fr <- generate_fourier_terms(80, 1)
  b <- fqardl:::bootstrap_nardl(d$y, X, c("x1", "x2"), "x1", fr, p = 2, q = 2,
                                n_boot = 9, seed = 1)
  expect_lt(max(b$engine$dgpcheck), 1e-8)
  expect_lt(b$engine$design_check, 1e-8)
  # Find tests the three level regressors (x1+, x1-, x2)
  D <- fqardl:::build_nardl_design(d$y, W, fr, 2, 2)
  u <- lm(D$y ~ D$X - 1)
  r <- lm(D$y ~ D$X[, !colnames(D$X) %in% D$level_names] - 1)
  expect_equal(b$orig_Find, anova(r, u)$F[2], tolerance = 1e-10)
})

test_that("with Fourier terms no analytical verdict is given", {
  data(macro_data, package = "fqardl")
  set.seed(1)
  f <- fqardl(gdp ~ inflation, data = macro_data, tau = c(0.25, 0.5),
              max_k = 1, verbose = FALSE)
  expect_match(f$bounds_test$decision, "No analytical verdict")
  expect_false(f$bounds_test$bounds_valid)
  expect_match(f$bounds_test$bounds_note, "not valid with Fourier terms")
  expect_true(all(grepl("No analytical verdict", f$bounds_test$summary$decision)))
  expect_output(summary(f), "not valid with Fourier terms")
  z <- fnardl(gdp ~ oil_price, data = macro_data, max_k = 1, max_p = 2,
              max_q = 2, verbose = FALSE)
  expect_match(z$bounds_test$decision, "No analytical verdict")
  expect_false(z$bounds_test$bounds_valid)
})

test_that("fnardl bounds F uses only the lagged levels", {
  d <- sim_rw(7, n = 100)
  dd <- decompose_variables(d, "x1")
  X <- cbind(x1_pos = dd$x1$positive, x1_neg = dd$x1$negative)
  fr <- generate_fourier_terms(100, 1)
  fit <- estimate_nardl(d$y, X, fr, p = 2, q = 2)
  b <- fqardl:::perform_nardl_bounds_test(fit, 100, 2, 3)
  D <- fit$design
  u <- lm(D$y ~ D$X - 1)
  r <- lm(D$y ~ D$X[, !colnames(D$X) %in% c("y_lag1", "x1_pos_lag1", "x1_neg_lag1")] - 1)
  expect_equal(b$F_stat, anova(r, u)$F[2], tolerance = 1e-10)
  expect_equal(b$k, 2L)
})

test_that("fqardl and fnardl run the bootstrap end to end", {
  d <- sim_rw(8, n = 70)
  f <- fqardl(y ~ x1, d, tau = 0.5, max_k = 1, max_p = 1, max_q = 1,
              bootstrap = TRUE, n_boot = 5, seed = 1, verbose = FALSE)
  expect_true(f$bootstrap$decision %in% c("COINTEGRATION", "NO_COINTEGRATION",
                                          "DEGENERATE_1", "DEGENERATE_2"))
  expect_output(summary(f), "Fov AND t AND Find")
  z <- fnardl(y ~ x1, d, max_k = 1, max_p = 1, max_q = 1, bootstrap = TRUE,
              n_boot = 5, fourier_k = 1, seed = 1, verbose = FALSE)
  expect_true(z$bootstrap$k_fixed)
  expect_identical(z$bootstrap$p_value, z$bootstrap$p_value_t)
})

test_that("the fast selection in the bootstrap equals the package selection", {
  for (s in 1:12) {
    set.seed(100 + s)
    n <- c(60, 90, 120)[s %% 3 + 1]
    y <- cumsum(rnorm(n))
    X <- matrix(cumsum(rnorm(2 * n)), n)
    for (cr in c("BIC", "AIC", "HQ")) {
      k <- select_fourier_frequency(y, X, 3, cr)$optimal_k
      expect_identical(fqardl:::.fq_select_k_fast(y, X, 3, cr), k)
      fr <- generate_fourier_terms(n, k)
      a <- select_optimal_lags(y, X, fr, 3, 3, cr)
      b <- fqardl:::.fq_select_lags_fast(y, X, fr, 3, 3, cr)
      expect_identical(c(b$optimal_p, b$optimal_q), c(a$optimal_p, a$optimal_q))
    }
  }
  d <- fqardl:::build_ardl_design(y, X, fr, 2, 2)
  m <- lm(d$y ~ d$X - 1)
  expect_equal(fqardl:::.fq_loglik(qr.resid(qr(d$X), d$y)), as.numeric(logLik(m)),
               tolerance = 1e-12)
})

test_that("summary.fnardl prints the long-run multipliers, not NA", {
  data(macro_data, package = "fqardl")
  z <- fnardl(gdp ~ oil_price, data = macro_data, max_k = 1, max_p = 1,
              max_q = 1, verbose = FALSE)
  out <- capture.output(summary(z))
  lp <- grep("^oil_price\\(\\+\\)", out, value = TRUE)
  expect_false(grepl("NA", lp))
  v <- unname(z$multipliers$long_run_positive)
  expect_match(lp, sprintf("%.4f", v), fixed = TRUE)
})

test_that("bootstrap partial sums continue from the data block (engine E1)", {
  d <- sim_rw(9, n = 80)
  X <- as.matrix(d["x1"])
  fr <- generate_fourier_terms(80, 1)
  b <- fqardl:::bootstrap_nardl(d$y, X, "x1", "x1", fr, p = 2, q = 1,
                                n_boot = 9, seed = 2)
  expect_lt(max(b$engine$dgpcheck), 1e-8)
  expect_identical(b$engine$engine_version, "1.1.0")
  expect_identical(b$engine$settings$xlags, 1)
  expect_identical(b$engine$settings$j0, 0)
})
