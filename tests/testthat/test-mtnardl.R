# Tests for Multi-Threshold NARDL functions

test_that("decompose_multi_threshold refuses thresholds other than 0", {
  # CHANGED in 1.1.0: the multi-threshold allocation was wrong (the regime
  # pieces did not sum to the change), so nonzero thresholds are an error.
  set.seed(123)
  x <- cumsum(rnorm(100, sd = 2))
  expect_error(decompose_multi_threshold(x, c(-1, 0, 1)), "other than 0")
  expect_error(decompose_multi_threshold(x, 2), "other than 0")
})

test_that("decompose_multi_threshold with single threshold equals NARDL", {
  set.seed(456)
  n <- 50
  x <- cumsum(rnorm(n))
  
  result <- decompose_multi_threshold(x, 0)
  
  # Should have 2 regimes (positive and negative)
  expect_equal(ncol(result$components), 2)
  # and they are the partial sums of the positive and negative changes
  dx <- c(0, diff(x))
  expect_equal(unname(result$components[, 1]), cumsum(pmin(dx, 0)))
  expect_equal(unname(result$components[, 2]), cumsum(pmax(dx, 0)))
  expect_equal(rowSums(result$components), x - x[1])
})

test_that("mtnardl estimates model correctly", {
  skip_on_cran()
  
  set.seed(789)
  n <- 150
  
  # Generate asymmetric data
  x <- cumsum(rnorm(n, sd = 2))
  y <- numeric(n)
  y[1] <- 100
  
  for (i in 2:n) {
    dx <- x[i] - x[i-1]
    if (dx > 1) {
      effect <- -0.5 * dx  # Large positive
    } else if (dx > 0) {
      effect <- -0.2 * dx  # Small positive
    } else if (dx > -1) {
      effect <- 0.3 * abs(dx)  # Small negative
    } else {
      effect <- 0.6 * abs(dx)  # Large negative
    }
    y[i] <- 0.9 * y[i-1] + effect + rnorm(1, sd = 0.3)
  }
  
  data <- data.frame(y = y, x = x)
  
  expect_error(suppressWarnings(mtnardl(
    formula = y ~ x, data = data, decompose = "x",
    thresholds = list(x = c(-1, 0, 1)), max_p = 2, max_q = 2,
    verbose = FALSE)), "other than 0")

  result <- suppressWarnings(mtnardl(
    formula = y ~ x,
    data = data,
    decompose = "x",
    max_p = 2,
    max_q = 2,
    verbose = FALSE
  ))
  
  expect_s3_class(result, "mtnardl")
  expect_s3_class(result, "fqardl_mtnardl")
  expect_length(result$regime_names$x, 2)
  expect_true(!is.null(result$multipliers))
  expect_true(!is.null(result$regime_tests))
})

test_that("mtnardl print method works", {
  skip_on_cran()
  
  set.seed(111)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.9 * c(100, head(cumsum(rnorm(n)), -1)) + 0.3 * x + rnorm(n)
  
  data <- data.frame(y = y, x = x)
  
  result <- suppressWarnings(mtnardl(
    formula = y ~ x,
    data = data,
    decompose = "x",
    thresholds = list(x = c(0)),
    max_p = 1,
    max_q = 1,
    verbose = FALSE
  ))
  
  expect_output(print(result), "Multi-Threshold NARDL")
})

test_that("mtnardl warns once per session that it is deprecated", {
  data(oil_gdp_data, package = "fqardl")
  st <- fqardl:::.fqardl_state
  old <- st$mtnardl_warned
  on.exit(st$mtnardl_warned <- old)
  st$mtnardl_warned <- FALSE
  expect_warning(mtnardl(gdp ~ oil_price, data = oil_gdp_data, max_p = 1,
                         max_q = 1, verbose = FALSE),
                 "deprecated and will be removed in fqardl 2.0.0")
  expect_no_warning(mtnardl(gdp ~ oil_price, data = oil_gdp_data, max_p = 1,
                            max_q = 1, verbose = FALSE))
  st$mtnardl_warned <- FALSE
  w <- tryCatch(mtnardl(gdp ~ oil_price, data = oil_gdp_data, max_p = 1,
                        max_q = 1, verbose = FALSE),
                warning = function(w) conditionMessage(w))
  expect_match(w, "ardlverse::mtnardl")
  expect_match(w, "p - 1")
  expect_match(w, "computed incorrectly when p or q exceeded 1")
})

test_that("mtnardl bounds F uses only y_lag1 and the level regressors (F1)", {
  set.seed(1)
  n <- 150
  x <- cumsum(rnorm(n)); y <- numeric(n)
  for (t in 3:n) y[t] <- y[t - 1] - 0.3 * (y[t - 1] - 0.5 * x[t - 1]) +
    0.2 * (y[t - 1] - y[t - 2]) + rnorm(1, sd = 0.5)
  dec <- decompose_multi_threshold(x, 0)
  X <- dec$components
  colnames(X) <- paste0("x_", dec$names)
  fit <- fqardl:::estimate_mtnardl(y, X, p = 2, q = 2, case = 3)
  b <- fqardl:::perform_mtnardl_bounds(fit, n, 2, 3)
  # hand computation on an independently built design: p = 2, q = 2,
  # sample t = 3..n, partial sums of the negative and positive changes
  dx <- c(0, diff(x))
  xn <- cumsum(pmin(dx, 0)); xp <- cumsum(pmax(dx, 0))
  tt <- 3:n
  dy  <- y[tt] - y[tt - 1]
  ly  <- y[tt - 1]; dy1 <- y[tt - 1] - y[tt - 2]
  lxn <- xn[tt - 1]; lxp <- xp[tt - 1]
  dxn0 <- xn[tt] - xn[tt - 1]; dxn1 <- xn[tt - 1] - xn[tt - 2]
  dxp0 <- xp[tt] - xp[tt - 1]; dxp1 <- xp[tt - 1] - xp[tt - 2]
  ur <- lm(dy ~ ly + dy1 + lxn + lxp + dxn0 + dxn1 + dxp0 + dxp1)
  rr <- lm(dy ~ dy1 + dxn0 + dxn1 + dxp0 + dxp1)
  expect_equal(b$F_stat, anova(rr, ur)$F[2], tolerance = 1e-10)
  expect_equal(unname(fit$coefficients["y_lag1"]), unname(coef(ur)["ly"]),
               tolerance = 1e-10)
  expect_equal(b$k, 2L)
  expect_equal(b$cv_5, c(3.79, 4.85))
  # the 1.0.6 statistic also tested dy_lag1 and the lagged differences
  old <- fqardl:::wald_bounds_F(fit$coefficients, fit$vcov,
                                grep("_lag1$", names(fit$coefficients), value = TRUE))
  expect_false(isTRUE(all.equal(b$F_stat, old)))
})

test_that("mtnardl refuses cases other than 3 (F3)", {
  data(oil_gdp_data, package = "fqardl")
  expect_error(suppressWarnings(mtnardl(gdp ~ oil_price, data = oil_gdp_data,
                                        case = 1, verbose = FALSE)),
               "only case = 3")
  expect_error(suppressWarnings(mtnardl(gdp ~ oil_price, data = oil_gdp_data,
                                        case = 4, verbose = FALSE)),
               "only case = 3")
})

test_that("print and summary are registered for fqardl_mtnardl only", {
  reg <- getNamespaceInfo("fqardl", "S3methods")
  cls <- reg[reg[, 1] %in% c("print", "summary"), 2]
  expect_true("fqardl_mtnardl" %in% cls)
  expect_false("mtnardl" %in% cls)
  expect_identical(sum(cls == "fqardl_mtnardl"), 2L)  # print and summary
  data(oil_gdp_data, package = "fqardl")
  m <- suppressWarnings(mtnardl(gdp ~ oil_price, data = oil_gdp_data,
                                max_p = 1, max_q = 1, verbose = FALSE))
  expect_identical(class(m), c("fqardl_mtnardl", "mtnardl"))
  expect_output(summary(m), "REGIME-SPECIFIC LONG-RUN MULTIPLIERS")
})

test_that("plot() on an fqardl mtnardl object stops with an informative error", {
  data(oil_gdp_data, package = "fqardl")
  m <- suppressWarnings(mtnardl(gdp ~ oil_price, data = oil_gdp_data,
                                max_p = 1, max_q = 1, verbose = FALSE))
  expect_error(plot(m), "ardlverse::mtnardl")
  reg <- getNamespaceInfo("fqardl", "S3methods")
  expect_true(any(reg[, 1] == "plot" & reg[, 2] == "fqardl_mtnardl"))
})
