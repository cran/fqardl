# Tests for FNARDL functions

test_that("decompose_variables works correctly", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  
  data <- data.frame(x = x)
  result <- decompose_variables(data, "x")
  
  # Check structure
  expect_true("x" %in% names(result))
  expect_true(all(c("positive", "negative", "d_positive", "d_negative") %in% names(result$x)))
  
  # Check lengths
  expect_equal(length(result$x$positive), n)
  expect_equal(length(result$x$negative), n)
  
  # Positive should be non-negative cumsum
  expect_true(all(diff(result$x$positive) >= 0))
  
  # Negative should be non-positive cumsum
  expect_true(all(diff(result$x$negative) <= 0))
  
  # Sum should equal original (approximately)
  reconstructed <- result$x$positive + result$x$negative
  expect_equal(reconstructed, x - x[1], tolerance = 1e-10)
})

test_that("fnardl handles basic estimation", {
  skip_on_cran()
  
  set.seed(456)
  n <- 100
  
  # Generate asymmetric relationship
  oil <- cumsum(rnorm(n, sd = 2))
  gdp <- numeric(n)
  gdp[1] <- 100
  for (i in 2:n) {
    d_oil <- oil[i] - oil[i-1]
    if (d_oil > 0) {
      gdp[i] <- 0.95 * gdp[i-1] - 0.3 * d_oil + rnorm(1, sd = 0.5)
    } else {
      gdp[i] <- 0.95 * gdp[i-1] - 0.7 * d_oil + rnorm(1, sd = 0.5)
    }
  }
  
  data <- data.frame(gdp = gdp, oil = oil)
  
  result <- fnardl(
    formula = gdp ~ oil,
    data = data,
    decompose = "oil",
    max_k = 1,
    max_p = 2,
    max_q = 2,
    bootstrap = FALSE
  )
  
  # Check class
  expect_s3_class(result, "fnardl")
  
  # Check asymmetry test exists
  expect_true(!is.null(result$asymmetry_tests))
  expect_true("oil" %in% names(result$asymmetry_tests))
  
  # Check multipliers
  expect_true(!is.null(result$multipliers$long_run_positive))
  expect_true(!is.null(result$multipliers$long_run_negative))
})

test_that("test_asymmetry returns valid results", {
  skip_on_cran()
  
  set.seed(789)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- cumsum(rnorm(n))
  
  data <- data.frame(y = y, x = x)
  decomposed <- decompose_variables(data, "x")
  X <- cbind(decomposed$x$positive, decomposed$x$negative)
  colnames(X) <- c("x_pos", "x_neg")
  
  fourier <- generate_fourier_terms(n, 1)
  nardl_result <- estimate_nardl(y, X, fourier, 2, 2, 3)
  
  asym_test <- test_asymmetry(nardl_result, "x")
  
  # Check structure
  expect_true("x" %in% names(asym_test))
  expect_true(!is.null(asym_test$x$wald_stat))
  expect_true(!is.null(asym_test$x$p_value))
  expect_true(asym_test$x$p_value >= 0 && asym_test$x$p_value <= 1)
})
