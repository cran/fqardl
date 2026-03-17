# Tests for Multi-Threshold NARDL functions

test_that("decompose_multi_threshold creates correct regimes", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n, sd = 2))
  
  result <- decompose_multi_threshold(x, c(-1, 0, 1))
  
  # Should have 4 regimes for 3 thresholds
  expect_equal(ncol(result$components), 4)
  expect_length(result$names, 4)
  expect_true(all(result$thresholds == c(-1, 0, 1)))
})

test_that("decompose_multi_threshold with single threshold equals NARDL", {
  set.seed(456)
  n <- 50
  x <- cumsum(rnorm(n))
  
  result <- decompose_multi_threshold(x, 0)
  
  # Should have 2 regimes (positive and negative)
  expect_equal(ncol(result$components), 2)
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
  
  result <- mtnardl(
    formula = y ~ x,
    data = data,
    decompose = "x",
    thresholds = list(x = c(-1, 0, 1)),
    max_p = 2,
    max_q = 2
  )
  
  expect_s3_class(result, "mtnardl")
  expect_length(result$regime_names$x, 4)
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
  
  result <- mtnardl(
    formula = y ~ x,
    data = data,
    decompose = "x",
    thresholds = list(x = c(0)),
    max_p = 1,
    max_q = 1
  )
  
  expect_output(print(result), "Multi-Threshold NARDL")
})
