# Tests for Fourier functions

test_that("generate_fourier_terms works correctly", {
  n <- 100
  k <- 2
  
  fourier <- generate_fourier_terms(n, k)
  
  # Check dimensions
expect_equal(nrow(fourier), n)
  expect_equal(ncol(fourier), 2)  # sin and cos
  
  # Check column names
  expect_true(all(c("sin_2", "cos_2") %in% colnames(fourier)))
  
  # Check values are bounded
  expect_true(all(fourier >= -1 & fourier <= 1))
})

test_that("generate_fourier_terms with cumulative = TRUE", {
  n <- 100
  k <- 3
  
  fourier <- generate_fourier_terms(n, k, cumulative = TRUE)
  
  # Should have 2*k columns
  expect_equal(ncol(fourier), 2 * k)
})

test_that("generate_fourier_terms rejects invalid k", {
  expect_error(generate_fourier_terms(100, 0))
  expect_error(generate_fourier_terms(100, -1))
  expect_error(generate_fourier_terms(100, 1.5))
})

test_that("select_fourier_frequency returns valid k", {
  set.seed(123)
  n <- 100
  y <- cumsum(rnorm(n))
  X <- matrix(cumsum(rnorm(n)), ncol = 1)
  
  result <- select_fourier_frequency(y, X, max_k = 3, criterion = "BIC")
  
  expect_true(result$optimal_k >= 1)
  expect_true(result$optimal_k <= 3)
  expect_equal(length(result$ic_values), 3)
})
