# Tests for FQARDL main function

test_that("fqardl handles basic estimation", {
  skip_on_cran()  # Skip on CRAN due to computation time
  
  set.seed(123)
  n <- 100
  
  # Generate simple cointegrated data
  x1 <- cumsum(rnorm(n))
  y <- 0.5 * x1 + rnorm(n, sd = 0.5)
  
  data <- data.frame(y = y, x1 = x1)
  
  result <- fqardl(
    formula = y ~ x1,
    data = data,
    tau = c(0.5),
    max_k = 2,
    max_p = 2,
    max_q = 2,
    criterion = "BIC",
    bootstrap = FALSE
  )
  
  # Check class
  expect_s3_class(result, "fqardl")
  
  # Check components exist
  expect_true(!is.null(result$optimal_k))
  expect_true(!is.null(result$optimal_p))
  expect_true(!is.null(result$optimal_q))
  expect_true(!is.null(result$long_run))
  expect_true(!is.null(result$bounds_test))
})

test_that("fqardl rejects invalid inputs", {
  data <- data.frame(y = 1:10, x = 1:10)
  
  # Invalid formula
  expect_error(fqardl("not a formula", data))
  
  # Invalid tau
  expect_error(fqardl(y ~ x, data, tau = c(0, 0.5)))
  expect_error(fqardl(y ~ x, data, tau = c(0.5, 1)))
})

test_that("fqardl print and summary work", {
  skip_on_cran()
  
  set.seed(123)
  n <- 80
  x1 <- cumsum(rnorm(n))
  y <- 0.5 * x1 + rnorm(n, sd = 0.5)
  data <- data.frame(y = y, x1 = x1)
  
  result <- fqardl(y ~ x1, data, tau = 0.5, max_k = 1, max_p = 1, max_q = 1)
  
  # Print should not error
  expect_output(print(result))
  expect_output(summary(result))
})
