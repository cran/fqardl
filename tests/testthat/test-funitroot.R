# Tests for Fourier Unit Root functions

test_that("fourier_adf_test works on stationary series", {
  set.seed(123)
  n <- 200
  t <- 1:n
  
  # Stationary with Fourier component
  y <- 5 + 3 * sin(2 * pi * t / n) + rnorm(n, sd = 0.5)
  
  result <- fourier_adf_test(y, model = "c", max_freq = 3)
  
  expect_s3_class(result, "fadf")
  expect_true(result$reject_null)  # Should reject unit root
  expect_true(result$optimal_frequency >= 1)
  expect_true(result$optimal_frequency <= 3)
})

test_that("fourier_adf_test works on unit root series", {
  set.seed(456)
  n <- 200
  
  # Random walk (unit root)
  y <- cumsum(rnorm(n))
  
  result <- fourier_adf_test(y, model = "c", max_freq = 3)
  
  expect_s3_class(result, "fadf")
  expect_false(result$reject_null)  # Should NOT reject unit root
})

test_that("fourier_kpss_test works correctly", {
  set.seed(789)
  n <- 200
  t <- 1:n
  
  # Stationary series
  y <- 5 + 2 * cos(2 * pi * t / n) + rnorm(n, sd = 0.3)
  
  result <- fourier_kpss_test(y, model = "c", max_freq = 3)
  
  expect_s3_class(result, "fkpss")
  expect_false(result$reject_null)  # Should NOT reject stationarity
  expect_true(result$statistic > 0)
})

test_that("fadf_f_test correctly identifies need for Fourier terms", {
  set.seed(111)
  n <- 200
  t <- 1:n
  
  # Series with strong Fourier component (breaks)
  y <- 5 + 5 * sin(2 * pi * t / n) + rnorm(n, sd = 0.5)
  
  result <- fadf_f_test(y, model = "c", k = 1, p = 2)
  
  expect_true(result$reject)  # Should reject linearity (Fourier terms needed)
  expect_true(result$f_stat > result$critical_values[2])  # Above 5% CV
})

test_that("fourier_unit_root_analysis provides joint conclusion", {
  set.seed(222)
  n <- 150
  
  # Stationary series
  y <- rnorm(n)
  
  result <- fourier_unit_root_analysis(y, name = "Test", max_freq = 2)
  
  expect_true("adf" %in% names(result))
  expect_true("kpss" %in% names(result))
  expect_true("conclusion" %in% names(result))
})

test_that("get_fadf_critical_values returns correct structure", {
  cv <- get_fadf_critical_values(200, "c", 1)
  
  expect_length(cv, 3)  # 1%, 5%, 10%
  expect_true(cv[1] < cv[2])  # More negative at 1%
  expect_true(cv[2] < cv[3])
})

test_that("get_fkpss_critical_values returns correct structure", {
  cv <- get_fkpss_critical_values("c", 1)
  
  expect_length(cv, 3)
  expect_true(cv[1] > cv[2])  # Higher at 1%
  expect_true(cv[2] > cv[3])
})
