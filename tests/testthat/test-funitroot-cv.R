test_that("Fourier ADF critical values follow Enders and Lee (2012) by k and T", {
  expect_equal(fqardl:::get_fadf_critical_values(100, "c", 1), c(-4.42, -3.81, -3.49))
  expect_equal(fqardl:::get_fadf_critical_values(100, "c", 2), c(-3.97, -3.27, -2.91))
  expect_equal(fqardl:::get_fadf_critical_values(200, "ct", 3), c(-4.38, -3.77, -3.43))
  expect_equal(fqardl:::get_fadf_critical_values(1000, "c", 5), c(-3.53, -2.93, -2.61))
  expect_equal(fqardl:::get_fadf_f_critical_values(100, "c"), c(10.35, 7.58, 6.35))
})

test_that("Fourier KPSS critical values follow Becker, Enders and Lee (2006)", {
  expect_equal(fqardl:::get_fkpss_critical_values("c", 2, 100), c(0.6671, 0.4152, 0.3150))
  expect_equal(fqardl:::get_fkpss_critical_values("ct", 1, 100), c(0.0716, 0.0546, 0.0471))
})

test_that("no approximate p-values are reported for the Fourier unit root tests", {
  set.seed(1)
  y <- cumsum(rnorm(120))
  r <- fourier_adf_test(y, max_freq = 2, verbose = FALSE)
  expect_true(is.na(r$p_value))
  expect_true(is.na(r$f_test$p_value))
  expect_equal(r$critical_values, fqardl:::get_fadf_critical_values(120, "c", r$optimal_frequency))
  k <- fourier_kpss_test(y, max_freq = 2, verbose = FALSE)
  expect_true(is.na(k$p_value))
})
