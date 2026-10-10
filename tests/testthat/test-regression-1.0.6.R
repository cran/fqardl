## Values computed with fqardl 1.0.6 that 1.1.0 must reproduce exactly:
## QARDL point estimates, long-run multipliers and the cross-quantile Wald
## statistic. The quantile standard errors use quantreg's bootstrap, so the
## Wald statistic depends on the seed set by fqardl(seed = ) (here 2024 via
## set.seed) and on quantreg (6.1 when pinned).

test_that("QARDL estimates and the CKS Wald test equal the 1.0.6 values", {
  data(macro_data, package = "fqardl")
  set.seed(2024)
  f <- fqardl(gdp ~ inflation + interest_rate, data = macro_data,
              tau = c(0.25, 0.5, 0.75), max_k = 2, max_p = 2, max_q = 2,
              verbose = FALSE)
  expect_identical(c(f$optimal_k, f$optimal_p, f$optimal_q), c(1L, 1L, 1L))
  b <- f$qardl_results[[2]]$coefficients
  expect_equal(unname(b), c(0.5384352795125634, -0.0992543504490277,
                            0.0837071848921190, -0.1030259942517435,
                            0.0216285767006430, -0.0345360426129133,
                            0.0252694471282358, 0.0915684467755738),
               tolerance = 1e-10)
  expect_equal(unname(f$long_run[, "inflation"]),
               c(1.524954290220998, 0.843360361671068, 1.133178007546106),
               tolerance = 1e-10)
  expect_equal(unname(f$long_run[, "interest_rate"]),
               c(-0.584907831829649, -1.037999783240259, -1.476627004719809),
               tolerance = 1e-10)
  skip_if_not(packageVersion("quantreg") == "6.1", "pinned with quantreg 6.1")
  w <- suppressWarnings(quantile_wald_test(f$qardl_results, "x1_lag1"))
  expect_equal(c(w$wald_stat, w$p_value), c(0.761750403588096, 0.683263154315200),
               tolerance = 1e-8)
})
