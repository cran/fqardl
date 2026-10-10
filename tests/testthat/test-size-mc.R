## Size of the recursive bootstrap under the null of no cointegration
## (independent random walks). Slow, therefore not run on CRAN. Fixed lags and
## Fourier frequency for speed; the Monte Carlo table in NEWS.md (MC 200,
## B 199, n 100) uses the defaults with re-selection in every replication.

test_that("fqardl bootstrap has size close to 5 percent", {
  skip_on_cran()
  MC <- 60
  n <- 100
  rej <- matrix(NA, MC, 4, dimnames = list(NULL, c("Fov", "t", "Find", "dec")))
  fr <- generate_fourier_terms(n, 1)
  for (r in seq_len(MC)) {
    set.seed(5000 + r)
    y <- cumsum(rnorm(n)); x <- cumsum(rnorm(n))
    b <- bootstrap_bounds_test(y, cbind(x1 = x), fr, p = 1, q = 1, n_boot = 99,
                               seed = r)
    rej[r, ] <- c(b$p_value_F < 0.05, b$p_value_t < 0.05,
                  b$p_value_Find < 0.05, b$decision == "COINTEGRATION")
  }
  size <- colMeans(rej)
  expect_lt(size["Fov"], 0.15)
  expect_lt(size["t"], 0.15)
  expect_lt(size["dec"], 0.10)
})

test_that("fnardl bootstrap has size close to 5 percent", {
  skip_on_cran()
  MC <- 60
  n <- 100
  rej <- matrix(NA, MC, 4, dimnames = list(NULL, c("Fov", "t", "Find", "dec")))
  fr <- generate_fourier_terms(n, 1)
  for (r in seq_len(MC)) {
    set.seed(6000 + r)
    y <- cumsum(rnorm(n)); x <- cumsum(rnorm(n))
    b <- fqardl:::bootstrap_nardl(y, cbind(x = x), "x", "x", fr, p = 1, q = 1,
                                  n_boot = 99, seed = r)
    rej[r, ] <- c(b$p_value_F < 0.05, b$p_value_t < 0.05,
                  b$p_value_Find < 0.05, b$decision == "COINTEGRATION")
  }
  size <- colMeans(rej)
  expect_lt(size["Fov"], 0.15)
  expect_lt(size["t"], 0.15)
  expect_lt(size["dec"], 0.10)
})
