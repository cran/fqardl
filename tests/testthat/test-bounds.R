# Tests for bounds testing functions

test_that("get_pss_critical_values returns correct structure", {
  cv <- get_pss_critical_values(k = 2, case = 3)
  
  # Check all components exist
  expect_true(!is.null(cv$F_lower))
  expect_true(!is.null(cv$F_upper))
  expect_true(!is.null(cv$t_upper))
  
  # Check lengths (1%, 5%, 10%)
  expect_equal(length(cv$F_lower), 3)
  expect_equal(length(cv$F_upper), 3)
  expect_equal(length(cv$t_upper), 3)

  # CHANGED in 1.0.4: only the verified 5 percent column is supplied; the 1 and
  # 10 percent points are NA rather than numbers of unknown provenance.
  expect_true(cv$F_upper[2] > cv$F_lower[2])
  expect_true(cv$F_lower[2] > 0)
  expect_true(cv$F_upper[2] > 0)
  expect_true(all(is.na(c(cv$F_lower[c(1,3)], cv$F_upper[c(1,3)]))))
  
  # t critical values should be negative.
  # CHANGED in 1.0.3: only the 5 percent I(1) bound has been verified against
  # Pesaran, Shin and Smith Table CII(iii). The 1 and 10 percent bounds are
  # returned as NA rather than guessed, so the assertion is on the verified
  # entry only.
  expect_true(cv$t_upper[2] < 0)
  expect_true(all(cv$t_lower < 0))
  expect_true(all(is.na(cv$t_upper[c(1, 3)])))
})

test_that("critical values vary with k", {
  cv_k1 <- get_pss_critical_values(k = 1, case = 3)
  cv_k3 <- get_pss_critical_values(k = 3, case = 3)
  
  # Critical values should generally decrease with k
  expect_true(cv_k1$F_lower[2] > cv_k3$F_lower[2])
  expect_true(cv_k1$F_upper[2] > cv_k3$F_upper[2])
})

test_that("cases other than 3 are refused rather than silently mis-served", {
  # CHANGED in 1.0.3. Up to 1.0.2 the `case` argument switched the critical
  # value table while the design matrix stayed Case III, so cases 4 and 5
  # compared a Case III model against Case IV and Case V bounds. Handling the
  # cases properly means changing the regressors, which is scheduled for 2.0.0;
  # until then an unsupported case is an error, not a wrong number.
  expect_error(get_pss_critical_values(k = 2, case = 4), "Only case = 3")
  expect_error(get_pss_critical_values(k = 2, case = 5), "Only case = 3")
  expect_silent(get_pss_critical_values(k = 2, case = 3))
})
