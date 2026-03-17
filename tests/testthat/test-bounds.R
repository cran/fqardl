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
  
  # Upper should be > lower
  expect_true(all(cv$F_upper > cv$F_lower))
  
  # Critical values should be positive for F
  expect_true(all(cv$F_lower > 0))
  expect_true(all(cv$F_upper > 0))
  
  # t critical values should be negative
  expect_true(all(cv$t_upper < 0))
})

test_that("critical values vary with k", {
  cv_k1 <- get_pss_critical_values(k = 1, case = 3)
  cv_k3 <- get_pss_critical_values(k = 3, case = 3)
  
  # Critical values should generally decrease with k
  expect_true(cv_k1$F_lower[2] > cv_k3$F_lower[2])
})

test_that("critical values vary with case", {
  cv_case3 <- get_pss_critical_values(k = 2, case = 3)
  cv_case4 <- get_pss_critical_values(k = 2, case = 4)
  
  # Case 4 (with trend) should have higher critical values
  expect_true(cv_case4$F_lower[2] > cv_case3$F_lower[2])
})
