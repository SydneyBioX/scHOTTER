test_that("approximate_expectation_effective matches special cases and handles NA", {
  # Independent equal-variance case: E[S^2] = v
  for (n in c(2, 5, 10)) {
    v <- 0.3
    Sigma <- diag(rep(v, n))
    got <- approximate_expectation_effective(Sigma)
    expect_equal(got, v, tolerance = 1e-12)
  }

  # Handles NA rows/cols by dropping them (still finite if >=2 remain)
  SigmaNA <- diag(c(0.2, 0.25, NA_real_, 0.3))
  got_na <- approximate_expectation_effective(SigmaNA)
  expect_true(is.finite(got_na) && got_na > 0)

  # Returns NA when <2 kernels remain after dropping
  SigmaTiny <- matrix(c(0.1, NA, NA, 0.2), nrow = 2)
  expect_true(is.na(approximate_expectation_effective(SigmaTiny)))
})
