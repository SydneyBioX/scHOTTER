test_that("approximate_variance_effective matches known special cases", {
  # Independent equal-variance case: Var(S^2) = 2 v^2 / (n-1)
  for (n in c(2, 5, 10)) {
    v <- 0.3
    Sigma <- diag(rep(v, n))
    got <- approximate_variance_effective(Sigma)
    expect_equal(got, 2 * v^2 / (n - 1), tolerance = 1e-12)
  }

  # Handles NA rows/cols by dropping them
  SigmaNA <- diag(c(0.2, 0.2, NA_real_, 0.2))
  got_na <- approximate_variance_effective(SigmaNA)
  expect_true(is.finite(got_na) && got_na > 0)

  # Returns NA when <2 kernels remain after dropping
  SigmaTiny <- matrix(c(0.1, NA, NA, 0.1), nrow = 2)
  expect_true(is.na(approximate_variance_effective(SigmaTiny)))
})
