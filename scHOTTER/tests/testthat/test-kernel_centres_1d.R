test_that("1D kernel centres align along the varying axis and fix others", {
  set.seed(1)
  t  <- sort(runif(50, 0, 10))
  y0 <- rep(3, 50)
  Z  <- cbind(t = t, y = y0)

  C <- generate_kernel_centres_by_density_1d(Z, span = 0.2)
  expect_true(is.matrix(C))
  expect_equal(colnames(C), colnames(Z))
  # varying axis is "t": centres monotone within [min,max]
  expect_true(all(diff(C[, "t"]) >= 0))
  expect_gte(min(C[, "t"]), min(t))
  expect_lte(max(C[, "t"]), max(t))
  # constant axis is at (robust) median
  expect_true(all(abs(C[, "y"] - median(y0)) < 1e-12))
})
