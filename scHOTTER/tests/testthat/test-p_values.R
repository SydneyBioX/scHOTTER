test_that("compute_p_values returns upper-tail one-sided p-values", {
  z <- c(a = -1, b = 0, c = 1, d = 2, e = NA_real_)
  p <- compute_p_values(z)

  expect_type(p, "double")
  expect_equal(names(p), names(z))

  # Monotone decreasing in z
  expect_true(p["a"]  > p["b"] && p["b"] > p["c"] && p["c"] > p["d"])

  # Known values
  expect_equal(unname(p["b"]), 0.5, tolerance = 1e-12)    # P(Z>0)=0.5
  expect_equal(unname(p["d"]), 1 - pnorm(2), tolerance = 1e-12)

  # NA propagation
  expect_true(is.na(p["e"]))
})
