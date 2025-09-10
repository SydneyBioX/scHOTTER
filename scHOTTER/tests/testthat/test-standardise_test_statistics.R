test_that("standardisation works for scalar and vector inputs, with name alignment", {
  s2 <- c(a_b = 3, a_c = 5, b_c = 4)

  # scalar mu/var
  z1 <- compute_test_statistics_standardised_effective(variances = 1, expectations = 4, vars_of_local_corrs = s2)
  expect_equal(unname(z1), c(-1, 1, 0))

  # vector mu/var with names (per-pair) and shuffled order
  mu_vec  <- c(b_c = 4, a_b = 4, a_c = 4)
  var_vec <- c(a_c = 1, b_c = 1, a_b = 1)
  z2 <- compute_test_statistics_standardised_effective(var_vec, mu_vec, s2)
  expect_equal(z2[names(s2)], z1[names(s2)])

  # non-positive variance -> NA
  z3 <- compute_test_statistics_standardised_effective(variances = 0, expectations = 4, vars_of_local_corrs = s2)
  expect_true(all(is.na(z3)))
})
