test_that("compute_var_of_local_corrs returns named vector and handles NAs", {
  set.seed(1)
  Z <- matrix(rnorm(5 * 3), nrow = 5, dimnames = list(NULL, c("a_b","a_c","b_c")))

  # Make column a_b have < 2 finite values -> variance should be NA
  Z[1:4, "a_b"] <- NA_real_

  s2 <- compute_var_of_local_corrs(Z)

  expect_type(s2, "double")
  expect_setequal(names(s2), colnames(Z))

  # Compare against manual var; drop names on the 1-element vector
  expect_equal(unname(s2["b_c"]), stats::var(Z[, "b_c"], na.rm = TRUE))

  # Column with fewer than 2 finite values -> NA
  expect_true(is.na(s2["a_b"]))
})
