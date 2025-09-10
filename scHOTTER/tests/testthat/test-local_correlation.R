test_that("get_local_correlation_vector_pair returns named numeric vector", {
  set.seed(1)
  n <- 50
  expr <- cbind(g1 = rnorm(n), g2 = rnorm(n))
  coords <- cbind(runif(n), runif(n))
  centres <- generate_kernel_centres_by_density(coords, span = 0.2)
  W <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "block")

  z <- get_local_correlation_vector_pair(expr, W, "g1_g2")
  expect_type(z, "double")
  expect_equal(length(z), ncol(W))
  expect_setequal(names(z), colnames(W))
  # Fisher z is finite only when |r| < 1
  expect_true(all(is.na(z) | is.finite(z)))
})
