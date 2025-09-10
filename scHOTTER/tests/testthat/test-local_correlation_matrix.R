test_that("matrix version matches vector version for each pair", {
  set.seed(42)
  n <- 60
  expr <- cbind(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  coords <- cbind(runif(n), runif(n))
  centres <- generate_kernel_centres_by_density(coords, span = 0.2)
  W <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")

  Zmat <- get_local_correlation_matrix(expr, W, fisher_transform = TRUE, drop_zero_weight = TRUE)

  pairs <- colnames(Zmat)
  # check a couple of random pairs (or all if you prefer)
  for (lbl in sample(pairs, min(3, length(pairs)))) {
    zv <- get_local_correlation_vector_pair(expr, W, lbl, fisher_transform = TRUE, drop_zero_weight = TRUE)
    expect_equal(Zmat[, lbl], zv, tolerance = 1e-12)
  }
})
