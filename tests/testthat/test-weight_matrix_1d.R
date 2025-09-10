test_that("1D weight matrix builds with expected dims and sparsity", {
  set.seed(1)
  n <- 60
  t  <- sort(runif(n, 0, 10))
  y0 <- rep(0, n)   # constant axis
  coords  <- cbind(t = t, y = y0)
  centres <- generate_kernel_centres_by_density_1d(coords, span = 0.2)

  Wb <- generate_weight_matrix_euclidean_1d(coords, centres, span = 0.2, type = "block")
  Wg <- generate_weight_matrix_euclidean_1d(coords, centres, span = 0.2, type = "gaussian")

  expect_s4_class(Wb, "dgCMatrix")
  expect_s4_class(Wg, "dgCMatrix")
  expect_equal(dim(Wb), c(nrow(coords), nrow(centres)))
  expect_equal(dim(Wg), c(nrow(coords), nrow(centres)))
  expect_true(Matrix::nnzero(Wb) > 0)
  expect_true(Matrix::nnzero(Wg) > 0)
})
