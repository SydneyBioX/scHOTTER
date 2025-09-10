test_that("weight matrix builds and has sensible dims", {
  set.seed(1)
  coords <- cbind(runif(60, 0, 10), runif(60, 0, 5))
  centres <- generate_kernel_centres_by_density(coords, span = 0.1)
  Wg <- generate_weight_matrix_euclidean(coords, centres, span = 0.1, type = "gaussian")
  Wb <- generate_weight_matrix_euclidean(coords, centres, span = 0.1, type = "block")

  expect_s4_class(Wg, "dgCMatrix")
  expect_equal(dim(Wg), c(nrow(coords), nrow(centres)))
  expect_true(Matrix::nnzero(Wg) > 0)

  expect_s4_class(Wb, "dgCMatrix")
  expect_true(all(Wb@x %in% c(0, 1))) # block weights are 0/1
})
