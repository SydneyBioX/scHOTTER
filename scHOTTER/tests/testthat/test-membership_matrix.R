test_that("membership matrix is binary, sparse, and preserves dims", {
  set.seed(1)
  coords  <- cbind(runif(20), runif(20))
  centres <- generate_kernel_centres_by_density(coords, span = 0.25)
  Wg <- generate_weight_matrix_euclidean(coords, centres, span = 0.25, type = "gaussian")

  Mb <- generate_membership_matrix(Wg)

  expect_s4_class(Mb, "dgCMatrix")
  expect_equal(dim(Mb), dim(Wg))
  expect_true(all(Mb@x %in% c(1)))  # only 1's stored; zeros are dropped
  expect_true(Matrix::nnzero(Mb) <= Matrix::nnzero(Wg))
})
