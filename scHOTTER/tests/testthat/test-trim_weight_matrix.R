test_that("trim_weight_matrix removes specified kernels and preserves dims", {
  set.seed(1)
  coords  <- cbind(runif(20), runif(20))
  centres <- generate_kernel_centres_by_density(coords, span = 0.25)
  W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.25, type = "block")

  # No trimming -> identical dims
  W0 <- trim_weight_matrix(W, character(0))
  expect_equal(dim(W0), dim(W))

  # Trim a couple of columns by name
  kn <- colnames(W)[1:min(2, ncol(W))]
  W1 <- trim_weight_matrix(W, kn)
  expect_equal(ncol(W1), ncol(W) - length(kn))
  expect_false(any(kn %in% colnames(W1)))

  # Unknown names are ignored
  W2 <- trim_weight_matrix(W, c("not_a_kernel", kn[1]))
  expect_equal(ncol(W2), ncol(W) - 1L)

  # Trimming all columns yields n x 0 matrix
  W_all <- trim_weight_matrix(W, colnames(W))
  expect_equal(dim(W_all)[2], 0L)
  expect_s4_class(W_all, "dgCMatrix")
})
