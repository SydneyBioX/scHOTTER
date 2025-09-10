test_that("overlap correlations have unit diag (nondegenerate) and symmetry", {
  set.seed(1)
  coords  <- cbind(runif(40), runif(40))
  centres <- generate_kernel_centres_by_density(coords, span = 0.2)
  W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "block")

  Memb  <- generate_membership_matrix(W)
  ovl   <- determine_overlaps(Memb)
  light <- find_light_kernels(ovl, Memb)
  Wt    <- trim_weight_matrix(W, light)

  C <- approximate_between_coefficient_correlations_effective(Wt)

  expect_true(Matrix::isSymmetric(C))
  d <- Matrix::diag(C)
  expect_true(all(is.na(d) | abs(d - 1) < 1e-12))

  C2 <- approximate_between_coefficient_correlations_effective(Wt)
  expect_true(Matrix::isSymmetric(C2))
  d2 <- Matrix::diag(C2)
  expect_true(all(is.na(d2) | abs(d2 - 1) < 1e-12))
})
