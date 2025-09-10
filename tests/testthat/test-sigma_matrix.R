test_that("sigma matrix has expected diagonal and symmetry", {
  set.seed(2)
  n <- 80
  coords  <- cbind(runif(n), runif(n))
  centres <- generate_kernel_centres_by_density(coords, span = 0.2)
  W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")

  Memb  <- generate_membership_matrix(W)
  ovl   <- determine_overlaps(Memb)
  light <- find_light_kernels(ovl, Memb)
  Wt    <- trim_weight_matrix(W, light)

  neff  <- get_effective_sample_sizes(Wt)
  corr  <- approximate_between_coefficient_correlations_effective(Wt)  # unattenuated
  Sigma <- get_sigma_matrix(corr, neff)

  expect_true(Matrix::isSymmetric(Sigma))
  # Check diagonal where neff > 3
  ok <- which(neff > 3 & is.finite(neff))
  if (length(ok)) {
    expect_equal(
      as.numeric(Matrix::diag(Sigma)[ok]),
      unname(1 / (neff[ok] - 3)),
      tolerance = 1e-12
    )
  }
  # Where neff <= 3, diagonal should be NA
  bad <- which(!(neff > 3 & is.finite(neff)))
  if (length(bad)) expect_true(all(is.na(Matrix::diag(Sigma)[bad])))
  # Names are preserved
  expect_setequal(colnames(Sigma), colnames(corr))
})
