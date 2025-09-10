test_that("Kish n_eff works after trimming (block, gaussian, empties)", {
  set.seed(123)
  n <- 50
  coords  <- cbind(runif(n), runif(n))
  centres <- generate_kernel_centres_by_density(coords, span = 0.2)

  ## --- BLOCK weights ---
  Wb <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "block")
  Memb_b  <- generate_membership_matrix(Wb)
  ovl_b   <- determine_overlaps(Memb_b)
  light_b <- find_light_kernels(ovl_b, Memb_b)
  Wb_trim <- trim_weight_matrix(Wb, light_b)

  neff_b <- get_effective_sample_sizes(Wb_trim)

  if (ncol(Wb_trim) > 0) {
    counts_b <- Matrix::colSums(generate_membership_matrix(Wb_trim))  # 0/1, exact match
    expect_equal(unname(neff_b), as.numeric(counts_b))
  } else {
    expect_length(neff_b, 0L)
  }

  ## --- GAUSSIAN weights ---
  Wg <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
  Memb_g  <- generate_membership_matrix(Wg)
  ovl_g   <- determine_overlaps(Memb_g)
  light_g <- find_light_kernels(ovl_g, Memb_g)
  Wg_trim <- trim_weight_matrix(Wg, light_g)

  neff_g <- get_effective_sample_sizes(Wg_trim)

  if (ncol(Wg_trim) > 0) {
    pos_counts_g <- Matrix::colSums(generate_membership_matrix(Wg_trim))
    expect_true(all(neff_g <= as.numeric(pos_counts_g) + 1e-8))
  } else {
    expect_length(neff_g, 0L)
  }

  ## --- Empty column returns 0 (post-trim) ---
  if (ncol(Wb_trim) >= 1) {
    W0 <- Wb_trim
    W0[, 1] <- 0
    neff0 <- get_effective_sample_sizes(W0)
    expect_equal(unname(neff0[1]), 0)
  }
})
