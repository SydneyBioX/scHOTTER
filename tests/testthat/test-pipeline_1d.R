test_that("1D pipeline runs end-to-end and returns sane shapes", {
  set.seed(7)
  n <- 60
  t  <- sort(runif(n, 0, 10))
  y0 <- rep(0, n)
  coords <- cbind(t = t, y = y0)
  rn <- apply(coords, 1, function(v) paste0(v[1], "_", v[2]))
  expr <- cbind(g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
  rownames(expr) <- rn

  out <- scHOTTER_pipeline_1d(expr, span = 0.2, kernel_type = "gaussian")

  expect_true(is.numeric(out$p_values))
  expect_setequal(names(out$p_values), out$summary$pair)

  ii <- out$intermediates
  expect_s4_class(ii$weight_matrix_initial, "dgCMatrix")
  expect_s4_class(ii$weight_matrix, "dgCMatrix")
  expect_true(is.matrix(ii$local_corr_matrix) || is.null(ii$local_corr_matrix))
  expect_true(is.numeric(ii$neff))
  expect_true(Matrix::isSymmetric(ii$corr))
  expect_true(Matrix::isSymmetric(ii$sigma))
})
