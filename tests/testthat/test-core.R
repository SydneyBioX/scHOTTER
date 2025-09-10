test_that("pipeline runs end-to-end and returns sane shapes", {
  testthat::skip_on_cran()

  set.seed(11)
  n <- 30
  coords <- cbind(runif(n, 0, 10), runif(n, 0, 5))
  rn <- apply(coords, 1, function(v) paste0(v[1], "_", v[2]))
  expr <- cbind(g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
  rownames(expr) <- rn

  out <- scHOTTER_pipeline(expr, span = 0.2, kernel_type = "gaussian")

  # p-values & summary line up
  expect_true(is.numeric(out$p_values))
  expect_setequal(names(out$p_values), out$summary$pair)

  # key intermediates present with expected types
  ii <- out$intermediates
  expect_s4_class(ii$weight_matrix_initial, "dgCMatrix")
  expect_s4_class(ii$weight_matrix, "dgCMatrix")
  expect_true(is.matrix(ii$local_corr_matrix) || is.null(ii$local_corr_matrix))
  expect_true(is.numeric(ii$neff))
  expect_true(is.matrix(ii$corr) || inherits(ii$corr, "Matrix"))
  expect_true(is.matrix(ii$sigma) || inherits(ii$sigma, "Matrix"))
})
