test_that("pipeline returns p-values and intermediates with expected shapes", {
  set.seed(7)
  n <- 40
  coords <- cbind(runif(n), runif(n))
  rn <- apply(coords, 1, function(v) paste0(v[1], "_", v[2]))
  expr <- cbind(g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
  rownames(expr) <- rn

  out <- scHOTTER_pipeline(expr, span = 0.2, kernel_type = "block")
  expect_true(is.numeric(out$p_values))
  expect_equal(names(out$p_values), out$summary$pair)
  expect_true(all(c("coords","centres","weight_matrix","local_corr_matrix","neff","sigma") %in%
                    names(out$intermediates)))
})
