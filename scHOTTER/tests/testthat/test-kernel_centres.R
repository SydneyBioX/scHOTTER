test_that("kernel centres shape and type are sensible", {
  set.seed(1)
  coords <- cbind(runif(50, 0, 10), runif(50, 0, 5))
  centres <- generate_kernel_centres_by_density(coords, span = 0.1)
  expect_true(is.matrix(centres))
  expect_identical(colnames(centres), c("x","y"))
  expect_true(nrow(centres) >= 1)
  expect_true(all(is.finite(centres)))
})
