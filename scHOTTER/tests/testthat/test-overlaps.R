test_that("determine_overlaps returns named list incl. self and is symmetric", {
  M <- Matrix::sparseMatrix(i = c(1,2,2,3), j = c(1,1,2,2), x = 1,
                            dims = c(3,2), dimnames = list(NULL, c("k_1","k_2")))
  ovl <- determine_overlaps(M)

  expect_type(ovl, "list")
  expect_equal(names(ovl), c("k_1","k_2"))
  # includes self
  expect_true(all(c("k_1","k_2") %in% ovl[["k_1"]]))
  expect_true(all(c("k_1","k_2") %in% ovl[["k_2"]]))
  # symmetry
  expect_true("k_2" %in% ovl[["k_1"]] && "k_1" %in% ovl[["k_2"]])
})
