test_that("find_light_kernels returns expected names", {
  # Build a sparse membership matrix:
  # rows = points, cols = k1 k2 k3
  # - k1 has 3 points (rows 1, 2, 5), k2 has 3 points (rows 3, 4, 6), k3 has 1 point (row 2)
  # - overlaps: k1<->k3 via row 2; k1 and k2 do NOT overlap; k2 and k3 do NOT overlap
  M <- Matrix::sparseMatrix(
    i = c(1, 2,      3, 4, 6,      2),   # rows
    j = c(1, 1,      2, 2, 2,      3),   # cols: 1=k1,2=k2,3=k3
    x = 1,
    dims = c(6, 3),
    dimnames = list(NULL, c("k1", "k2", "k3"))
  )
  # Add a unique point to k1 so it has 3 points
  M[5, "k1"] <- 1

  ovl <- determine_overlaps(M)  # includes self

  # Points per kernel: k1=3, k2=3, k3=1
  # Overlaps per kernel (incl. self): k1=2 (k1,k3), k2=1 (k2), k3=2 (k1,k3)
  # With strict rule "<": only k3 is light (1 < 2)
  expect_setequal(find_light_kernels(ovl, M), "k3")

  # Now give k3 one additional unique point (so points == overlaps == 2)
  new_row <- Matrix::sparseMatrix(
    i = 1, j = 3, x = 1,
    dims = c(1, 3),
    dimnames = list(NULL, colnames(M))
  )
  M2 <- rbind(M, new_row)   # base rbind works with Matrix objects
  ovl2 <- determine_overlaps(M2)

  # k3 now has 2 points and still overlaps with exactly k1 (and itself) -> not light under "<"
  expect_setequal(find_light_kernels(ovl2, M2), character(0))
})
