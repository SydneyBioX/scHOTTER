#' Convert a weight matrix to a 0/1 membership matrix
#'
#' Turns any nonzero weight into \code{1} and zeros into \code{0}, preserving the
#' sparse structure and dimension names. This is useful for assigning points to
#' kernels based on whether their weight is positive (regardless of magnitude).
#' This is the intermediate step required to generate the list of overlapping
#' kernels.
#'
#' @param weight_matrix A numeric matrix or a \code{Matrix::dgCMatrix} of weights
#'   with rows as points and columns as kernels (e.g., from
#'   \code{\link{generate_weight_matrix_euclidean}}).
#'
#' @return A \code{Matrix::dgCMatrix} with the same dimensions and dimnames as
#'   \code{weight_matrix}, containing only \code{0/1} entries.
#'
#' @examples
#' set.seed(1)
#' coords  <- cbind(runif(30), runif(30))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
#' Mb      <- generate_membership_matrix(W)
#' Matrix::nnzero(Mb)      # number of assigned point-kernel memberships
#' unique(Mb@x)            # should be c(1) (no stored zeros)
#'
#' @export
generate_membership_matrix <- function(weight_matrix) {
  # Coerce to a compressed sparse column matrix (dgCMatrix)
  if (is.data.frame(weight_matrix)) weight_matrix <- as.matrix(weight_matrix)
  if (is.matrix(weight_matrix)) {
    weight_matrix <- Matrix::Matrix(weight_matrix, sparse = TRUE)
    weight_matrix <- methods::as(weight_matrix, "dgCMatrix")
  } else if (!inherits(weight_matrix, "dgCMatrix")) {
    weight_matrix <- methods::as(weight_matrix, "dgCMatrix")
  }

  # Binary-ize the stored values: any nonzero -> 1, zero -> 0
  Wb <- weight_matrix
  if (length(Wb@x)) {
    Wb@x <- as.numeric(Wb@x > 0)
    Wb <- Matrix::drop0(Wb)  # remove any explicit zeros from the sparse slots
  }
  Wb
}
