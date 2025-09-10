#' Trim kernels flagged as "light" from a weight matrix
#'
#' Removes columns (kernels) listed in \code{light_kernels} from a point×kernel
#' weight matrix. Returns the original matrix if no kernels are specified.
#'
#' @param weight_matrix A numeric matrix or \code{Matrix::dgCMatrix} with rows =
#'   points, cols = kernels (e.g., from \code{generate_weight_matrix_euclidean()}).
#'   Must have column names.
#' @param light_kernels Character vector of kernel names to remove. If empty,
#'   the input matrix is returned unchanged.
#'
#' @return A \code{Matrix::dgCMatrix} with the specified columns removed (possibly
#'   with zero columns if all kernels are removed).
#'
#' @examples
#' set.seed(1)
#' coords  <- cbind(runif(40), runif(40))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
#' Memb    <- generate_membership_matrix(W)
#' ovl     <- determine_overlaps(Memb)
#' light   <- find_light_kernels(ovl, Memb)
#' W_trim  <- trim_weight_matrix(W, light)
#' dim(W_trim)
#'
#' @export
trim_weight_matrix <- function(weight_matrix, light_kernels){
  # Coerce to sparse dgCMatrix if needed
  if (is.data.frame(weight_matrix)) weight_matrix <- as.matrix(weight_matrix)
  if (is.matrix(weight_matrix)) {
    weight_matrix <- Matrix::Matrix(weight_matrix, sparse = TRUE)
    weight_matrix <- methods::as(weight_matrix, "dgCMatrix")
  } else if (!inherits(weight_matrix, "dgCMatrix")) {
    weight_matrix <- methods::as(weight_matrix, "dgCMatrix")
  }

  if (length(light_kernels) == 0L) return(weight_matrix)
  if (is.null(colnames(weight_matrix))) {
    stop("weight_matrix must have column names to match against light_kernels.")
  }

  lk <- unique(as.character(light_kernels[!is.na(light_kernels)]))
  keep <- is.na(match(colnames(weight_matrix), lk))

  weight_matrix[, keep, drop = FALSE]
}
