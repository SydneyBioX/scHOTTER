#' Determine overlapping kernels (share \eqn{\geq} 1 point)
#'
#' From a binary membership matrix \eqn{M} (\eqn{\text{points} \times \text{kernels}}), compute which
#' kernels overlap, where overlap means at least one shared point. Internally
#' uses \eqn{A = M^\top M}; \eqn{A_{ij} > 0} indicates that kernels \eqn{i} and \eqn{j}
#' share at least one point. This is an intermediate step required to determine
#' which kernels are "light" and shouldn't be considered in downstream analysis.
#'
#' @param membership_matrix A 0/1 matrix (preferably \code{Matrix::dgCMatrix})
#'   with rows = points and cols = kernels. You can obtain this via
#'   \code{\link{generate_membership_matrix}}.
#'
#' @return A named list of length \eqn{k} (number of kernels). Element \code{[[i]]}
#'   is a character vector of kernel names that overlap with kernel \code{i}
#'   **including itself**.
#'
#' @examples
#' # Minimal controlled example
#' M <- Matrix::sparseMatrix(i = c(1,2,2,3), j = c(1,1,2,2), x = 1,
#'                           dims = c(3,2), dimnames = list(NULL, c("k_1","k_2")))
#' determine_overlaps(M)
#'
#' # Typical workflow
#' set.seed(1)
#' coords  <- cbind(runif(40), runif(40))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "block")
#' Memb    <- generate_membership_matrix(W)
#' ovl     <- determine_overlaps(Memb)
#' str(ovl[1:2])
#'
#' @export
determine_overlaps <- function(membership_matrix){
  # Coerce to sparse dgCMatrix if needed
  if (is.data.frame(membership_matrix)) membership_matrix <- as.matrix(membership_matrix)
  if (is.matrix(membership_matrix)) {
    membership_matrix <- Matrix::Matrix(membership_matrix, sparse = TRUE)
    membership_matrix <- methods::as(membership_matrix, "dgCMatrix")
  } else if (!inherits(membership_matrix, "dgCMatrix")) {
    membership_matrix <- methods::as(membership_matrix, "dgCMatrix")
  }

  k <- ncol(membership_matrix)
  if (k == 0L) return(stats::setNames(vector("list", 0L), character(0)))

  kernel_names <- colnames(membership_matrix)
  if (is.null(kernel_names)) kernel_names <- paste0("k_", seq_len(k))

  # A = t(M) %*% M; logical overlap when A > 0
  overlap_indicator <- Matrix::crossprod(membership_matrix) > 0

  overlapping_kernels_list <- vector("list", k)
  for (i in seq_len(k)) {
    idx <- which(overlap_indicator[i, , drop = TRUE])
    overlapping_kernels_list[[i]] <- kernel_names[idx]
  }
  names(overlapping_kernels_list) <- kernel_names
  overlapping_kernels_list
}
