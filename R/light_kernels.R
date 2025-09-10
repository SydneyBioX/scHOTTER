#' Identify "light" kernels with too few unique points
#'
#' Flags kernels that don't have enough assigned points relative to how many
#' other kernels they overlap with. The criterion used is
#' \code{num_points < num_overlaps}, where \code{num_points} is the column sum
#' of the (0/1) membership matrix and \code{num_overlaps} is the length of the
#' overlap set for that kernel (typically including itself).
#'
#' @param overlapping_list A named list where each element \code{overlapping_list[[k]]}
#'   is a character vector of kernel names that overlap kernel \code{k}. See
#'   \code{\link{determine_overlaps}}.
#' @param membership_matrix A 0/1 matrix (preferably \code{Matrix::dgCMatrix})
#'   of point–kernel memberships with rows = points, cols = kernels.
#'
#' @return A character vector of kernel names deemed "light".
#'
#' @details
#' This function assumes the overlap list includes each kernel itself (as produced
#' by \code{determine_overlaps()}). If you supply an overlap list that excludes
#' self, the threshold will be one smaller; adjust upstream if needed.
#'
#' @examples
#' set.seed(1)
#' coords  <- cbind(runif(40), runif(40))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "block")
#' Memb    <- generate_membership_matrix(W)
#' ovl     <- determine_overlaps(Memb)
#' find_light_kernels(ovl, Memb)
#'
#' @export
find_light_kernels <- function(overlapping_list, membership_matrix){
  # Coerce membership to sparse dgCMatrix if needed
  if (is.data.frame(membership_matrix)) membership_matrix <- as.matrix(membership_matrix)
  if (is.matrix(membership_matrix)) {
    membership_matrix <- Matrix::Matrix(membership_matrix, sparse = TRUE)
    membership_matrix <- methods::as(membership_matrix, "dgCMatrix")
  } else if (!inherits(membership_matrix, "dgCMatrix")) {
    membership_matrix <- methods::as(membership_matrix, "dgCMatrix")
  }

  kernel_names <- colnames(membership_matrix)
  if (is.null(kernel_names)) kernel_names <- paste0("k_", seq_len(ncol(membership_matrix)))

  # Points per kernel (fast for sparse matrices)
  num_points <- Matrix::colSums(membership_matrix)

  # Overlaps per kernel (align to kernel_names; missing names -> length 0)
  num_overlaps <- lengths(overlapping_list[kernel_names])

  # "Light" kernels: points <= overlaps
  kernel_names[num_points < num_overlaps]
}
