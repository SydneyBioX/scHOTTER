#' Sample variance of kernel-wise Fisher z per gene pair
#'
#' Given a matrix of local Fisher z correlations (rows = kernels, columns = gene
#' pairs like \code{"GENE1_GENE2"}), compute the unbiased sample variance
#' \eqn{S^2} down each column, ignoring \code{NA}s. Columns with fewer than two
#' finite values return \code{NA_real_}.
#'
#' @param local_corr_matrix Numeric matrix with rows = kernels and columns =
#'   gene-pair labels (e.g., from \code{get_local_correlation_matrix()} with
#'   \code{fisher_transform = TRUE}).
#'
#' @return A named numeric vector of length \code{ncol(local_corr_matrix)}
#'   containing the sample variance for each gene pair.
#'
#' @examples
#' set.seed(1)
#' expr <- cbind(a = rnorm(60), b = rnorm(60), c = rnorm(60))
#' coords <- cbind(runif(60), runif(60))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
#' Z <- get_local_correlation_matrix(expr, W, fisher_transform = TRUE)
#' s2 <- compute_var_of_local_corrs(Z)
#' head(s2)
#'
#' @export
compute_var_of_local_corrs <- function(local_corr_matrix){
  # Coerce & validate
  M <- local_corr_matrix
  if (is.data.frame(M)) M <- as.matrix(M)
  if (!is.matrix(M)) stop("local_corr_matrix must be a numeric matrix.")
  if (is.null(colnames(M))) {
    colnames(M) <- paste0("pair_", seq_len(ncol(M)))
  }

  safe_var <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    stats::var(x)  # unbiased sample variance (denominator n-1)
  }

  out <- vapply(seq_len(ncol(M)), function(j) safe_var(M[, j]), numeric(1))
  names(out) <- colnames(M)
  out
}
