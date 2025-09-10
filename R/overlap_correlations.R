#' Approximate correlation between kernel-wise Fisher coefficients
#'
#' Computes an overlap-based correlation between kernels from a point×kernel
#' weight matrix \eqn{W}. First, columns are normalized to
#' \eqn{A = W D^{-1}}, where \eqn{D = \mathrm{diag}(\sum_\ell w_{ k \ell})}.
#' Then \eqn{\Omega = A^{\top} A} measures weighted overlap, and we form the
#' cosine-similarity-style correlation
#' \deqn{C = D_\Omega^{-1/2}\, \Omega \, D_\Omega^{-1/2},}
#' where \eqn{D_\Omega = \mathrm{diag}(\Omega)}.
#'
#' Optionally, if effective sample sizes \code{neff} are provided, we apply the
#' Fisher-\eqn{z} attenuation factor \eqn{s_k=\sqrt{(n_{\mathrm{eff},k}-3)/(n_{\mathrm{eff},k}-1)}}
#' (for \eqn{n_{\mathrm{eff},k}>3}) by
#' \deqn{C \leftarrow S\, C \, S, \quad S=\mathrm{diag}(s_1,\dots,s_K).}
#' The diagonal is set to 1 for kernels with nonzero overlap norm and \code{NA}
#' for empty/degenerate kernels.
#'
#' @param weight_matrix \eqn{\text{Point} \times \text{kernel}} weight matrix (preferably \code{Matrix::dgCMatrix}).
#'
#' @return A symmetric numeric matrix \eqn{K \times K} with row/column names equal
#'   to kernel names. Diagonal entries are \code{1} when the kernel has positive
#'   norm (nonempty) and \code{NA} otherwise.
#'
#' @examples
#' set.seed(1)
#' coords  <- cbind(runif(50), runif(50))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
#' C1      <- approximate_between_coefficient_correlations_effective(W)
#' neff    <- get_effective_sample_sizes(W)
#' C2      <- approximate_between_coefficient_correlations_effective(W)
#' range(C1, na.rm = TRUE)
#'
#' @export
approximate_between_coefficient_correlations_effective <- function(weight_matrix){
  # Coerce to sparse dgCMatrix if needed
  W <- weight_matrix
  if (is.data.frame(W)) W <- as.matrix(W)
  if (is.matrix(W)) {
    W <- Matrix::Matrix(W, sparse = TRUE)
    W <- methods::as(W, "dgCMatrix")
  } else if (!inherits(W, "dgCMatrix")) {
    W <- methods::as(W, "dgCMatrix")
  }

  k <- ncol(W)
  kernel_names <- colnames(W)
  if (is.null(kernel_names)) kernel_names <- paste0("k_", seq_len(k))
  if (k == 0L) {
    C_empty <- matrix(numeric(0), nrow = 0, ncol = 0,
                      dimnames = list(character(0), character(0)))
    return(C_empty)
  }

  # Column-normalize: A = W %*% diag(1/sum_j w_jk)
  s1 <- Matrix::colSums(W)
  inv_s1 <- ifelse(s1 > 0, 1 / s1, 0)
  A <- W %*% Matrix::Diagonal(x = inv_s1)

  # Overlap matrix and its diagonal
  omega <- Matrix::crossprod(A)               # K x K, symmetric
  d <- Matrix::diag(omega)
  d[d < 0] <- 0
  d_sqrt <- sqrt(d)
  inv_d <- ifelse(d_sqrt > 0, 1 / d_sqrt, 0)

  # Cosine-similarity style correlation
  C <- Matrix::Diagonal(x = inv_d) %*% omega %*% Matrix::Diagonal(x = inv_d)

  # Set diagonal: 1 for nondegenerate kernels, NA for zero-norm columns
  di <- Matrix::diag(C)
  nondeg <- d_sqrt > 0
  di[nondeg] <- 1
  di[!nondeg] <- NA_real_
  Matrix::diag(C) <- di

  dimnames(C) <- list(kernel_names, kernel_names)
  C
}
