#' Expectation of the sample variance of Fisher z's under dependence
#'
#' Computes \eqn{\mathbb{E}[S^2]} for the unbiased sample variance
#' \eqn{S^2 = \frac{1}{n-1}\sum_{i=1}^n (Z_i - \bar{\boldsymbol{Z}})^2} when the
#' kernel-wise Fisher z vector \eqn{\boldsymbol{Z}} has covariance matrix
#' \eqn{\boldsymbol{\Sigma}}. The formula used is
#' \deqn{
#' \mathbb{E}[S^2]=\frac{1}{n-1}\left(
#' \mathrm{tr}(\boldsymbol{H} \boldsymbol{\Sigma})
#' \right), \quad
#' \boldsymbol{H} = \boldsymbol{I}_{n} -
#' \frac{1}{n}\boldsymbol{1}\boldsymbol{1}^{\top}
#' }
#' where \eqn{\boldsymbol{H}} is the centreing matrix.
#' For independent components with common variance \eqn{\sigma^2}, this reduces
#' to \eqn{\mathbb{E}[S^2]=\sigma^2}.
#'
#' @param sigma_matrix A symmetric \eqn{n\times n} covariance matrix for the
#'   kernel-wise Fisher z coefficients (e.g., from \code{get_sigma_matrix()}).
#'   May be a base matrix or a \pkg{Matrix} object. Rows/columns containing
#'   \code{NA}s are dropped prior to computation.
#'
#' @return A single numeric: \eqn{\mathbb{E}[S^2]}. Returns \code{NA_real_}
#'   if, after dropping \code{NA} rows/cols, fewer than 2 kernels remain.
#'
#' @examples
#' # Sanity check: independent equal-variance case
#' n <- 8; v <- 0.2
#' Sigma <- diag(rep(v, n))
#' approximate_expectation_effective(Sigma)  # = v
#'
#' @export
approximate_expectation_effective <- function(sigma_matrix){
  # Coerce to a Matrix object (dense) for stable ops
  S <- sigma_matrix
  if (is.data.frame(S)) S <- as.matrix(S)
  if (!inherits(S, "Matrix")) S <- Matrix::Matrix(S, sparse = FALSE)

  # Drop rows/cols with any NA (e.g., from neff <= 3 kernels)
  Smat <- as.matrix(S)
  keep <- which(rowSums(!is.na(Smat)) == ncol(Smat))
  if (length(keep) < nrow(S)) {
    if (length(keep) < 2L) return(NA_real_)
    S <- Matrix::Matrix(Smat[keep, keep, drop = FALSE], sparse = FALSE)
  }

  # Ensure symmetry numerically
  S <- (S + Matrix::t(S)) / 2

  n <- nrow(S)
  if (n < 2L) return(NA_real_)

  tr_S <- sum(Matrix::diag(S))
  ones <- rep.int(1, n)
  s    <- sum(as.numeric(S %*% ones))   # 1' S 1

  exp_S2 <- (1 / (n - 1)) * (tr_S - s / n)
  as.numeric(exp_S2)
}
