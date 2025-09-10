#' Variance of the sample variance of Fisher z's under dependence
#'
#' Computes \eqn{\mathrm{Var}(S^2)} for the unbiased sample variance
#' \eqn{S^2 = \frac{1}{n-1}\sum_{i=1}^n (Z_i - \bar{\boldsymbol{Z}})^2} when the
#' kernel-wise Fisher z vector \eqn{Z} has covariance matrix \eqn{\boldsymbol{\Sigma}}.
#' The formula used is
#' \deqn{
#' \mathrm{Var}(S^2)=\frac{2}{(n-1)^2}\left(
#' \mathrm{tr}(\boldsymbol{\Sigma}^{2}) - \frac{2}{n}\,\mathbf{1}^{\top}\boldsymbol{\Sigma}^{2}\mathbf{1}
#' + \frac{(\mathbf{1}^{\top}\boldsymbol{\Sigma}\mathbf{1})^2}{n^2}
#' \right),
#' }
#' which reduces to \eqn{2\sigma^4/(n-1)} for independent components with common
#' variance \eqn{\sigma^2}.
#'
#' @param sigma_matrix A symmetric \eqn{n\times n} covariance matrix for the
#'   kernel-wise Fisher z coefficients (e.g., from \code{get_sigma_matrix()}).
#'   May be a base matrix or a \pkg{Matrix} object. Rows/columns containing
#'   \code{NA}s are dropped prior to computation.
#'
#' @return A single numeric: \eqn{\mathrm{Var}(S^2)}. Returns \code{NA_real_}
#'   if, after dropping \code{NA} rows/cols, fewer than 2 kernels remain.
#'
#' @examples
#' # Sanity check: independent equal-variance case
#' set.seed(1)
#' n <- 10
#' v <- 0.2
#' Sigma <- diag(rep(v, n))
#' approximate_variance_effective(Sigma)  # ~ 2 * v^2 / (n-1)
#'
#' @export
approximate_variance_effective <- function(sigma_matrix){
  # Coerce to a Matrix object (dense) and symmetrize lightly
  S <- sigma_matrix
  if (is.data.frame(S)) S <- as.matrix(S)
  if (!inherits(S, "Matrix")) S <- Matrix::Matrix(S, sparse = FALSE)
  # Drop rows/cols with any NA (e.g., neff <= 3 kernels)
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

  # tr(S^2) == ||S||_F^2 for symmetric S
  tr_S2 <- Matrix::norm(S, "F")^2

  # u = S * 1, so u'u = 1' S^2 1 and s = 1' S 1
  ones <- rep.int(1, n)
  u  <- as.numeric(S %*% ones)
  u2 <- sum(u * u)
  s  <- sum(u)

  trace_term <- tr_S2 - (2 / n) * u2 + (s * s) / (n * n)
  var_S2 <- (2 / (n - 1)^2) * trace_term

  as.numeric(var_S2)
}
