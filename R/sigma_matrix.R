#' Covariance matrix of kernel-wise Fisher z coefficients
#'
#' Given a kernel–kernel correlation matrix \code{corr} (for the Fisher
#' coefficients) and effective sample sizes \code{neff}, returns the covariance
#' matrix \eqn{\boldsymbol{\Sigma} = \boldsymbol{D}\,\mathrm{corr}\,\boldsymbol{D}}
#' where \eqn{\boldsymbol{D}=\mathrm{diag}(\sqrt{\boldsymbol{v}})}
#' and \eqn{v_k = 1/(n_{\mathrm{eff},k}-3)} (variance of Fisher z for a correlation).
#'
#' @param corr A symmetric \eqn{n \times n} correlation matrix between kernels
#'   (e.g., from \code{approximate_between_coefficient_correlations_effective(W)}).
#'   Preferably a \code{Matrix} object; a base matrix is also accepted.
#'   **Pass the unattenuated correlation** (call with \code{neff = NULL}) so you
#'   don't double-apply any scaling.
#' @param neff Numeric vector of length \eqn{n} with effective sample sizes per
#'   kernel (e.g., from \code{get_effective_sample_sizes(W)}).
#'
#' @return A \code{Matrix} object of size \eqn{n \times n} with row/column names
#'   taken from \code{corr}. Entries corresponding to kernels with
#'   \eqn{n_{\mathrm{eff}} \leq 3} are returned as \code{NA_real_}.
#'
#' @details
#' For kernels with \eqn{n_{\mathrm{eff}} \leq 3}, the Fisher z variance is undefined.
#' In that case, we set the corresponding diagonal multiplier to zero while
#' forming \eqn{\boldsymbol{D}\,\mathrm{corr}\,\boldsymbol{D}}, and then mark the affected rows/columns of
#' \eqn{\boldsymbol{\Sigma}} as \code{NA}.
#'
#' @examples
#' set.seed(1)
#' coords  <- cbind(runif(60), runif(60))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W       <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
#' neff    <- get_effective_sample_sizes(W)
#' corr    <- approximate_between_coefficient_correlations_effective(W)  # unattenuated
#' Sigma   <- get_sigma_matrix(corr, neff)
#' Matrix::diag(Sigma)[1:5]
#'
#' @export
get_sigma_matrix <- function(corr, neff){
  # Coerce corr to a Matrix object (dense) for stable ops
  C <- corr
  if (is.data.frame(C)) C <- as.matrix(C)
  if (!inherits(C, "Matrix")) C <- Matrix::Matrix(C, sparse = FALSE)

  k <- ncol(C)
  if (length(neff) != k)
    stop("length(neff) must equal ncol(corr).")

  kn <- colnames(C)
  if (is.null(kn)) kn <- paste0("k_", seq_len(k))

  # v_k = 1/(neff_k - 3), undefined when neff <= 3
  v <- rep(NA_real_, k)
  ok <- is.finite(neff) & (neff > 3)
  v[ok] <- 1 / (neff[ok] - 3)

  s <- sqrt(v)  # may contain NAs where neff <= 3

  # Build D with zeros where s is NA to avoid propagating NaNs in the product
  s_fill <- ifelse(is.na(s), 0, s)
  D <- Matrix::Diagonal(x = s_fill)

  Sigma <- D %*% C %*% D

  # Mark rows/cols with undefined variance as NA
  if (anyNA(s)) {
    bad <- which(is.na(s))
    if (length(bad)) {
      Sigma[bad, ] <- NA_real_
      Sigma[, bad] <- NA_real_
    }
  }

  dimnames(Sigma) <- list(kn, kn)
  Sigma
}
