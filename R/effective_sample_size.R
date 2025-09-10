#' Effective sample size per kernel (Kish, 1965)
#'
#' Computes Kish's effective sample size for each kernel/column of a weight
#' matrix: \deqn{n_{\mathrm{eff},k} = \frac{(\sum_\ell w_{ k \ell})^2}
#' {\sum_\ell w_{k \ell }^2}.}
#'
#' @param weight_mat_trimmed A point×kernel weight matrix (preferably
#'   \code{Matrix::dgCMatrix}); can be the output of
#'   \code{\link{trim_weight_matrix}}.
#' @param tol Numeric tolerance. Columns with \eqn{\sum w^2 \leq } \code{tol}
#'   are treated as empty and return \code{0}.
#'
#' @return A named numeric vector of length \eqn{k} (kernels), giving
#'   \eqn{n_{\mathrm{eff}}} for each kernel (0 for empty kernels).
#'
#' @details
#' For block weights (0/1), \eqn{n_{\mathrm{eff},k}} equals the number of
#' positively weighted points in kernel \eqn{k}. For general nonnegative weights,
#' \eqn{n_{\mathrm{eff},k} \leq \sum_\ell 1\{w_{k \ell}>0\}}.
#'
#' @references
#' Kish, L. (1965). \emph{Survey Sampling}. Wiley.
#'
#' @examples
#' set.seed(1)
#' coords  <- cbind(runif(40), runif(40))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' Wg <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
#' neff <- get_effective_sample_sizes(Wg)
#' head(neff)
#'
#' @export
get_effective_sample_sizes <- function(weight_mat_trimmed, tol = 1e-12){
  # Coerce to sparse dgCMatrix if needed
  W <- weight_mat_trimmed
  if (is.data.frame(W)) W <- as.matrix(W)
  if (is.matrix(W)) {
    W <- Matrix::Matrix(W, sparse = TRUE)
    W <- methods::as(W, "dgCMatrix")
  } else if (!inherits(W, "dgCMatrix")) {
    W <- methods::as(W, "dgCMatrix")
  }

  # Sums per column
  s1 <- Matrix::colSums(W)

  # Sum of squares per column, done sparsely by squaring the nonzero slots
  if (length(W@x)) {
    W2 <- W
    W2@x <- W2@x^2
    s2 <- Matrix::colSums(W2)
  } else {
    s2 <- rep.int(0, ncol(W))
  }

  # Handle empties robustly:
  # - if s1 ~ 0 and s2 ~ 0 -> n_eff = 0 (empty kernel)
  # - otherwise n_eff = (s1^2)/s2
  neff <- rep(NA_real_, ncol(W))
  empty <- s2 <= tol | s1 <= tol
  neff[empty] <- 0
  ok <- !empty
  neff[ok] <- (s1[ok]^2) / s2[ok]

  names(neff) <- colnames(W)
  neff
}
