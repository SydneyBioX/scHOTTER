#' Upper-tail p-values for standardized test statistics
#'
#' Computes one-sided p-values \eqn{P(Z > z)} assuming the input test statistics
#' are standard normal under the null. This is appropriate when larger values of
#' the statistic (e.g., a standardized sample variance) indicate stronger evidence
#' against the null.
#'
#' @param test_stats Named numeric vector of Z-scores (e.g., from
#'   \code{compute_test_statistics_standardised_effective()}).
#'
#' @return A named numeric vector of upper-tail p-values in \eqn{[0,1]}.
#'   Non-finite inputs (NA/Inf) return \code{NA_real_}.
#'
#' @examples
#' s2  <- c(a_b = 3, a_c = 5, b_c = 4)
#' mu  <- 4
#' var <- 1
#' z   <- compute_test_statistics_standardised_effective(var, mu, s2)
#' p   <- compute_p_values(z)
#' p
#'
#' @export
compute_p_values <- function(test_stats){
  if (!is.numeric(test_stats)) stop("test_stats must be numeric.")
  nm <- names(test_stats)
  z  <- as.numeric(test_stats)
  # p = P(Z > z) under N(0,1)
  p_values <- stats::pnorm(z, lower.tail = FALSE)
  # preserve names if present
  if (!is.null(nm)) names(p_values) <- nm
  p_values[!is.finite(z)] <- NA_real_
  p_values
}
