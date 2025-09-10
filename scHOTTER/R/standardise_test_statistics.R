#' Standardize sample-variance test statistics into Z-scores
#'
#' Converts the raw sample variance of local Fisher z correlations
#' (one value per gene pair) into standardized Z-scores using the
#' null mean and variance derived from the kernel dependence structure. The resulting
#' standardised statistic is distributed approximately according to the standard
#' Normal distribution.
#'
#' @param variances Numeric scalar or vector: \eqn{\mathrm{Var}(S^2)} under the null.
#'   Typically the scalar output of \code{approximate_variance_effective(Sigma)}.
#'   Can also be a named vector matching \code{vars_of_local_corrs}.
#' @param expectations Numeric scalar or vector: \eqn{\mathbb{E}[S^2]} under the null.
#'   Typically the scalar output of \code{approximate_expectation_effective(Sigma)}.
#'   Can also be a named vector matching \code{vars_of_local_corrs}.
#' @param vars_of_local_corrs Named numeric vector: observed sample variances
#'   (one per gene pair), e.g. from \code{compute_var_of_local_corrs()}.
#'
#' @return A named numeric vector of Z-scores aligned with
#'   \code{names(vars_of_local_corrs)}. Entries with non-finite or non-positive
#'   variances are returned as \code{NA_real_}.
#'
#' @examples
#' # Simple algebraic example
#' s2  <- c(a_b = 3, a_c = 5, b_c = 4)
#' mu  <- 4      # E[S^2]
#' var <- 1.0    # Var(S^2)
#' compute_test_statistics_standardised_effective(var, mu, s2)
#'
#' @export
compute_test_statistics_standardised_effective <- function(variances, expectations, vars_of_local_corrs){
  # Basic checks
  if (!is.numeric(vars_of_local_corrs)) stop("vars_of_local_corrs must be numeric.")
  gp <- names(vars_of_local_corrs)
  if (is.null(gp)) gp <- paste0("pair_", seq_len(length(vars_of_local_corrs)))

  n <- length(vars_of_local_corrs)

  # Helper to align scalar/vector to the gene-pair names
  align_vec <- function(x, target_names, label){
    if (length(x) == 1L) {
      rep_len(as.numeric(x), length.out = length(target_names))
    } else if (length(x) == length(target_names)) {
      if (!is.null(names(x))) {
        m <- match(target_names, names(x))
        if (anyNA(m)) stop(sprintf("%s has names that do not cover all gene pairs.", label))
        as.numeric(x[m])
      } else {
        as.numeric(x)
      }
    } else {
      stop(sprintf("%s must be length 1 or length(vars_of_local_corrs).", label))
    }
  }

  mu  <- align_vec(expectations, gp, "expectations")
  v   <- align_vec(variances,    gp, "variances")

  denom <- sqrt(v)
  bad   <- !is.finite(denom) | denom <= 0

  z <- (as.numeric(vars_of_local_corrs) - mu) / denom
  z[bad] <- NA_real_
  names(z) <- gp
  z
}
