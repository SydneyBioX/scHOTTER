#' End-to-end pipeline: local co-expression p-values (stores intermediates)
#'
#' Takes a gene expression matrix whose \strong{rownames encode 2D coordinates}
#' (e.g. "12.3x45.6", "12.3_45.6", "12.3,45.6", or "12.3 45.6") and returns one-sided
#' upper-tail p-values for the standardized sample-variance statistic of
#' kernel-wise Fisher z correlations, \emph{plus} all intermediate objects.
#'
#' @param expr Numeric matrix (\eqn{\text{points} \times  \text{genes}}). \strong{Row names must contain}
#'   "x y" coordinates separated by `x`, `_`, `,`, or whitespace (e.g. `"10.2_5.7"`).
#'   Column names are gene symbols.
#' @param span Optional numeric in (0,1]; target proportion of points per kernel
#'   (passed to kernel generation and weighting). Default uses the package
#'   heuristic if `NULL`.
#' @param kernel_type `"block"` or `"gaussian"`; passed to
#'   [generate_weight_matrix_euclidean()].
#' @param fisher_transform Logical; Fisher-transform local correlations (`atanh`).
#'   Must be `TRUE` to match downstream variance formulas (default `TRUE`).
#' @param drop_zero_weight Logical; restrict each kernel's correlation to rows
#'   with \code{w>0} for that kernel (default `TRUE`).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{p_values}: named numeric vector of p-values (one per gene pair).
#'   \item \code{summary}: data.frame with columns \code{pair}, \code{s2}, \code{z}, \code{p}.
#'   \item \code{intermediates}: a list containing
#'     \code{coords}, \code{centres}, \code{weight_matrix_initial}, \code{membership},
#'     \code{overlaps}, \code{light_kernels}, \code{weight_matrix}, \code{local_corr_matrix},
#'     \code{neff}, \code{corr}, \code{sigma}, \code{expectation}, \code{variance},
#'     \code{s2}, \code{z_scores}, \code{p_values}.
#'   \item \code{params}: the input parameters actually used.
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 80
#' coords <- cbind(runif(n, 0, 10), runif(n, 0, 5))
#' rn <- apply(coords, 1, function(v) paste0(v[1], "_", v[2]))
#' expr <- cbind(g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
#' rownames(expr) <- rn
#' out <- scHOTTER_pipeline(expr, span = 0.2, kernel_type = "gaussian")
#' head(out$summary)
#'
#' @export
scHOTTER_pipeline <- function(expr,
                              span = NULL,
                              kernel_type = c("block", "gaussian"),
                              fisher_transform = TRUE,
                              drop_zero_weight = TRUE) {
  kernel_type <- match.arg(kernel_type)

  # ---- validate & parse inputs ----
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || is.null(colnames(expr)))
    stop("expr must be a numeric matrix with gene column names.")
  if (is.null(rownames(expr)))
    stop("expr must have rownames encoding coordinates like 'x_y' or 'x,y'.")

  coords <- .parse_coords_from_rownames(rownames(expr))
  if (!is.matrix(coords) || ncol(coords) != 2 || any(!is.finite(coords))) {
    stop("Failed to parse numeric 2D coordinates from rownames(expr). ",
         "Expected formats include 'x_y', 'x,y', or 'x y'.")
  }

  # ---- pipeline ----
  centres <- generate_kernel_centres_by_density(coords, span = span)

  W0 <- generate_weight_matrix_euclidean(
    coords, centres, span = span, type = kernel_type
  )

  membership <- generate_membership_matrix(W0)
  overlaps   <- determine_overlaps(membership)
  light      <- find_light_kernels(overlaps, membership)
  W          <- trim_weight_matrix(W0, light)

  # guard: no kernels left
  if (ncol(W) == 0L) {
    empty <- list(
      p_values = stats::setNames(numeric(0), character(0)),
      summary  = data.frame(pair = character(0), s2 = numeric(0), z = numeric(0), p = numeric(0)),
      intermediates = list(
        coords = coords, centres = centres, weight_matrix_initial = W0,
        membership = membership, overlaps = overlaps, light_kernels = light,
        weight_matrix = W, local_corr_matrix = NULL, neff = NULL, corr = NULL,
        sigma = NULL, expectation = NA_real_, variance = NA_real_,
        s2 = NULL, z_scores = NULL, p_values = NULL
      ),
      params = list(span = span, kernel_type = kernel_type,
                    fisher_transform = fisher_transform,
                    drop_zero_weight = drop_zero_weight)
    )
    return(empty)
  }

  Z  <- get_local_correlation_matrix(
    expr, W, fisher_transform = fisher_transform, drop_zero_weight = drop_zero_weight
  )

  s2 <- compute_var_of_local_corrs(Z)

  neff <- get_effective_sample_sizes(W)

  # Use unattenuated overlap correlation for Sigma construction
  corr  <- approximate_between_coefficient_correlations_effective(W)
  Sigma <- get_sigma_matrix(corr, neff)

  mu    <- approximate_expectation_effective(Sigma)
  varS2 <- approximate_variance_effective(Sigma)

  z <- compute_test_statistics_standardised_effective(
    variances = varS2, expectations = mu, vars_of_local_corrs = s2
  )
  p <- compute_p_values(z)

  # tidy summary
  summary_df <- data.frame(
    pair = names(s2),
    s2   = as.numeric(s2),
    z    = as.numeric(z),
    p    = as.numeric(p),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    p_values = p,
    summary  = summary_df,
    intermediates = list(
      coords = coords,
      centres = centres,
      weight_matrix_initial = W0,
      membership = membership,
      overlaps = overlaps,
      light_kernels = light,
      weight_matrix = W,
      local_corr_matrix = Z,
      neff = neff,
      corr = corr,
      sigma = Sigma,
      expectation = mu,
      variance = varS2,
      s2 = s2,
      z_scores = z,
      p_values = p
    ),
    params = list(
      span = span,
      kernel_type = kernel_type,
      fisher_transform = fisher_transform,
      drop_zero_weight = drop_zero_weight
    )
  )
}

# ---- internal helper ----
#' @keywords internal
#' @noRd
.parse_coords_from_rownames <- function(rn) {
  # Split on x, _, comma, semicolon, or whitespace
  parts <- strsplit(rn, "[,_;[:space:]x]+", perl = TRUE)
  # keep only rows that split into at least 2 tokens; take first two
  xy <- t(vapply(parts, function(p) {
    if (length(p) < 2) return(c(NA_real_, NA_real_))
    c(as.numeric(p[1]), as.numeric(p[2]))
  }, numeric(2)))
  colnames(xy) <- c("x", "y")
  xy
}
