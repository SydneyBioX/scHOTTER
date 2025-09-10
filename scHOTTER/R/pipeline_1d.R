#' End-to-end pipeline (1D trajectory): local co-expression p-values
#'
#' A 1D analogue of \code{scHOTTER_pipeline()} for trajectory-like data where
#' \emph{exactly one} coordinate axis varies across points. Row names of
#' \code{expr} must encode coordinates (e.g. "xxy", "x_y", "x,y" or "x y"). The varying
#' axis is used to place kernel centres and compute 1D weights; all subsequent
#' steps (local Fisher z, analytic null, Z, p) are identical to the 2D pipeline.
#'
#' @param expr Numeric matrix (\eqn{\text{points} \times  \text{genes}}). Row names must encode coordinates.
#'   Column names are gene symbols.
#' @param span Optional numeric in (0, 1]; target proportion of points per kernel.
#'   If \code{NULL}, a heuristic \eqn{13/n} capped at 0.5 is used.
#' @param kernel_type \code{"block"} or \code{"gaussian"}; passed to
#'   \code{\link{generate_weight_matrix_euclidean_1d}}.
#' @param fisher_transform Logical; apply Fisher \code{atanh} to local correlations
#'   (default \code{TRUE}). Keep \code{TRUE} to match the null variance model.
#' @param drop_zero_weight Logical; restrict each kernel’s computation to rows
#'   with weight \code{> 0} for that kernel (default \code{TRUE}).
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
#' t  <- sort(runif(n, 0, 10))
#' y0 <- rep(0, n)                           # constant second axis
#' coords <- cbind(t = t, y = y0)
#' rn <- apply(coords, 1, function(v) paste0(v[1], "_", v[2]))
#' expr <- cbind(g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
#' rownames(expr) <- rn
#'
#' out1d <- scHOTTER_pipeline_1d(expr, span = 0.2, kernel_type = "gaussian")
#' head(out1d$summary)
#'
#' @export
scHOTTER_pipeline_1d <- function(expr,
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
    stop("expr must have rownames encoding coordinates like 'x_y', 'x,y', or 'x y'.")

  coords <- .parse_coords_from_rownames(rownames(expr))
  if (!is.matrix(coords) || ncol(coords) < 1L || any(!is.finite(coords))) {
    stop("Failed to parse numeric coordinates from rownames(expr).")
  }

  # Check data are effectively 1D: exactly one varying axis
  rng <- apply(coords, 2, function(x) range(x, na.rm = TRUE))
  widths <- rng[2, ] - rng[1, ]
  eff <- widths > 0
  dim_effective <- sum(eff)
  if (dim_effective != 1L) {
    stop("Data are not effectively 1D: expected exactly 1 varying coordinate; found ", dim_effective, ". ",
         "Use scHOTTER_pipeline() for 2D/spatial data.")
  }

  # ---- pipeline (1D centres & weights) ----
  centres <- generate_kernel_centres_by_density_1d(coords, span = span)

  W0 <- generate_weight_matrix_euclidean_1d(
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
                    drop_zero_weight = drop_zero_weight,
                    mode = "trajectory_1d")
    )
    return(empty)
  }

  Z  <- get_local_correlation_matrix(
    expr, W, fisher_transform = fisher_transform, drop_zero_weight = drop_zero_weight
  )

  s2 <- compute_var_of_local_corrs(Z)
  neff <- get_effective_sample_sizes(W)

  # Use unattenuated overlap correlation; Sigma applies the 1/(neff-3) scaling
  corr  <- approximate_between_coefficient_correlations_effective(W)
  Sigma <- get_sigma_matrix(corr, neff)

  mu    <- approximate_expectation_effective(Sigma)
  varS2 <- approximate_variance_effective(Sigma)

  z <- compute_test_statistics_standardised_effective(
    variances = varS2, expectations = mu, vars_of_local_corrs = s2
  )
  p <- compute_p_values(z)

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
      drop_zero_weight = drop_zero_weight,
      mode = "trajectory_1d"
    )
  )
}
