#' Generate kernel centres for 1D (trajectory) coordinates
#'
#' For data where exactly one coordinate dimension varies across points
#' (i.e. an effectively 1D trajectory embedded in \eqn{\mathbb{R}^p}),
#' place kernel centres along the varying axis so that, on average, each kernel
#' collects about \code{span * n_points} points (under a uniform-density
#' approximation). Non-varying dimensions are fixed at their column medians.
#'
#' @param coords Numeric matrix or data frame of coordinates with \code{n} rows
#'   (points) and \code{p} columns (axes). Exactly one column must have
#'   non-zero range; other columns may be constant (or nearly so). The nonzero
#'   coordinate may correspond to a pseudotimepoint or similar.
#' @param span Optional numeric in \eqn{(0,1]}. Target proportion of points per
#'   kernel. If \code{NULL}, a heuristic default \eqn{13 / n} is used (capped at
#'   \code{0.5}); a message is emitted.
#'
#' @return A numeric matrix with the same columns as \code{coords}. The varying
#'   axis is filled with an evenly spaced grid over its range; constant axes are
#'   filled with their (NA-robust) medians. Row count equals the number of
#'   kernel centres.
#'
#' @details
#' Let \eqn{L} be the range length of the varying axis and \eqn{n} the number of
#' points. We approximate the 1D density by \eqn{n / L}, target the number of
#' points per kernel as \eqn{\lceil \mathrm{span} \cdot n \rceil}, set the kernel
#' length to \eqn{\mathrm{target\_n} / (n/L)}, and place
#' \eqn{\lceil L / \mathrm{kernel\_len} \rceil} centres evenly over the axis.
#'
#' @examples
#' set.seed(1)
#' t  <- sort(runif(100, 0, 10))   # pseudotime
#' y0 <- rep(0, 100)               # constant second dim
#' centres <- generate_kernel_centres_by_density_1d(cbind(t, y0), span = 0.2)
#' head(centres)
#'
#' @export
generate_kernel_centres_by_density_1d <- function(coords, span = NULL){
  # Coerce & validate
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  if (!is.matrix(coords) || ncol(coords) < 1L)
    stop("coords must be a numeric matrix/data.frame with >=1 column.")
  n_points <- nrow(coords)
  if (n_points < 2L) stop("Need at least 2 points.")

  # Default span (heuristic)
  if (is.null(span)) {
    span <- 13 / n_points
    if (span > 1) span <- 0.5
    message("span not specified, defaulting to ", round(span, 2))
  }
  if (!is.numeric(span) || length(span) != 1L || span <= 0 || span > 1)
    stop("span must be a single numeric value in (0, 1].")

  target_n <- max(1L, ceiling(span * n_points))

  # Effective dimension: columns with non-zero range
  rng <- apply(coords, 2, function(x) range(x, na.rm = TRUE))
  widths <- rng[2, ] - rng[1, ]
  eff <- widths > 0
  dim_effective <- sum(eff)

  if (dim_effective != 1L) {
    stop("Data is not effectively 1D: expected exactly 1 varying axis; found ", dim_effective, ".")
  }

  # Identify the varying axis and its length
  x_idx <- which(eff)[1]
  x_rng <- rng[, x_idx]
  L <- diff(x_rng)
  if (!is.finite(L) || L <= 0) stop("Degenerate range on the varying axis.")

  # 1D density & number of kernel positions
  density_1d <- n_points / L
  kernel_len <- target_n / density_1d
  n_x <- max(1L, ceiling(L / kernel_len))

  # Even grid over the varying axis
  grid_x <- if (n_x <= 1L) rep(x_rng[1], 1L) else seq(x_rng[1], x_rng[2], length.out = n_x)

  # Fixed values (medians) for non-varying axes
  fixed_vals <- if (any(!eff)) {
    apply(coords[, !eff, drop = FALSE], 2, function(z) stats::median(z, na.rm = TRUE))
  } else {
    numeric(0)
  }

  # Build centres matrix
  centres <- matrix(NA_real_, nrow = length(grid_x), ncol = ncol(coords))
  colnames(centres) <- colnames(coords)
  centres[, x_idx] <- grid_x
  if (any(!eff)) {
    for (j in which(!eff)) {
      # name-safe lookup into fixed_vals
      val <- if (!is.null(names(fixed_vals))) fixed_vals[colnames(coords)[j]] else fixed_vals
      centres[, j] <- rep(val, length(grid_x))
    }
  }

  storage.mode(centres) <- "double"
  centres
}
