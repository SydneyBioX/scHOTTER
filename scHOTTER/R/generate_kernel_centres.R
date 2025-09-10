#' Generate approximate kernel centres on a grid to hit a target span
#'
#' Places kernel centres on a regular grid over the 2D coordinate bounds so that,
#' on average, each kernel would contain roughly \code{span * n_points} data points
#' (under a uniform-density approximation). This is useful for building spatial
#' kernels whose size scales with the observed point density.
#'
#' @param coords A numeric matrix or data frame with two columns giving x,y
#'   coordinates of points (rows are points). Column order is \code{c(x, y)}.
#' @param span Optional numeric in \eqn{(0, 1]}. Interpreted as the target
#'   proportion of data points per kernel. If \code{NULL}, a heuristic default
#'   \eqn{13 / n_{\mathrm{points}}} is used (capped at \code{0.5}); a message is emitted.
#'
#' @details
#' Let \eqn{n} be the number of points and \eqn{A} the area of the bounding box
#' of \code{coords}. We approximate point density by \eqn{n / A} and choose a
#' square kernel with side length \eqn{\sqrt{(\mathrm{target\_n}) / (n/A)}} where
#' \eqn{\mathrm{target\_n} = \lceil \mathrm{span} \times n \rceil}. Grid spacing
#' equals this side length, producing \eqn{n_x \times n_y} centres that tile the
#' bounding box.
#'
#' @return A numeric matrix with two columns \code{c("x","y")} giving kernel
#'   centre coordinates.
#'
#' @examples
#' set.seed(1)
#' coords <- cbind(runif(100, 0, 10), runif(100, 0, 5))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.1)
#' dim(centres)       # number of centres
#' head(centres)
#'
#' @export
generate_kernel_centres_by_density <- function(coords, span = NULL){
  # basic input checks (lightweight; keep examples fast)
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  if (!is.matrix(coords) || ncol(coords) != 2)
    stop("coords must be a numeric matrix/data.frame with 2 columns (x, y).")
  if (!is.null(span) && (!is.numeric(span) || length(span) != 1 || span <= 0 || span > 1))
    stop("span must be NULL or a single numeric value in (0, 1].")

  n_points <- nrow(coords)
  if (n_points < 2) stop("Need at least 2 points.")

  if (is.null(span)) {
    span <- 13 / n_points
    if (span > 1) span <- 0.5
    message("span not specified, defaulting to ", round(span, 2))
  }

  target_n <- ceiling(span * n_points)

  bounds <- apply(coords, 2, range)
  area <- prod(bounds[2, ] - bounds[1, ])
  if (area <= 0) stop("Degenerate bounds (zero area). Check coords.")
  density <- n_points / area
  approx_kernel_area <- target_n / density
  kernel_side <- sqrt(approx_kernel_area)

  x_range <- bounds[, 1]
  y_range <- bounds[, 2]
  n_x <- max(1L, ceiling((x_range[2] - x_range[1]) / kernel_side))
  n_y <- max(1L, ceiling((y_range[2] - y_range[1]) / kernel_side))

  grid_x <- seq(x_range[1], x_range[2], length.out = n_x)
  grid_y <- seq(y_range[1], y_range[2], length.out = n_y)

  centres <- as.matrix(expand.grid(x = grid_x, y = grid_y))
  storage.mode(centres) <- "double"
  centres
}
