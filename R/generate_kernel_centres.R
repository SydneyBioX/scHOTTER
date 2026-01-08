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
  dims <- ncol(coords)
  total_area <- prod(bounds[2, ] - bounds[1, ])
  if (total_area <= 0) stop("Degenerate bounds (zero total area). Check coords.")
  density <- n_points / total_area
  approx_kernel_area <- target_n / density
  kernel_side <- sqrt(approx_kernel_area)

  ranges <- lapply(seq_len(dims), function(i) bounds[,i])

  n_per_dim <- c()
  for(i in 1:dims){
    n_per_dim[i] <- max(1L,
                    ceiling((ranges[[i]][2]-ranges[[i]][1])/kernel_side)
    )}

  grids <- list()
  for(i in 1:dims){grids[[i]] <- seq(ranges[[i]][1], ranges[[i]][2],
                                     length.out = n_per_dim[i])}

  centres <- as.matrix(expand.grid(grids))
  storage.mode(centres) <- "double"
  centres
}
