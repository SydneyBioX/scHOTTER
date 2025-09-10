#' Build a point-by-kernel weight matrix (Euclidean neighborhoods)
#'
#' Computes a sparse \eqn{n_{\mathrm{points}} \times n_{\mathrm{kernels}}} weight
#' matrix \code{W} using Euclidean distance from each kernel centre to each point.
#' Two schemes are supported:
#' \itemize{
#'   \item \strong{block}: weights are 1 for neighbors within the search radius, 0 otherwise.
#'   \item \strong{gaussian}: weights are \eqn{\exp\{ -d^2 / (2\sigma^2) \}}, with
#'         \eqn{\sigma = \mathrm{gap}/2}, where \code{gap} is the minimum grid spacing,
#'         and \eqn{d} is the distance from the kernel centre to the point of interest.
#' }
#'
#' @param coords Numeric matrix or data frame with two columns \code{c(x, y)} giving
#'   point coordinates (rows are points).
#' @param grid_centres Numeric matrix or data frame with two columns \code{c(x, y)}
#'   giving kernel-centre coordinates (e.g., from
#'   \code{\link{generate_kernel_centres_by_density}}).
#' @param span Optional numeric in \eqn{(0, 1]}. Target proportion of points to associate
#'   with each kernel (used to choose \code{k = ceiling(span * n_points)} neighbors).
#'   If \code{NULL}, defaults to \eqn{13 / n_{\mathrm{points}}} (capped at 0.5).
#' @param type Character; either \code{"block"} or \code{"gaussian"}.
#'
#' @details
#' \strong{Grid spacing (gap).} The search radius is set to the smallest spacing between
#' unique \code{x}- or \code{y}-coordinates of \code{grid_centres}. If only a single
#' centre exists on an axis, we fall back to a coarse spacing based on the data range.
#'
#' \strong{Neighbor search.} We query up to \eqn{k = \lceil \mathrm{span}\, n\rceil}
#' nearest neighbors per kernel centre using \pkg{RANN}. If your \pkg{RANN} version
#' does not support radius-restricted search, a standard k-NN query still works; the
#' \code{"block"} scheme then effectively selects the returned neighbors.
#'
#' @return A \code{Matrix::dgCMatrix} sparse matrix \code{W} with rownames taken from
#'   \code{coords} (if present) and column names \code{"k_1", ..., "k_m"} for m kernels.
#'
#' @seealso \code{\link{generate_kernel_centres_by_density}}
#'
#' @examples
#' set.seed(1)
#' coords <- cbind(runif(100, 0, 10), runif(100, 0, 5))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.1)
#' Wb <- generate_weight_matrix_euclidean(coords, centres, span = 0.1, type = "block")
#' Wg <- generate_weight_matrix_euclidean(coords, centres, span = 0.1, type = "gaussian")
#' dim(Wb); Matrix::nnzero(Wb)
#'
#' @export
generate_weight_matrix_euclidean <- function(coords, grid_centres,
                                             span = NULL,
                                             type = c("block", "gaussian")) {
  type <- match.arg(type)

  # Coerce & validate inputs
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  if (is.data.frame(grid_centres)) grid_centres <- as.matrix(grid_centres)
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("coords must be a numeric matrix/data.frame with 2 columns (x, y).")
  }
  if (!is.matrix(grid_centres) || ncol(grid_centres) != 2) {
    stop("grid_centres must be a numeric matrix/data.frame with 2 columns (x, y).")
  }
  n_points  <- nrow(coords)
  n_kernels <- nrow(grid_centres)
  if (n_points < 1L || n_kernels < 1L) stop("Need at least 1 point and 1 kernel centre.")

  # Default span
  if (is.null(span)) {
    span <- 13 / n_points
    if (span > 1) span <- 0.5
    message("span not specified, defaulting to ", round(span, 2))
  }
  if (!is.numeric(span) || length(span) != 1 || span <= 0 || span > 1) {
    stop("span must be a single numeric value in (0, 1].")
  }
  target_n <- ceiling(span * n_points)

  # Robust grid spacing (gap): smallest positive spacing in x/y of centres
  ux <- sort(unique(grid_centres[, 1]))
  uy <- sort(unique(grid_centres[, 2]))
  dx <- if (length(ux) >= 2) min(diff(ux)) else diff(range(coords[, 1])) / max(1, sqrt(n_kernels))
  dy <- if (length(uy) >= 2) min(diff(uy)) else diff(range(coords[, 2])) / max(1, sqrt(n_kernels))
  gap <- min(dx, dy)
  if (!is.finite(gap) || gap <= 0) stop("Could not determine a positive grid spacing (gap).")

  # k-NN (optionally radius-restricted) from centres -> points
  # Use RANN explicitly to avoid import ambiguity.
  nn <- RANN::nn2(data = coords, query = grid_centres, k = target_n,
                  searchtype = "radius", radius = gap)

  i <- c(t(nn$nn.idx))
  j <- rep(seq_len(n_kernels), each = target_n)
  distances <- c(t(nn$nn.dists))

  # keep only valid neighbor slots
  keep <- which(i > 0L & is.finite(distances))
  i <- i[keep]
  j <- j[keep]
  distances <- distances[keep]

  # Build sparse weight matrix
  if (type == "block") {
    x <- rep.int(1.0, length(i))
    W <- Matrix::sparseMatrix(i = i, j = j, x = x,
                              dims = c(n_points, n_kernels),
                              dimnames = list(rownames(coords),
                                              paste0("k_", seq_len(n_kernels))))
    # defensive: cap at 1 (should already be 1s)
    if (length(W@x)) W@x[W@x > 1] <- 1
  } else { # gaussian
    sigma <- gap / 2
    x <- exp(-(distances^2) / (2 * sigma^2))
    W <- Matrix::sparseMatrix(i = i, j = j, x = x,
                              dims = c(n_points, n_kernels),
                              dimnames = list(rownames(coords),
                                              paste0("k_", seq_len(n_kernels))))
  }

  W
}
