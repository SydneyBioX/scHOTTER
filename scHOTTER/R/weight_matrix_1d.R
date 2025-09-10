#' Build a point-by-kernel weight matrix for 1D (trajectory) data
#'
#' Computes a sparse \eqn{n_{\mathrm{points}} \times n_{\mathrm{kernels}}} weight matrix
#' using Euclidean distance along the single varying coordinate (pseudotime).
#' Two schemes are supported:
#' \itemize{
#'   \item \strong{block}: weights are 1 for neighbors within the search radius, 0 otherwise.
#'   \item \strong{gaussian}: weights are \eqn{\exp\{-d^2 / (2\sigma^2)\}}, with
#'         \eqn{\sigma = \mathrm{gap}/2}, where \code{gap} is the kernel spacing (or a
#'         fallback radius computed from density) and \eqn{d} is the distance from the
#'         kernel centre to the point of interest.
#' }
#'
#' @param coords Numeric matrix/data frame of coordinates. Exactly one column must
#'   vary (effective 1D). Other columns may be constant and are ignored for distances.
#' @param grid_centres Numeric matrix/data frame of kernel-centre coordinates
#'   (e.g., from \code{\link{generate_kernel_centres_by_density_1d}}).
#' @param span Optional numeric in \eqn{(0, 1]}. Target proportion of points to associate
#'   with each kernel (\code{k = ceiling(span * n_points)}). If \code{NULL}, defaults to
#'   \eqn{13 / n_{\mathrm{points}}} (capped at 0.5).
#' @param type Character; either \code{"block"} or \code{"gaussian"}.
#'
#' @return A \code{Matrix::dgCMatrix} with rows = points and columns = kernels
#'   (\code{"k_1", ..., "k_m"}).
#'
#' @seealso \code{\link{generate_kernel_centres_by_density_1d}},
#'   \code{\link{generate_weight_matrix_euclidean}}
#'
#' @examples
#' set.seed(1)
#' t  <- sort(runif(80, 0, 10))
#' y0 <- rep(0, 80)
#' coords <- cbind(t = t, y = y0)
#' centres <- generate_kernel_centres_by_density_1d(coords, span = 0.2)
#' Wb <- generate_weight_matrix_euclidean_1d(coords, centres, span = 0.2, type = "block")
#' Wg <- generate_weight_matrix_euclidean_1d(coords, centres, span = 0.2, type = "gaussian")
#' dim(Wb); Matrix::nnzero(Wb)
#'
#' @export
generate_weight_matrix_euclidean_1d <- function(coords, grid_centres,
                                                span = NULL,
                                                type = c("block", "gaussian")){
  type <- match.arg(type)

  # ---- coerce & validate ----
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  if (is.data.frame(grid_centres)) grid_centres <- as.matrix(grid_centres)
  if (!is.matrix(coords) || !is.matrix(grid_centres))
    stop("coords and grid_centres must be numeric matrices/data.frames.")
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

  # Effective 1D check (exactly one varying axis)
  rng <- apply(coords, 2, function(x) range(x, na.rm = TRUE))
  widths <- rng[2, ] - rng[1, ]
  eff <- widths > 0
  dim_effective <- sum(eff)
  if (dim_effective != 1L) {
    stop("Data is not effectively 1D: expected exactly 1 varying axis; found ", dim_effective, ".")
  }
  var_idx <- which(eff)[1]

  # ---- determine 'gap' (radius / spacing) ----
  min_pos_diff <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 2) return(Inf)
    u <- sort(unique(v))
    if (length(u) < 2) return(Inf)
    d <- diff(u)
    d <- d[d > 0]
    if (length(d) == 0) return(Inf)
    min(d)
  }

  # Prefer spacing between centres along the varying axis
  gap <- min_pos_diff(grid_centres[, var_idx])

  # Fallback: pick radius so about 'target_n' points fall within +/- gap
  if (!is.finite(gap) || gap == 0) {
    L <- widths[var_idx]
    dens <- n_points / L
    # want approx: 2 * gap * dens ~= target_n  -> gap = target_n / (2 * dens)
    gap <- target_n / (2 * dens)
  }
  if (!is.finite(gap) || gap <= 0) stop("Could not determine a positive search radius 'gap'.")

  # ---- 1D neighbor search (use only the varying axis) ----
  data_1d  <- matrix(coords[, var_idx], ncol = 1)
  query_1d <- matrix(grid_centres[, var_idx], ncol = 1)
  nn <- RANN::nn2(data = data_1d, query = query_1d, k = target_n,
                  searchtype = "radius", radius = gap)

  i <- c(t(nn$nn.idx))
  j <- rep(seq_len(n_kernels), each = target_n)
  distances <- c(t(nn$nn.dists))

  # keep valid matches
  keep <- which(i > 0L & is.finite(distances))
  i <- i[keep]; j <- j[keep]; distances <- distances[keep]

  # ---- build sparse weights ----
  if (type == "block") {
    x <- rep.int(1.0, length(i))
    W <- Matrix::sparseMatrix(i = i, j = j, x = x,
                              dims = c(n_points, n_kernels),
                              dimnames = list(rownames(coords),
                                              paste0("k_", seq_len(n_kernels))))
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
