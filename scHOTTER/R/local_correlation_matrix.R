#' Local (kernel-wise) Fisher correlations for all gene pairs
#'
#' For every unordered gene pair in \code{expr}, compute a vector of
#' weighted Pearson correlations \eqn{\hat{r}}—one per kernel—using \code{weight_matrix}.
#' By default the correlations are Fisher-transformed (\eqn{z = \mathrm{atanh}(\hat{r})}).
#' Columns are labeled as \code{"GENE1_GENE2"}.
#'
#' @param expr Numeric matrix (rows = points/cells, cols = genes).
#' @param weight_matrix Point×kernel weight matrix (preferably
#'   \code{Matrix::dgCMatrix}) with \code{nrow(weight_matrix) == nrow(expr)}.
#' @param fisher_transform Logical; if \code{TRUE} (default), return Fisher-\eqn{z}
#'   values; otherwise return raw correlations.
#' @param drop_zero_weight Logical; if \code{TRUE} (default), restrict each kernel's
#'   calculation to rows with strictly positive weight for that kernel.
#'
#' @return A numeric matrix with \code{nrow = ncol(weight_matrix)} (kernels) and
#'   \code{ncol = choose(p, 2)} (gene pairs, where \code{p = ncol(expr)}).
#'   Row names are kernel names; column names are \code{"GENE1_GENE2"}.
#'   Entries are Fisher-\eqn{z} if \code{fisher_transform=TRUE}, else raw \eqn{\hat{r}}.
#'
#' @seealso \code{\link{get_local_correlation_vector_pair}}
#'
#' @examples
#' set.seed(1)
#' expr <- cbind(g1 = rnorm(60), g2 = rnorm(60), g3 = rnorm(60))
#' coords <- cbind(runif(60), runif(60))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
#' Zmat <- get_local_correlation_matrix(expr, W)
#' dim(Zmat)   # kernels x choose(3,2) = kernels x 3
#'
#' @export
get_local_correlation_matrix <- function(expr, weight_matrix,
                                         fisher_transform = TRUE,
                                         drop_zero_weight = TRUE) {
  # Coerce & validate inputs
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr)) stop("expr must be a numeric matrix with gene columns.")
  if (is.data.frame(weight_matrix)) weight_matrix <- as.matrix(weight_matrix)
  if (is.matrix(weight_matrix)) {
    weight_matrix <- Matrix::Matrix(weight_matrix, sparse = TRUE)
    weight_matrix <- methods::as(weight_matrix, "dgCMatrix")
  } else if (!inherits(weight_matrix, "dgCMatrix")) {
    weight_matrix <- methods::as(weight_matrix, "dgCMatrix")
  }
  if (nrow(expr) != nrow(weight_matrix)) {
    stop("nrow(expr) must equal nrow(weight_matrix).")
  }

  gene_names <- colnames(expr)
  if (is.null(gene_names) || length(gene_names) < 2L) {
    stop("expr must have >= 2 named gene columns.")
  }

  # All unordered pairs and their labels
  gene_pairs <- utils::combn(gene_names, 2, simplify = FALSE)
  pair_labels <- vapply(gene_pairs, function(p) paste0(p[1], "_", p[2]), character(1))

  kernel_names <- colnames(weight_matrix)
  if (is.null(kernel_names)) kernel_names <- paste0("k_", seq_len(ncol(weight_matrix)))
  n_kernels <- ncol(weight_matrix)

  out <- matrix(NA_real_, nrow = n_kernels, ncol = length(gene_pairs),
                dimnames = list(kernel_names, pair_labels))

  # inline weighted Pearson (identical to vector version)
  weighted_cor <- function(x, y, w) {
    w_sum <- sum(w, na.rm = TRUE)
    if (!is.finite(w_sum) || w_sum == 0) return(NA_real_)
    mx <- sum(w * x, na.rm = TRUE) / w_sum
    my <- sum(w * y, na.rm = TRUE) / w_sum
    cov_xy <- sum(w * (x - mx) * (y - my), na.rm = TRUE) / w_sum
    var_x  <- sum(w * (x - mx)^2,         na.rm = TRUE) / w_sum
    var_y  <- sum(w * (y - my)^2,         na.rm = TRUE) / w_sum
    denom <- sqrt(var_x * var_y)
    if (!is.finite(denom) || denom == 0) return(NA_real_)
    cov_xy / denom
  }

  # Pre-extract gene columns as a list for a tiny speedup
  X <- expr

  for (i in seq_len(n_kernels)) {
    w <- as.numeric(weight_matrix[, i])

    idx <- if (drop_zero_weight) which(w > 0) else seq_len(length(w))
    if (length(idx) < 2L) next

    for (j in seq_along(gene_pairs)) {
      g1 <- gene_pairs[[j]][1]
      g2 <- gene_pairs[[j]][2]
      r  <- weighted_cor(X[idx, g1], X[idx, g2], w[idx])
      if (isTRUE(fisher_transform)) {
        out[i, j] <- if (is.finite(r) && abs(r) < 1) atanh(r) else NA_real_
      } else {
        out[i, j] <- r
      }
    }
  }

  out
}
