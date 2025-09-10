#' Local (kernel-wise) weighted correlation for a gene pair
#'
#' For a gene pair encoded as \code{"GENE1_GENE2"}, compute a vector of
#' weighted Pearson correlations \eqn{\hat{r}}—one per kernel—from a point×kernel weight
#' matrix. By default, correlations are Fisher-transformed (\eqn{z = \mathrm{atanh}(\hat{r})}).
#'
#' @param expr Numeric matrix (rows = points/cells, cols = genes) with column
#'   names containing the requested genes.
#' @param weight_matrix A point×kernel weight matrix (preferably
#'   \code{Matrix::dgCMatrix}) with the same number of rows as \code{expr}.
#' @param gene_pair Character scalar like \code{"GENE1_GENE2"} indicating the
#'   two gene column names in \code{expr}.
#' @param fisher_transform Logical; if \code{TRUE} (default), apply the Fisher
#'   \eqn{\mathrm{atanh}} transform to each correlation. Values with \eqn{|r| \geq 1}
#'   or non-finite are returned as \code{NA_real_}.
#' @param drop_zero_weight Logical; if \code{TRUE} (default), restrict each
#'   kernel's computation to rows with strictly positive weight for that kernel.
#'   If \code{FALSE}, all rows are considered (zeros contribute nothing).
#'
#' @return A named numeric vector of length \eqn{k} (number of kernels), where
#'   names are the kernel column names in \code{weight_matrix}. Entries are
#'   Fisher-\eqn{z} values if \code{fisher_transform=TRUE}, otherwise raw
#'   correlations. Kernels with fewer than 2 positively weighted points (or
#'   undefined variance) return \code{NA_real_}.
#'
#' @examples
#' set.seed(1)
#' expr <- cbind(g1 = rnorm(60), g2 = rnorm(60))
#' coords <- cbind(runif(60), runif(60))
#' centres <- generate_kernel_centres_by_density(coords, span = 0.2)
#' W <- generate_weight_matrix_euclidean(coords, centres, span = 0.2, type = "gaussian")
#' zvec <- get_local_correlation_vector_pair(expr, W, "g1_g2")
#' head(zvec)
#'
#' @export
get_local_correlation_vector_pair <- function(expr, weight_matrix, gene_pair,
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

  # Parse "GENE1_GENE2"
  parts <- strsplit(gene_pair, "_", fixed = TRUE)[[1]]
  if (length(parts) != 2L) {
    stop("gene_pair must be a single string like 'GENE1_GENE2'.")
  }
  g1 <- parts[1]; g2 <- parts[2]

  if (!(g1 %in% colnames(expr) && g2 %in% colnames(expr))) {
    missing <- setdiff(c(g1, g2), colnames(expr))
    stop(sprintf("Gene(s) not found in expr: %s",
                 paste(missing, collapse = ", ")))
  }

  kernel_names <- colnames(weight_matrix)
  if (is.null(kernel_names)) kernel_names <- paste0("k_", seq_len(ncol(weight_matrix)))
  n_kernels <- ncol(weight_matrix)

  out <- rep(NA_real_, n_kernels)
  names(out) <- kernel_names

  x_all <- expr[, g1]
  y_all <- expr[, g2]

  # inline weighted Pearson (na-tolerant)
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

  for (i in seq_len(n_kernels)) {
    # extract weights for kernel i (as dense numeric)
    w <- as.numeric(weight_matrix[, i])

    idx <- if (drop_zero_weight) which(w > 0) else seq_len(length(w))
    if (length(idx) < 2L) next

    r <- weighted_cor(x_all[idx], y_all[idx], w[idx])

    if (isTRUE(fisher_transform)) {
      out[i] <- if (is.finite(r) && abs(r) < 1) atanh(r) else NA_real_
    } else {
      out[i] <- r
    }
  }

  out
}
