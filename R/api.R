#' Estimate count overdispersion for PFlog calibration
#'
#' Estimates the negative-binomial overdispersion parameter used by `target =
#' "auto"`. The input matrix is expected to be features by cells, as in a Seurat
#' assay layer.
#'
#' @param counts Sparse or dense count matrix with features in rows and cells in
#'   columns.
#' @return A list with `alpha`, `mean_depth`, and `k`, where `k = 4 * alpha *
#'   mean_depth`.
#' @export
estimate_overdispersion <- function(counts) {
  s <- matrix_slots_as_csr(counts)
  .scclr_overdispersion(s$x, s$i, s$p, s$n_cells, s$n_features)
}

#' Normalize a count matrix with sparse PFlog
#'
#' Returns the sparse shifted-log values and the per-cell centering vector. The
#' dense normalized matrix is represented implicitly as `sparse - center`.
#'
#' @param counts Count matrix with features in rows and cells in columns.
#' @param target `"auto"`, `"mean"`, `"median"`, or a numeric fixed `K`.
#' @param alpha Optional overdispersion value. If supplied, it overrides
#'   `target` and uses `K = 4 * alpha * mean_depth`.
#' @param log1p Whether to use the shifted logarithm.
#' @param center Whether to return and use the row/cell centering vector.
#' @return A list with `sparse`, `center`, `k`, and `alpha`.
#' @export
normalize_matrix <- function(counts, target = "auto", alpha = NULL,
                             log1p = TRUE, center = TRUE) {
  target <- resolve_target(target)
  s <- matrix_slots_as_csr(counts)
  res <- .scclr_normalize(
    s$x, s$i, s$p, s$n_cells, s$n_features,
    target$name, target$fixed, alpha, isTRUE(log1p), isTRUE(center)
  )
  sparse <- matrix_from_csr_as_features_by_cells(
    res,
    n_features = s$n_features,
    n_cells = s$n_cells,
    feature_names = s$feature_names,
    cell_names = s$cell_names
  )
  list(
    sparse = sparse,
    center = as.numeric(res$row_center),
    k = as.numeric(res$k),
    alpha = if (is.null(res$alpha)) NULL else as.numeric(res$alpha)
  )
}

#' Run sparse PCA on a PFlog matrix representation
#'
#' @param sparse Sparse shifted-log matrix with features in rows and cells in
#'   columns.
#' @param center Optional per-cell centering vector. If supplied, PCA is run on
#'   `t(sparse) - center` without densifying.
#' @param n.components Number of principal components.
#' @param ncv,maxiter,seed,tol Parameters forwarded to the Rust eigensolver.
#' @return A list containing `scores`, `loadings`, explained variance fields,
#'   singular values, and PCA dimensions.
#' @export
pca_matrix <- function(sparse, center = NULL, n.components = 50,
                       ncv = NULL, maxiter = NULL, seed = 0, tol = 0) {
  s <- matrix_slots_as_csr(sparse)
  res <- .scclr_pca(
    s$x, s$i, s$p, s$n_cells, s$n_features,
    center %||% NULL, n.components, ncv, maxiter, seed, tol
  )
  scores <- matrix(
    as.numeric(res$scores),
    nrow = as.integer(res$n_samples),
    ncol = as.integer(res$n_components),
    byrow = TRUE
  )
  components <- matrix(
    as.numeric(res$components),
    nrow = as.integer(res$n_components),
    ncol = as.integer(res$n_features),
    byrow = TRUE
  )
  rownames(scores) <- s$cell_names
  colnames(scores) <- paste0("PC_", seq_len(ncol(scores)))
  loadings <- t(components)
  rownames(loadings) <- s$feature_names
  colnames(loadings) <- colnames(scores)
  list(
    scores = scores,
    loadings = loadings,
    mean = as.numeric(res$mean),
    explained_variance = as.numeric(res$explained_variance),
    explained_variance_ratio = as.numeric(res$explained_variance_ratio),
    singular_values = as.numeric(res$singular_values),
    noise_variance = as.numeric(res$noise_variance),
    n_samples = as.integer(res$n_samples),
    n_features = as.integer(res$n_features),
    n_components = as.integer(res$n_components)
  )
}
