#' Add sparse PFlog normalization to a Seurat object
#'
#' `pflog()` reads a Seurat assay layer, computes PFlog using Rust, stores
#' the sparse shifted-log values as a new assay layer, and stores the per-cell
#' centering vector in `object[[]]`. Downstream PCA can then use the implicit
#' matrix `layer - center` without densifying.
#'
#' PFlog is the centered log-ratio of the counts shifted by a uniform pseudocount
#' `1/(4 * alpha)`, `center(log(x + 1/(4 * alpha)))` (computed as the equivalent
#' sparsity-preserving `center(log1p(4 * alpha * x))`).
#'
#' @param object A Seurat object.
#' @param assay Assay name. Defaults to `DefaultAssay(object)`.
#' @param layer Input layer, usually `"counts"`.
#' @param target `"auto"` to estimate alpha and set `K = 4 * alpha *
#'   mean_depth`, `"mean"`, `"median"`, or a numeric fixed `K`.
#' @param alpha Optional overdispersion value, overriding `target`.
#' @param key.added Name of the assay layer that will hold sparse shifted-log
#'   values.
#' @param center.key Metadata column for the per-cell centering vector.
#' @param log1p Whether to use shifted log.
#' @return The modified Seurat object.
#' @export
pflog <- function(object, assay = NULL, layer = "counts", target = "auto",
                  alpha = NULL, key.added = "pflog",
                  center.key = paste0(key.added, "_center"),
                  log1p = TRUE) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("The SeuratObject package is required for Seurat integration.", call. = FALSE)
  }
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  counts <- seurat_assay_matrix(object, assay = assay, layer = layer)
  norm <- normalize_matrix(counts, target = target, alpha = alpha, log1p = log1p, center = TRUE)

  object@misc$scclrR <- object@misc$scclrR %||% list()
  object <- set_seurat_assay_matrix(object, norm$sparse, assay = assay, layer = key.added)
  object[[center.key]] <- norm$center
  object@misc$scclrR[[key.added]] <- list(
    assay = assay,
    input_layer = layer,
    layer = key.added,
    center_key = center.key,
    k = norm$k,
    alpha = norm$alpha,
    target = target,
    log1p = isTRUE(log1p),
    representation = "sparse shifted-log layer plus per-cell center"
  )
  object
}

#' Run sparse PFlog PCA on a Seurat object
#'
#' Reads the sparse PFlog layer and centering vector written by `pflog()`,
#' runs PCA via Rust on the implicit centered matrix, and stores a Seurat
#' dimensional reduction.
#'
#' @param object A Seurat object.
#' @param assay Assay name. Defaults to `DefaultAssay(object)`.
#' @param layer PFlog sparse layer.
#' @param center.key Metadata column containing the centering vector.
#' @param reduction.name Name of the Seurat dimensional reduction to write.
#' @param reduction.key Prefix for component names.
#' @param n.components Number of PCs.
#' @param ncv,maxiter,seed,tol Parameters forwarded to the Rust eigensolver.
#' @return The modified Seurat object.
#' @export
run_pca <- function(object, assay = NULL, layer = "pflog",
                    center.key = paste0(layer, "_center"),
                    reduction.name = "scclr_pca",
                    reduction.key = "sclr_",
                    n.components = 50, ncv = NULL, maxiter = NULL,
                    seed = 0, tol = 0) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("The SeuratObject package is required for Seurat integration.", call. = FALSE)
  }
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  sparse <- seurat_assay_matrix(object, assay = assay, layer = layer)
  if (!center.key %in% colnames(object[[]])) {
    stop("Could not find centering vector metadata column: ", center.key, call. = FALSE)
  }
  pca <- pca_matrix(
    sparse,
    center = object[[center.key, drop = TRUE]],
    n.components = n.components,
    ncv = ncv,
    maxiter = maxiter,
    seed = seed,
    tol = tol
  )
  object@misc$scclrR <- object@misc$scclrR %||% list()
  dr <- SeuratObject::CreateDimReducObject(
    embeddings = pca$scores,
    loadings = pca$loadings,
    key = reduction.key,
    assay = assay
  )
  object[[reduction.name]] <- dr
  object@misc$scclrR[[reduction.name]] <- list(
    assay = assay,
    layer = layer,
    center_key = center.key,
    n_components = n.components,
    explained_variance = pca$explained_variance,
    explained_variance_ratio = pca$explained_variance_ratio,
    singular_values = pca$singular_values
  )
  object
}
