`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

as_dgC <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    x
  } else {
    methods::as(x, "dgCMatrix")
  }
}

resolve_target <- function(target) {
  if (isTRUE(target) || isFALSE(target)) {
    stop("target must be 'mean', 'median', 'auto', or a numeric K", call. = FALSE)
  }
  if (is.numeric(target) && length(target) == 1L && is.finite(target)) {
    return(list(name = "fixed", fixed = as.numeric(target)))
  }
  if (is.character(target) && length(target) == 1L &&
      target %in% c("mean", "median", "auto")) {
    return(list(name = target, fixed = NULL))
  }
  stop("target must be 'mean', 'median', 'auto', or a numeric K", call. = FALSE)
}

matrix_slots_as_csr <- function(mat) {
  mat <- as_dgC(mat)
  list(
    x = as.numeric(mat@x),
    i = as.integer(mat@i),
    p = as.integer(mat@p),
    n_cells = ncol(mat),
    n_features = nrow(mat),
    cell_names = colnames(mat),
    feature_names = rownames(mat)
  )
}

seurat_assay_matrix <- function(object, assay = NULL, layer = "counts") {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      exists("LayerData", envir = asNamespace("SeuratObject"), mode = "function")) {
    SeuratObject::LayerData(object = object[[assay]], layer = layer)
  } else {
    SeuratObject::GetAssayData(object = object, assay = assay, slot = layer)
  }
}

set_seurat_assay_matrix <- function(object, value, assay = NULL, layer = "pflogpf") {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      exists("LayerData<-", envir = asNamespace("SeuratObject"), mode = "function")) {
    SeuratObject::LayerData(object = object[[assay]], layer = layer) <- value
    object
  } else {
    SeuratObject::SetAssayData(object = object, assay = assay, slot = layer, new.data = value)
  }
}

matrix_from_csr_as_features_by_cells <- function(res, n_features, n_cells,
                                                 feature_names = NULL,
                                                 cell_names = NULL) {
  mat <- methods::new(
    "dgCMatrix",
    x = as.numeric(res$data),
    i = as.integer(res$indices),
    p = as.integer(res$indptr),
    Dim = as.integer(c(n_features, n_cells))
  )
  dimnames(mat) <- list(feature_names, cell_names)
  mat
}
