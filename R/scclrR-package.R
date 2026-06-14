#' scclrR: Sparse PFlogPF normalization and PCA for Seurat
#'
#' `scclrR` stores PFlogPF as a sparse shifted-log matrix plus a per-cell
#' centering vector, and runs PCA on the implicit centered matrix through Rust.
#'
#' @keywords internal
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods as new
#' @useDynLib scclrR
"_PACKAGE"
