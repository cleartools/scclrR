# These wrappers match the exported extendr symbols from src/rust/src/lib.rs.
# They are intentionally thin; user-facing validation lives in R/api.R.

.scclr_overdispersion <- function(x, i, p, n_cells, n_features) {
  .Call(
    "wrap__scclr_overdispersion",
    x, i, p, as.integer(n_cells), as.integer(n_features),
    PACKAGE = "scclrR"
  )
}

.scclr_normalize <- function(x, i, p, n_cells, n_features, target, fixed, alpha,
                             log1p, center) {
  .Call(
    "wrap__scclr_normalize",
    x, i, p, as.integer(n_cells), as.integer(n_features),
    target, fixed, alpha, log1p, center,
    PACKAGE = "scclrR"
  )
}

.scclr_pca <- function(x, i, p, n_cells, n_features, row_center,
                       n_components, ncv, maxiter, seed, tol) {
  .Call(
    "wrap__scclr_pca",
    x, i, p, as.integer(n_cells), as.integer(n_features), row_center,
    as.integer(n_components), ncv, maxiter, as.double(seed), as.double(tol),
    PACKAGE = "scclrR"
  )
}

.scclr_normalize_pca <- function(x, i, p, n_cells, n_features, n_components,
                                 target, fixed, alpha, ncv, maxiter, seed, tol) {
  .Call(
    "wrap__scclr_normalize_pca",
    x, i, p, as.integer(n_cells), as.integer(n_features), as.integer(n_components),
    target, fixed, alpha, ncv, maxiter, as.double(seed), as.double(tol),
    PACKAGE = "scclrR"
  )
}
