# scclrR

`scclrR` is the R/Seurat companion to
[`cleartools/scclr`](https://github.com/cleartools/scclr). It computes sparse
PFlogPF, also called shifted centered log-ratio normalization, and runs sparse
PCA through Rust bindings to `runorm` and `rupca`.

Rendered documentation is published at <https://cleartools.github.io/scclrR/>.
The Angelidis pseudobulk vignette is rendered at
<https://cleartools.github.io/scclrR/articles/seurat-pflogpf.html>.

The key representation is:

- a sparse matrix containing the shifted log values, stored as a Seurat assay
  layer, and
- a per-cell centering vector, stored in Seurat metadata and package metadata.

The dense PFlogPF value is never materialized unless the user asks for it:

```text
PFlogPF[i, j] = sparse_log_values[i, j] - center[i]
```

This is the same sparse-plus-centering representation used by `scclr` for
AnnData/Scanpy.

## Install

You need R, a Rust toolchain with `cargo`, and the usual Seurat dependencies.

```r
remotes::install_github("cleartools/scclrR")
```

For local development:

```r
devtools::document()
devtools::load_all()
```

## Seurat use

```r
library(scclrR)

pbmc <- pflogpf(pbmc, assay = "RNA", layer = "counts", target = "auto")
pbmc <- run_pca(pbmc, assay = "RNA", layer = "pflogpf", n.components = 50)

Embeddings(pbmc, "scclr_pca")[1:5, 1:5]
```

`target = "auto"` estimates the overdispersion `alpha` and sets
`K = 4 * alpha * mean_depth`. Numeric `target` values are interpreted as a
fixed `K`.

See the
[rendered Angelidis pseudobulk vignette](https://cleartools.github.io/scclrR/articles/seurat-pflogpf.html)
for a complete Seurat-style workflow. The source is
`vignettes/seurat-pflogpf.Rmd`.
