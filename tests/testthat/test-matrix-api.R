library(scclrR)

test_that("matrix API returns a sparse layer and center vector", {
  counts <- Matrix::Matrix(
    c(
      5, 0, 1, 0,
      0, 3, 0, 2,
      1, 0, 4, 0
    ),
    nrow = 3,
    ncol = 4,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  norm <- normalize_matrix(counts, target = "mean")
  expect_s4_class(norm$sparse, "dgCMatrix")
  expect_equal(dim(norm$sparse), dim(counts))
  expect_length(norm$center, ncol(counts))
  expect_true(is.numeric(norm$k))
})
