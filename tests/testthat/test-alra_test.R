library(testthat)

set.seed(42)
synthetic_counts <- matrix(rpois(400, lambda = 2), nrow = 40, ncol = 10)
synthetic_counts[synthetic_counts <= 1] <- 0
rownames(synthetic_counts) <- paste0('cell', seq_len(nrow(synthetic_counts)))
colnames(synthetic_counts) <- paste0('gene', seq_len(ncol(synthetic_counts)))

A <- synthetic_counts

# Use a small k/q to keep the test fast while exercising the core pipeline
A_norm <- normalize_data(A)
alra_completed <- alra(A_norm, k = 5, q = 1)

small_counts <- Matrix::Matrix(t(A), sparse = TRUE)

get_data_layer <- function(obj) {
  getter <- SeuratObject::GetAssayData
  call_args <- list(object = obj, assay = SeuratObject::DefaultAssay(obj), layer = 'data')
  tryCatch(
    do.call(getter, call_args),
    error = function(err) {
      call_args$layer <- NULL
      call_args$slot <- 'data'
      do.call(getter, call_args)
    }
  )
}

test_that('alra returns completed matrices of expected shape', {
  expect_type(alra_completed, 'list')
  expect_equal(length(alra_completed), 3L)
  expect_equal(dim(alra_completed[[1]]), dim(A_norm))
  expect_equal(dim(alra_completed[[3]]), dim(A_norm))
  expect_true(all(alra_completed[[3]] >= 0))
})

test_that('alraSeurat updates Seurat objects', {
  skip_if_not_installed('Seurat')
  skip_if_not_installed('SeuratObject')

  seurat_obj <- Seurat::CreateSeuratObject(counts = small_counts)
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)

  baseline <- get_data_layer(seurat_obj)
  updated <- alraSeurat(seurat_obj, k = 5)
  completed <- get_data_layer(updated)

  expect_s4_class(updated, 'Seurat')
  expect_true(all(dim(completed) == dim(baseline)))
  expect_false(identical(as.matrix(completed), as.matrix(baseline)))
})
