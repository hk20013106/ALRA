#' Run ALRA on a Seurat object
#'
#' Applies \code{alra()} to the selected assay data of a Seurat object and stores the
#' imputed matrix back into the object. The helper adapts to the Seurat v5 layer
#' API while remaining backward compatible with the slot-based interface that was
#' used in earlier releases.
#'
#' @param object A Seurat object.
#' @param assay Name of the assay to modify. Defaults to \code{SeuratObject::DefaultAssay(object)}.
#' @param input_layer Name of the layer (or slot in pre-v5 Seurat) to read. Defaults to \code{"data"}.
#' @param output_layer Name of the layer/slot where the imputed matrix should be stored.
#'   The default (\code{NULL}) overwrites the same layer that was read.
#' @param make_sparse Logical flag indicating whether the imputed matrix should be
#'   converted to a sparse \code{dgCMatrix} before being stored. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to \code{alra()}.
#'
#' @return The updated Seurat object.
#' @export
#' @examples
#' \dontrun{
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small <- NormalizeData(pbmc_small)
#' pbmc_small <- alraSeurat(pbmc_small)
#' }
alraSeurat <- function(object,
                       assay = NULL,
                       input_layer = "data",
                       output_layer = NULL,
                       make_sparse = TRUE,
                       ...) {
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject package is required but not installed")
  }

  assay <- assay %||% SeuratObject::DefaultAssay(object)
  if (!assay %in% names(object@assays)) {
    stop(sprintf("Assay '%s' not found in the provided Seurat object", assay))
  }

  reader <- .get_seurat_reader()
  writer <- .get_seurat_writer()

  data_layer <- input_layer %||% "data"
  target_layer <- output_layer %||% data_layer

  expression_matrix <- reader(object = object, assay = assay, slot_or_layer = data_layer)

  cell_names <- colnames(expression_matrix)
  gene_names <- rownames(expression_matrix)

  completion <- alra(t(as.matrix(expression_matrix)), ...)[[3]]
  alra_matrix <- t(completion)

  if (make_sparse) {
    alra_matrix <- Matrix::Matrix(alra_matrix, sparse = TRUE)
  }

  rownames(alra_matrix) <- gene_names
  colnames(alra_matrix) <- cell_names

  object <- writer(object = object, assay = assay, slot_or_layer = target_layer, new_data = alra_matrix)
  object
}

._call_with_layer <- function(fun, args_layer, args_slot) {
  tryCatch(
    do.call(fun, args_layer),
    error = function(err) {
      if (grepl('layer', conditionMessage(err), fixed = TRUE)) {
        do.call(fun, args_slot)
      } else {
        stop(err)
      }
    }
  )
}

.get_seurat_reader <- function() {
  get_fun <- SeuratObject::GetAssayData
  function(object, assay, slot_or_layer) {
    args_layer <- list(object = object, assay = assay, layer = slot_or_layer)
    args_slot <- list(object = object, assay = assay, slot = slot_or_layer)
    ._call_with_layer(get_fun, args_layer, args_slot)
  }
}

.get_seurat_writer <- function() {
  set_fun <- SeuratObject::SetAssayData
  function(object, assay, slot_or_layer, new_data) {
    args_layer <- list(object = object, assay = assay, layer = slot_or_layer, new.data = new_data)
    args_slot <- list(object = object, assay = assay, slot = slot_or_layer, new.data = new_data)
    if ('replace' %in% names(formals(set_fun))) {
      args_layer$replace <- TRUE
      args_slot$replace <- TRUE
    }
    ._call_with_layer(set_fun, args_layer, args_slot)
  }
}

`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) {
    lhs
  } else {
    rhs
  }
}
