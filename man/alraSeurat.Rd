\name{alraSeurat}
\alias{alraSeurat}
\title{Run ALRA on a Seurat object}
\usage{
alraSeurat(object, assay = NULL, input_layer = "data", output_layer = NULL,
  make_sparse = TRUE, ...)
}
\arguments{
\item{object}{A Seurat object.}

\item{assay}{Name of the assay to modify. Defaults to \code{SeuratObject::DefaultAssay(object)}.}

\item{input_layer}{Name of the layer (or slot in legacy Seurat) to read. Defaults to \code{"data"}.}

\item{output_layer}{Name of the layer or slot where the imputed matrix should be stored. The default (\code{NULL}) overwrites the same layer that was read.}

\item{make_sparse}{Logical flag indicating whether the imputed matrix should be converted to a sparse \code{dgCMatrix} before being stored.}

\item{...}{Additional arguments passed to \code{alra()}.}
}
\value{
An updated Seurat object with ALRA-completed data stored in the selected layer.
}
\description{
Applies \code{alra()} to the selected assay data of a Seurat object and stores the imputed matrix using the Seurat v5 layer API while remaining backward compatible with the slot-based interface from earlier releases.
}
\examples{
\dontrun{
library(Seurat)
data("pbmc_small")
pbmc_small <- NormalizeData(pbmc_small)
pbmc_small <- alraSeurat(pbmc_small)
}
}
