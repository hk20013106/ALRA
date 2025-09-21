# R script to run ALRA on a PBMC toy dataset

# 1. Install and load necessary packages
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
if (!requireNamespace("gdown", quietly = TRUE)) {
    remotes::install_github("wleepang/gdown")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
    install.packages("Seurat")
}
if (!requireNamespace("reticulate", quietly = TRUE)) {
    install.packages("reticulate")
}
if (!requireNamespace("ALRA", quietly = TRUE)) {
    remotes::install_github("hk20013106/ALRA")
}

library(gdown)
library(Seurat)
library(reticulate)
library(ALRA)

# Check for python dependencies for reading .h5ad
if (!py_module_available("scanpy")) {
    py_install("scanpy")
}
if (!py_module_available("anndata")) {
    py_install("anndata")
}
if (!py_module_available("h5py")) {
    py_install("h5py")
}


# 2. Download the pbmc.h5ad file
file_id <- "1LaYOadbotGC6gXAlo-aKfHz-spoFnawk"
file_name <- "pbmc.h5ad"
if (!file.exists(file_name)) {
    cat("Downloading pbmc.h5ad...\n")
    gdown::download(id = file_id, output = file_name, quiet = FALSE)
    cat("Download complete.\n")
} else {
    cat("pbmc.h5ad already exists.\n")
}

# 3. Load the data and prepare for ALRA
cat("Loading data...\n")
ad <- reticulate::import("anndata", convert = FALSE)
adata <- ad$read_h5ad(file_name)
so <- SeuratObject::as.Seurat(adata, counts = "X", data = NULL)

# From the notebook, the counts are in a layer called "counts"
# Let's use that.
so <- SeuratObject::as.Seurat(adata, counts = "counts", data = NULL)


# We use library and log normalization
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
A_norm <- as.matrix(GetAssayData(so, slot = "data"))

# 4. Run ALRA
cat("Running ALRA...\n")
# Note: ALRA expects genes as rows and cells as columns, which is Seurat's format.
# The user prompt said "cells are rows and genes are columns", so we might need to transpose.
# Let's check the ALRA function documentation. The original implementation expects genes x cells.
# Seurat's GetAssayData returns genes x cells, so we don't need to transpose A_norm.

# Let's choose the k for ALRA automatically
k_choice <- choose_k(A_norm)
A_alra <- alra(A_norm, k = k_choice$k)

# The result is a list with the imputed matrix
A_alra_imputed <- A_alra[[3]]

# 5. Confirm analysis ran
cat("ALRA analysis complete.\n")
cat("Original matrix dimensions:", dim(A_norm), "\n")
cat("Imputed matrix dimensions:", dim(A_alra_imputed), "\n")
cat("A few values from the imputed matrix:\n")
print(A_alra_imputed[1:5, 1:5])
