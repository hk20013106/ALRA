# Example script for README

# 1. Load necessary packages
if (!requireNamespace("Seurat", quietly = TRUE)) {
    install.packages("Seurat")
}
if (!requireNamespace("SeuratData", quietly = TRUE)) {
    install.packages("SeuratData")
}
if (!requireNamespace("ALRA", quietly = TRUE)) {
    remotes::install_github("hk20013106/ALRA")
}

library(Seurat)
library(SeuratData)
library(ALRA)

# 2. Load and prepare pbmc3k dataset
if (!"pbmc3k" %in% installed.packages()) {
    InstallData("pbmc3k")
}
data("pbmc3k")

# The pbmc3k dataset from SeuratData might be an older Seurat object.
# Update it to the latest structure.
pbmc3k <- UpdateSeuratObject(pbmc3k)


# 3. Perform ALRA imputation
# ALRA works on the log-normalized data.
pbmc3k <- NormalizeData(pbmc3k, verbose = FALSE)

# Run ALRA on the normalized data matrix
A_norm <- as.matrix(GetAssayData(pbmc3k, slot = "data"))
alra_results <- alra(A_norm)

# Create a new assay with the imputed data
imputed_matrix <- alra_results[[3]]
rownames(imputed_matrix) <- rownames(A_norm)
colnames(imputed_matrix) <- colnames(A_norm)
pbmc3k[["alra"]] <- CreateAssayObject(counts = imputed_matrix)
pbmc3k_alra <- pbmc3k


# 4. Test the validity of the analysis
# A simple test is to check if a known marker gene shows higher expression
# in the relevant cell type after imputation.
# For example, 'MS4A1' is a marker for B cells.

# Get the original and imputed expression of MS4A1
original_ms4a1 <- GetAssayData(pbmc3k, assay = "RNA", slot = "data")["MS4A1", ]
imputed_ms4a1 <- GetAssayData(pbmc3k_alra, assay = "alra", slot = "data")["MS4A1", ]

# Check the number of cells expressing MS4A1 before and after
cells_expressing_before <- sum(original_ms4a1 > 0)
cells_expressing_after <- sum(imputed_ms4a1 > 0)

cat("Cells expressing MS4A1 before ALRA:", cells_expressing_before, "\n")
cat("Cells expressing MS4A1 after ALRA:", cells_expressing_after, "\n")

# A successful imputation should increase the number of expressing cells
if (cells_expressing_after > cells_expressing_before) {
    cat("Validation successful: ALRA increased the detection of MS4A1.\n")
} else {
    cat("Validation failed or no change: ALRA did not increase detection of MS4A1.\n")
}

# Another test: check the average expression in B cells
b_cell_ident <- "B" # This might need adjustment based on pbmc3k annotations
if (b_cell_ident %in% levels(pbmc3k)) {
    b_cells <- WhichCells(pbmc3k, idents = b_cell_ident)
    avg_expr_before <- mean(original_ms4a1[b_cells])
    avg_expr_after <- mean(imputed_ms4a1[b_cells])
    cat("Average MS4A1 in B cells before:", avg_expr_before, "\n")
    cat("Average MS4A1 in B cells after:", avg_expr_after, "\n")
    if (avg_expr_after > avg_expr_before) {
        cat("Validation successful: ALRA increased average expression of MS4A1 in B cells.\n")
    }
} else {
    cat("Could not find B cell identity to run further validation.\n")
}
