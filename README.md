# Adaptively-thresholded Low Rank Approximation (ALRA)

ALRA provides zero-preserving low-rank completion for single-cell RNA-seq count matrices. The package exposes the core `alra()` pipeline, utilities for rank selection and normalisation, and seamless integration with modern Seurat workflows.

## Key Features
- **Zero-preserving matrix completion** using randomized SVD and adaptive thresholding.
- **Automatic rank heuristics** via `choose_k()` to pick an informative latent dimension.
- **Convenience normalisation** with `normalize_data()` for library-size scaling and log conversion.
- **Optional MKL acceleration** by enabling `use.mkl = TRUE` when the [`fastRPCA`](https://github.com/KlugerLab/rpca-mkl) backend is installed.

## Installation

ALRA targets R \>= 4.1 and requires the `Matrix`, `rsvd`, and `SeuratObject` packages. To install the development version from this repository:

```r
# install.packages("remotes")  # if remotes is missing
remotes::install_github("hk20013106/ALRA", dependencies = TRUE)
```

For Windows users leveraging MKL acceleration, install `fastRPCA` before calling `alra()` with `use.mkl = TRUE`:

```r
remotes::install_github("KlugerLab/rpca-mkl/fastRPCA")
```

After installation, load the package as usual:

```r
library(ALRA)
```

## Quick Start: Core Workflow

1. **Prepare a genes-by-cells matrix.** ALRA expects rows = genes, columns = cells.
2. **Normalise the matrix.** Use `Seurat::NormalizeData` or a custom method.
3. **Select rank `k`.** Call `choose_k()` for an automated heuristic or pass a fixed value.
4. **Run ALRA.** Execute `alra()` and extract the completed matrix (`[[3]]`).

```r
library(ALRA)
library(Seurat)

# Assume pbmc is a Seurat object with raw counts
pbmc <- NormalizeData(pbmc, verbose = FALSE)
A_norm <- as.matrix(GetAssayData(pbmc, slot = "data"))

k_choice <- choose_k(A_norm)
result <- alra(A_norm, k = k_choice$k)
completed_matrix <- result[[3]]
```

`completed_matrix` preserves the original dimensionality while imputing low-abundance transcripts in a zero-aware manner.

## Seurat Integration Example with pbmc3k

This example demonstrates a full workflow: loading the `pbmc3k` dataset, running a standard Seurat analysis, performing ALRA imputation, and visualizing the results side-by-side.

The complete, runnable code for this analysis is available in the `full_analysis_example.R` script.

### Visual Comparison of Original vs. ALRA Imputed Data

After running the analysis, we can generate `FeaturePlots` to compare the expression of key marker genes on the UMAP projections derived from both the original and the ALRA-imputed data. The figure below shows that ALRA effectively enhances the expression signal of marker genes, making cell type-specific patterns more distinct.

![Comparison of Original and ALRA Imputed Data](alra_marker_comparison.png)

*Figure: Feature plots of canonical marker genes before (top row, RNA assay) and after (bottom row, ALRA assay) imputation. The UMAP for each set was calculated independently. ALRA clarifies the expression patterns for markers like MS4A1 (B cells) and GNLY (NK cells), which are sparse in the original data.*


## Testing and Reproducibility

Developers can run the included unit tests (requires `pkgload`, `testthat`, `Seurat`, and `SeuratObject`):

```r
pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)
testthat::test_dir("tests/testthat", reporter = "summary")
```

Set a random seed before calling `alra()` to mirror published results, as the randomized SVD initialisation introduces stochasticity:

```r
set.seed(42)
```

## Troubleshooting

- **`fastRPCA` missing**: install via GitHub or run with `use.mkl = FALSE` (default) to use the pure R backend.
- **Layer vs. slot errors**: Ensure your Seurat objects are updated to the v5 structure using `UpdateSeuratObject()`.

## Citation

If you use ALRA in your research, please cite the accompanying preprint: *Zero-preserving imputation of scRNA-seq data using low-rank approximation* (Linderman et al., bioRxiv 2018).

## License

ALRA is distributed under the MIT license. See `LICENSE.txt` for details.
