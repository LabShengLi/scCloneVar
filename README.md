<h1>
  scCloneVar
  <img src="Figures/figure_v1.png" align="right" width="180">
</h1>

Clone-wise heterogeneity and differential expression analysis toolkit for single-cell RNA-seq data.

## Functions Overview

scCloneVar is an integrated toolkit for studying clonal structure and transcriptional heterogeneity in single-cell RNA-seq data. It supports clone-level differential expression analysis, enabling comparisons between user-defined clone groups using multiple statistical approaches and gene sets. The package quantifies transcriptional variability within and between clones in PCA space and identifies differential variability genes (DVGs), capturing heterogeneity that may not be detected through mean expression analysis alone. In addition, scCloneVar performs Output Activity (OA) analysis to assess lineage bias across clones and generates publication-ready visualizations and summary reports. Functional interpretation is supported through pathway enrichment analysis using MSigDB-based GSEA.


## Installation

```r
install.packages("devtools")
devtools::install_github("LabShengLi/scCloneVar")
library(scCloneVar)
```

## Required Packages

scCloneVar depends on several widely used single-cell and statistical R packages. These dependencies are typically installed automatically, but manual installation may help resolve installation issues on some systems.

**Core single-cell and data wrangling**

```r
install.packages(c("Seurat", "Matrix", "dplyr", "tidyr", "tibble",
                   "purrr", "magrittr", "stringr", "glue", "scales",
                   "openxlsx", "progress", "progressr"))
```

**Plotting and visualization**

```r
install.packages(c("ggplot2", "patchwork", "ggrepel", "ggsignif",
                   "ggvenn", "RColorBrewer", "ggsci", "colorspace",
                   "ggplotify"))
```

**Pathway analysis (Bioconductor)**

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "enrichplot", "msigdbr",
                       "org.Mm.eg.db", "org.Hs.eg.db"))
```

A few quick notes from running this:

- `org.Mm.eg.db` and `org.Hs.eg.db` are only needed if you actually run the GSEA step — install whichever one matches your species.
- `MAST` is used internally by `Seurat::FindMarkers()` when you select `test.use = "MAST"`. If you only ever use the Wilcoxon test you can skip it, but most of the DEG workflow assumes MAST is available (`BiocManager::install("MAST")`).
- The `msigdbr` API changed somewhere around v10 (`collection` / `subcollection` instead of `category` / `subcategory`).

## Functions

### `run_clonewise_DEG_suite()`

The main differential expression analysis function. Given a Seurat object and two clone groups, it performs DEG testing across three gene sets (all genes, top HVGs, and filtered HVGs) using Wilcoxon and MAST. Results include BH/BY-adjusted p-values, group-level average expression, and counts of significant genes based on user-defined log2FC and FDR thresholds.

Outputs are saved to a timestamped directory as a single Excel workbook with summary tables and comparison-specific results. A labeled volcano plot is generated for the MAST analysis on the top HVG gene set.

### `compute_pca_clone_heterogeneity()`

This function quantifies transcriptional heterogeneity in PCA space. After subsetting the data to selected cell types (default: LT-HSC and ST-HSC), it performs highly variable gene selection, scaling, and PCA on the subsetted cells. Two distance metrics are then calculated for each clone:

* **Intra-clone distance**: the average Euclidean distance between cells and their clone centroid in the top *n* principal components, reflecting the transcriptional dispersion within a clone.
* **Inter-clone distance**: the average distance between cells from one clone and the centroids of all other clones, providing a measure of transcriptional separation between clones.

The function performs within-sample comparisons of intra- versus inter-clone distances using Wilcoxon tests, as well as between-sample comparisons (e.g., young versus old). Results are returned as a faceted violin plot with significance annotations, along with summary statistics and the underlying pairwise distance matrix.

### `compute_variance_metrics()`

Identifies genes with differential transcriptional variability between two groups using Brown–Forsythe and Levene’s tests. The function also calculates mean-adjusted variance by modeling the global mean–variance relationship with LOESS, helping distinguish true variability changes from expression-level effects. Outputs include raw and adjusted variance estimates, log2 fold-changes, adjusted p-values, and summary statistics for downstream DVG analysis and pathway enrichment.


### `plot_variance_summaries()`

Generates a visualization report from the output of `compute_variance_metrics()`, including variance boxplots, mean-adjusted variance distributions, raw and adjusted variance volcano plots, and summary tables of significant genes across multiple log2FC thresholds.


### `plot_density_gene()`

This function visualizes the expression distribution of a selected gene across two groups using density plots, excluding cells with zero expression for clearer interpretation. When provided with results from `compute_variance_metrics()`, it also displays mean-adjusted variance for each group, allowing users to assess differences in both expression level and transcriptional variability.

### `run_clone_distribution_engine()`

This function visualizes clone composition across samples using stacked bar plots of relative clone frequency. Clones below a user-defined threshold are grouped into an "Other" category. The output includes individual plots, a combined figure, and summary statistics for clone frequencies within each sample.

### `run_low_high_OA_analysis_for_a_single_sample()`

This function performs a complete Output Activity (OA) analysis workflow, including OA score calculation, UMAP visualization, differential expression analysis between Low- and High-OA clones, marker gene enrichment, and clone density mapping. All results, figures, and summary tables are saved to a sample-specific output directory.

### `run_msigdb_gsea_all_collections()`

This function performs pathway enrichment analysis on ranked DVG or DEG results using MSigDB gene sets and `clusterProfiler::GSEA()`. Gene symbols are mapped to Entrez IDs, and enrichment is evaluated across selected MSigDB collections. For each collection, the top significantly enriched pathways are visualized with enrichment plots, and leading-edge genes are converted back to gene symbols. The function returns pathway-level results, summary visualizations, and a consolidated GSEA results table.

## Demo Dataset

The package ships with `scCloneVar_test_demo`, a small Seurat object containing 20 Young and 20 Old clones from an in vitro cross-age dataset (with RNA assay, PCA, UMAP, and `CloneID` / `donor_age` metadata):

```r
data(scCloneVar_test_demo)
```

## Environment

Developed and tested under **R 4.4.1** on macOS (x86_64-apple-darwin20).

### Package dependencies

| Package | Version | Role |
|---|---|---|
| Seurat | `5.3.0` | Single-cell object, HVGs, PCA, FindMarkers |
| Matrix | 1.7-3 | Sparse matrix operations |
| dplyr | 1.1.4 | Data manipulation |
| tidyr | 1.3.1 | Reshaping |
| tibble | 3.2.1 | Data frames |
| purrr | 1.0.4 | Functional iteration |
| magrittr | 2.0.3 | Pipe |
| stringr | 1.5.1 | String handling |
| ggplot2 | 4.0.0 | Plotting |
| patchwork | 1.3.0 | Plot composition |
| ggrepel | 0.9.6 | Volcano labels |
| ggplotify | 0.1.2 | Plot conversion |
| ggsignif | `0.6.4` | Significance brackets |
| ggsci | `3.2.0` | Color palettes |
| ggvenn | `0.1.19` | Venn diagrams |
| scales | 1.4.0 | Axis scaling |
| RColorBrewer | 1.1-3 | Palettes |
| colorspace | 2.1-1 | Color manipulation |
| clusterProfiler | 4.12.6 | GSEA engine |
| enrichplot | 1.24.4 | Enrichment plots |
| msigdbr | 25.1.1 | MSigDB gene sets |
| org.Hs.eg.db | 3.19.1 | Human annotation |
| org.Mm.eg.db | `3.19.1` | Mouse annotation |
| openxlsx | 4.2.8 | Excel export |
| progress | 1.2.3 | Progress bars |
| progressr | `0.15.1` | Progress signaling |
| glue | 1.8.0 | String interpolation |



