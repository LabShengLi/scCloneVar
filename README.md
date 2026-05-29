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

The main differential expression analysis function. Given a Seurat object and two groups of clone IDs (Group 1 vs Group 2), it performs differential expression testing across three gene sets: all genes, the top *n* highly variable genes (HVGs), and a filtered HVG set based on minimum detection rate, mean normalized expression, and mean raw count thresholds. Two statistical methods are supported: Wilcoxon and MAST. For each analysis, the function calculates both Benjamini–Hochberg (BH) and Benjamini–Yekutieli (BY) adjusted p-values, adds average expression values for each group, and summarizes the number of significant genes based on user-defined log2 fold-change and FDR thresholds.

All results are saved to a single Excel workbook containing a `Summary` sheet and separate sheets for each comparison. A labeled volcano plot is generated for the MAST analysis on the top HVG gene set, highlighting the top 20 upregulated and downregulated genes. All output files are organized into a timestamped subdirectory for easy tracking and reproducibility.

### `compute_pca_clone_heterogeneity()`

This function quantifies transcriptional heterogeneity in PCA space. After subsetting the data to selected cell types (default: LT-HSC and ST-HSC), it performs highly variable gene selection, scaling, and PCA on the subsetted cells. Two distance metrics are then calculated for each clone:

* **Intra-clone distance**: the average Euclidean distance between cells and their clone centroid in the top *n* principal components, reflecting the transcriptional dispersion within a clone.
* **Inter-clone distance**: the average distance between cells from one clone and the centroids of all other clones, providing a measure of transcriptional separation between clones.

The function performs within-sample comparisons of intra- versus inter-clone distances using Wilcoxon tests, as well as between-sample comparisons (e.g., young versus old). Results are returned as a faceted violin plot with significance annotations, along with summary statistics and the underlying pairwise distance matrix.

### `compute_variance_metrics()`

This function identifies genes with differential transcriptional variability between two groups (e.g., young vs. old, or low- versus high-OA clones). For each highly variable gene (HVG), variability is assessed using both the **Brown–Forsythe test**, which is based on deviations from the median and is more robust to outliers, and **Levene’s test**, which is based on deviations from the mean.

To account for the strong relationship between mean expression and variance in single-cell RNA-seq data, the function also models the global mean–variance trend using LOESS. A **mean-adjusted variance** is then calculated as the residual variance after removing the expected variance associated with a gene’s expression level. This adjustment helps identify genes with unusually high or low variability beyond what would be expected from their mean expression alone.

The function returns a gene-level results table containing raw and mean-adjusted variance estimates for each group, log2 fold-changes of both variance measures, and Benjamini–Hochberg (BH) adjusted p-values for both statistical tests. These results serve as the primary input for downstream differential variability gene (DVG) visualization and pathway enrichment analyses.


### `plot_variance_summaries()`

This function generates a comprehensive visualization report from the output of `compute_variance_metrics()`. The report includes:

* A per-gene boxplot of raw variance for each group, with significance assessed using a Wilcoxon test.
* A boxplot of mean-adjusted variance, allowing comparison of variability after accounting for the mean–variance relationship.
* Density plots showing the distribution of mean-adjusted variance across groups.
* Two volcano plots, one based on raw variance and the other on mean-adjusted variance, with the top 10 genes showing the largest increases and decreases in variability labeled.
* A summary table reporting the number of significant genes at multiple log2 fold-change thresholds (0, 0.1, 0.25, and 0.5).

Together, these visualizations provide an overview of differential transcriptional variability between groups and help identify genes that contribute most strongly to changes in cellular heterogeneity.


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
