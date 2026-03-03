# scCloneVar

Clone-wise heterogeneity and differential expression analysis toolkit for single-cell RNA-seq data.

## Installation

```r
install.packages("devtools")
devtools::install_github("ChrisXC25/scCloneVar")
library(scCloneVar)

## Functions Overview
run_clonewise_DEG_suite()

Performs clone-wise differential expression analysis between two clone sets using multiple gene universes and statistical methods; input: Seurat object and clone ID sets; output: DEG tables, summary statistics, and volcano plots.

run_low_high_OA_analysis_for_a_single_sample()

Computes clone-level Output Activity (OA), classifies clones into low/high output groups, performs DEG analysis, and generates integrated visualization reports; input: Seurat object with clone metadata; output: OA tables, DEG results, UMAPs, volcano plots, Venn diagrams, and saved report files.

run_clone_distribution_engine()

Generates clone size distribution plots and descriptive statistics across user-defined sample comparisons; input: list of Seurat objects and comparison settings; output: stacked bar plots and summary statistics tables.

compute_pca_clone_heterogeneity()

Quantifies intra- and inter-clonal distances in PCA space to assess transcriptional heterogeneity; input: Seurat object with PCA reduction; output: clone-level distance metrics and summary tables.

compute_variance_and_cv()

Calculates gene-level variance, coefficient of variation (CV), and mean-adjusted variance across groups; input: expression matrix or Seurat object; output: statistical test results and heterogeneity metrics.

plot_variance_cv_summaries()

Visualizes differential variance genes (DVGs) and CV distributions with statistical annotations; input: DVG results table; output: heterogeneity summary plots.

run_msigdb_gsea_all_collections()

Performs GSEA across multiple MSigDB collections using ranked gene statistics; input: ranked gene table; output: enrichment results tables and pathway summaries.
