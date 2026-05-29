<h1>
  scCloneVar
  <img src="Figures/figure_v1.png" align="right" width="180">
</h1>

Clone-wise heterogeneity and differential expression analysis toolkit for single-cell RNA-seq data.

## Functions Overview

scCloneVar provides an integrated framework for analyzing clonal structure and transcriptional heterogeneity in single-cell RNA-seq data. The toolkit enables clone-wise differential expression analysis, allowing users to compare predefined clone groups across multiple gene universes and statistical models. It quantifies intra- and inter-clonal transcriptional dispersion in PCA space and identifies differential variance genes (DVGs) to capture heterogeneity beyond mean expression changes. The package also supports Output Activity (OA) analysis to classify clones by lineage bias and generate comprehensive visualization reports. Finally, pathway-level interpretation is facilitated through MSigDB-based GSEA, enabling biological contextualization of clone-associated signatures.

## Installation

```r
install.packages("devtools")
devtools::install_github("LabShengLi/scCloneVar")
library(scCloneVar)
```

## Required Packages

scCloneVar builds on a number of established single-cell and statistical R packages. Most of these will be pulled in automatically as dependencies, but if you run into installation issues it's usually faster to install them yourself first.

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

The main differential expression workhorse. Given a Seurat object and two sets of clone IDs (Group1 vs Group2), this runs DEG testing across **three gene universes** — all genes, the top *n* HVGs, and a filtered HVG set (genes passing minimum detection rate, mean normalized expression, and mean raw count thresholds) — and **two methods** (Wilcoxon and MAST). For each of the six resulting comparisons it computes BH and BY adjusted p-values, appends per-group average expression, and tabulates the number of significant genes at a user-defined log2FC × FDR cutoff.

Everything gets dumped into a single Excel workbook with a `Summary` sheet plus one sheet per comparison, and the MAST × Top-HVGs combination is rendered as a labeled volcano plot (top 20 genes per side). Outputs go into a timestamped subfolder.

### `compute_pca_clone_heterogeneity()`

This is the function for quantifying transcriptional dispersion in PCA space. After subsetting to cell types of interest (default: LT-HSC and ST-HSC) and re-running HVG selection → scaling → PCA on the subset, it computes two distance measures per clone:

- **Intra-clone distance** — average Euclidean distance from each cell to its own clone centroid in the top *n* PCs. This is essentially "how spread out is this clone in transcriptional space?"
- **Inter-clone distance** — average distance from cells of one clone to the centroids of all *other* clones. This gives you a sense of how distinct clones are from each other.

It then runs within-sample Wilcoxon tests (intra vs inter) and between-sample tests (e.g., young vs old) and returns a faceted violin plot with significance brackets, plus the underlying summary tables and pairwise distance matrix.

### `compute_variance_metrics()`

For each HVG, this tests whether transcriptional variance differs between two groups (e.g., young vs old, or low- vs high-OA clones) using both **Brown–Forsythe** (median-centered, more robust) and **Levene's test** (mean-centered). On top of the raw variance, it also fits a LOESS curve to the global mean–variance relationship and returns a **mean-adjusted variance** — basically the residual after regressing out the expected variance for a gene's expression level. This matters because in scRNA-seq, variance scales strongly with mean expression, so without adjusting you'll just keep rediscovering highly expressed genes.

Returns a per-gene table with raw and adjusted variance for each group, log2 fold-changes of both, and BH-adjusted p-values for both tests. The output of this function is what feeds directly into the DVG (differential variance gene) volcano plots and GSEA.

### `plot_variance_summaries()`

Takes the data frame produced by `compute_variance_metrics()` and generates a full visual report:

- Per-gene variance boxplot (raw variance, both groups, Wilcoxon-tested)
- Mean-adjusted variance boxplot
- Density plot of mean-adjusted variance distributions
- Two volcano plots — one for raw variance, one for mean-adjusted — with the top 10 genes on each side labeled
- A summary count table reporting how many genes pass at several log2FC thresholds (0, 0.1, 0.25, 0.5)

### `plot_density_gene()`

A focused per-gene plot for sanity-checking individual hits. Given a gene name and a grouping column, it draws expression density curves for the two groups (cells with zero expression are dropped so the distributions are interpretable). Optionally, if you pass in the DVG table from `compute_variance_metrics()`, it appends a small bar plot of mean-adjusted variance next to the density — useful for showing simultaneously that a gene differs in both *level* and *spread*.

### `run_clone_distribution_engine()`

A wrapper for visualizing how clones are distributed across samples and replicates. You give it a list of comparisons (each specifying a Seurat object, the samples to include, a minimum frequency threshold, and a title) and it produces stacked bar plots of relative clone size, lumping all sub-threshold clones into a single "Other" category. A custom palette builder makes sure colors are visually distinct (it filters out low-saturation colors from a pooled palette of Brewer, NPG, and D3 schemes) and "Other" is always grey. Returns the individual plots, a combined `patchwork` panel, and descriptive statistics (number of clones, min/max/mean frequency) per sample.

### `run_low_high_OA_analysis_for_a_single_sample()`

1. **Cell type UMAP** — a labeled `DimPlot` saved as PDF.
2. **Output Activity (OA) computation** — for each clone, OA = (relative non-HSC frequency) / (relative HSC frequency). Clones in the bottom 30% are tagged as Low-OA (HSC-biased, low differentiation output) and the top 30% as High-OA (differentiation-biased). Per-clone and per-cell tables are saved as Excel.
3. **OA UMAP** — cells colored by their clone's OA on a diverging red→blue scale, with cell-type centroids overlaid.
4. **Clone-wise DEG** — internally calls `run_clonewise_DEG_suite()` on Low-OA vs High-OA clones, restricted to specified cell types (default: LT-HSC, ST-HSC, MPP).
5. **Reference-marker volcano** — overlays user-supplied low/high-output marker genes on the MAST volcano and labels the ones that come out significant.
6. **Universe-filtered Venn diagram** — intersects the significant Low-OA DEGs with a reference low-output marker list, restricted to the shared gene universe (so the overlap is statistically meaningful), and reports a hypergeometric p-value.
7. **Contour UMAP** — density contours overlaid on the cell type UMAP, useful for showing where clones concentrate.

Everything gets written to `output_dir/sample_label/` along with an `.rds` of the full results object.

### `run_msigdb_gsea_all_collections()`

For pathway-level interpretation of the DVG (or DEG) signal. Takes a results data frame with a `gene` column and a ranking variable (default: `log2FC_mean_adjusted_variance_YO`), maps SYMBOL → ENTREZID via `org.Mm.eg.db` or `org.Hs.eg.db`, and runs `clusterProfiler::GSEA()` across a configurable list of MSigDB collections (defaults cover H, C2 subcollections, C3, C4, C5, C6, C7, C8).

For each collection, the top *n* up- and down-regulated pathways (by NES, at FDR < 0.05) get rendered as `gseaplot2` enrichment plots, and leading-edge ENTREZ IDs are converted back to gene symbols for the summary table. Returns per-collection results, combined up/down panel plots across all collections, and a single tidy GSEA summary table.

## Demo Dataset

The package ships with `scCloneVar_test_demo`, a small Seurat object containing 20 Young and 20 Old clones from an in vitro cross-age dataset (with RNA assay, PCA, UMAP, and `CloneID` / `donor_age` metadata). Use it to test the workflow before running on your own data:

```r
data(scCloneVar_test_demo)
```
