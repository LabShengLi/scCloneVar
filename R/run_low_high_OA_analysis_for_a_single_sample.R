#' OA function (per-sample input)
#'
#' @importFrom magrittr %>%
#' @export
run_low_high_OA_analysis_for_a_single_sample <- function(
    seurat_obj,
    sample_label = "Sample",
    output_dir = "OA_Analysis",
    clone_col = "CloneID",
    HSC_types = c("LT-HSC","ST-HSC"),
    nonHSC_types = NULL,
    deg_celltypes = c("LT-HSC","ST-HSC","MPP"),
    top_hvgs = 2000,
    min_detect_prop = 0.10,
    min_mean_norm  = 0.10,
    min_mean_raw   = 2,
    fc_cutoff      = 0.25,
    fdr_cutoff     = 0.05,
    low_output_markers,
    high_output_markers,
    celltype_colors,
    full_ref_deg_genes_list
) {
  
  message("\n==============================")
  message(glue::glue("🔷 Starting OA Analysis: {sample_label}"))
  message("==============================\n")
  
  # -----------------------------------------
  # 0. Create output folder
  # -----------------------------------------
  sample_dir <- file.path(output_dir, sample_label)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  message(glue::glue("📁 Output folder created: {sample_dir}\n"))
  
  # -----------------------------------------
  # 1. DimPlot
  # -----------------------------------------
  message("🔹 Step 1: Creating DimPlot …")
  
  p_dim <- Seurat::DimPlot(
    seurat_obj,
    reduction = "umap",
    cols      = celltype_colors,
    pt.size   = 0.4,
    label     = TRUE,
    label.size = 6
  ) +
    ggplot2::ggtitle(glue::glue("{sample_label} — Celltype UMAP")) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = 16)
    )
  
  ggplot2::ggsave(
    filename = file.path(sample_dir, glue::glue("{sample_label}_DimPlot.pdf")),
    plot = p_dim, width = 8.5, height = 7.5
  )
  
  message("Done: DimPlot saved.\n")
  
  # -----------------------------------------
  # 2. Compute OA per clone
  # -----------------------------------------
  message("🔹 Step 2: Computing Output Activity (OA)…")
  
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_ID") %>%
    dplyr::rename(CloneID = !!clone_col)
  
  if (is.null(nonHSC_types)) {
    nonHSC_types <- setdiff(unique(meta$celltype), HSC_types)
  }
  
  message(glue::glue("   HSC types: {paste(HSC_types, collapse=', ')}"))
  message(glue::glue("   non-HSC types: {paste(nonHSC_types, collapse=', ')}"))
  
  meta <- meta %>%
    dplyr::filter(CloneID != "0", !is.na(CloneID))
  
  clone_counts <- meta %>%
    dplyr::mutate(group = ifelse(celltype %in% HSC_types, "HSC", "nonHSC")) %>%
    dplyr::count(CloneID, group) %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = n,
      values_fill = 0
    )
  
  clone_freq <- clone_counts %>%
    dplyr::mutate(
      total_HSC_all = sum(HSC),
      total_nonHSC_all = sum(nonHSC),
      HSC_freq = HSC / total_HSC_all,
      nonHSC_freq = nonHSC / total_nonHSC_all,
      OA = (nonHSC_freq + 1e-6) / (HSC_freq + 1e-6),
      log2OA = log2(OA)
    )
  
  q_low  <- stats::quantile(clone_freq$OA, 0.30)
  q_high <- stats::quantile(clone_freq$OA, 0.70)
  
  low_clones  <- clone_freq %>% dplyr::filter(OA <= q_low)  %>% dplyr::pull(CloneID)
  high_clones <- clone_freq %>% dplyr::filter(OA >= q_high) %>% dplyr::pull(CloneID)
  
  message(glue::glue("   - Low Output clones: {length(low_clones)}"))
  message(glue::glue("   - High Output clones: {length(high_clones)}"))
  
  openxlsx::write.xlsx(
    clone_freq,
    file = file.path(sample_dir, glue::glue("{sample_label}_OA_by_clone.xlsx")),
    overwrite = TRUE
  )
  
  cell_OA_df <- meta %>%
    dplyr::left_join(clone_freq %>% dplyr::select(CloneID, OA), by = "CloneID") %>%
    dplyr::select(cell_ID, CloneID, OA)
  
  openxlsx::write.xlsx(
    cell_OA_df,
    file = file.path(sample_dir, glue::glue("{sample_label}_OA_by_cell.xlsx")),
    overwrite = TRUE
  )
  
  message("Done: OA tables saved.\n")
  
  # -----------------------------------------
  # 3. OA UMAP plot
  # -----------------------------------------
  message("🔹 Step 3: Creating OA UMAP …")
  
  umap_df <- Seurat::Embeddings(seurat_obj, "umap") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_ID")
  
  meta_umap <- meta %>%
    dplyr::left_join(umap_df, by = "cell_ID") %>%
    dplyr::left_join(clone_freq %>% dplyr::select(CloneID, OA), by = "CloneID") %>%
    dplyr::filter(!is.na(OA))
  
  meta_umap$OA <- scales::squish(meta_umap$OA, c(0, 2))
  
  cell_centroids <- meta_umap %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  
  p_oa <- ggplot2::ggplot(meta_umap, ggplot2::aes(UMAP_1, UMAP_2, color = OA)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.55) +
    ggplot2::scale_color_gradientn(
      colors = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"),
      values = scales::rescale(c(0, 0.7, 1, 1.3, 2)),
      limits = c(0, 2),
      name = "OA"
    ) +
    ggplot2::geom_text(
      data = cell_centroids,
      ggplot2::aes(label = celltype),
      size = 5.5, fontface = "bold", color = "grey10"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 16) +
    ggplot2::ggtitle(glue::glue("{sample_label} — OA (Ai) UMAP"))
  
  ggplot2::ggsave(
    filename = file.path(sample_dir, glue::glue("{sample_label}_OA_UMAP.pdf")),
    plot = p_oa, width = 6.5, height = 5.5
  )
  
  message("Done: OA UMAP saved.\n")
  
  # -----------------------------------------
  # 4. DEG
  # -----------------------------------------
  message("🔹 Step 4: Running clonewise DEG (HSC-only)…")
  
  seurat_deg <- Seurat::subset(seurat_obj, subset = celltype %in% deg_celltypes)
  
  deg_out <- run_clonewise_DEG_suite(
    seurat_obj       = seurat_deg,
    clone_set1_ids   = low_clones,
    clone_set2_ids   = high_clones,
    clone_col        = clone_col,
    group1_label     = "Low_output",
    group2_label     = "High_output",
    top_hvgs         = top_hvgs,
    min_detect_prop  = min_detect_prop,
    min_mean_norm    = min_mean_norm,
    min_mean_raw     = min_mean_raw,
    fc_cutoff        = fc_cutoff,
    fdr_cutoff       = fdr_cutoff,
    excel_name       = glue::glue("{sample_label}_DEG.xlsx"),
    volcano_name     = glue::glue("{sample_label}_volcano.pdf"),
    output_dir       = sample_dir
  )
  
  message("Done: DEG analysis completed.\n")
  message("   Notes: Adding reference-marker volcano plot ...")
  
  if (!exists("plots")) plots <- list()
  
  mast_top <- deg_out$results$DEG2_TopHVGs_MAST
  
  if (!is.null(mast_top) && nrow(mast_top) > 0) {
    
    reference_genes <- unique(c(low_output_markers, high_output_markers))
    
    df_ref <- mast_top %>%
      dplyr::mutate(
        neg_log10_fdr = -log10(pmax(p_val_adj, 1e-300)),
        significance = dplyr::case_when(
          p_val_adj < fdr_cutoff & avg_log2FC >  fc_cutoff ~ "↑Low_OA",
          p_val_adj < fdr_cutoff & avg_log2FC < -fc_cutoff ~ "↑High_OA",
          TRUE ~ "Not significant"
        ),
        is_reference = gene %in% reference_genes
      )
    
    label_genes <- df_ref %>%
      dplyr::filter(is_reference & significance != "Not significant")
    
    p_ref_volcano <- ggplot2::ggplot(df_ref, ggplot2::aes(avg_log2FC, neg_log10_fdr)) +
      ggplot2::geom_point(ggplot2::aes(color = significance), alpha = 0.85, size = 2.5) +
      ggrepel::geom_text_repel(
        data = label_genes,
        ggplot2::aes(label = gene),
        size = 5.3,
        color = "black",
        box.padding = 0.2,
        point.padding = 0.1,
        segment.size = 0.25
      ) +
      ggplot2::scale_color_manual(values = c(
        "↑Low_OA"  = "#D85B59",
        "↑High_OA" = "#5271AE",
        "Not significant" = "grey80"
      )) +
      ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                          linetype = "dashed", linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = -log10(fdr_cutoff),
                          linetype = "dashed", linewidth = 0.4) +
      ggplot2::coord_cartesian(xlim = c(-5, 5)) +
      ggplot2::theme_classic(base_size = 16) +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18),
        axis.title.y = ggplot2::element_text(size = 18)
      ) +
      ggplot2::labs(
        x = expression(log[2]("Fold change")),
        y = expression(-log[10]("adj. P value"))
      )
    
    ref_vol_path <- file.path(
      sample_dir,
      glue::glue("{sample_label}_reference_marker_volcano.pdf")
    )
    
    ggplot2::ggsave(ref_vol_path, p_ref_volcano,
                    width = 8, height = 6, dpi = 300)
    
    message(glue::glue("      ✔ Reference volcano saved: {ref_vol_path}"))
    plots$reference_volcano <- p_ref_volcano
    
  } else {
    message("  No MAST Top-HVGs DEG available — skipping reference volcano.")
  }
  
  # -----------------------------------------
  # 5. Universe-filtered Venn
  # -----------------------------------------
  message("🔹 Step 5: Creating universe-filtered Venn diagram …")
  
  mast_genes <- mast_top$gene
  ref_genes  <- full_ref_deg_genes_list
  gene_universe <- intersect(mast_genes, ref_genes)
  
  sig_low_DEGs <- mast_top %>%
    dplyr::filter(avg_log2FC > fc_cutoff, p_val_adj < fdr_cutoff) %>%
    dplyr::pull(gene)
  
  reference_set_universe <- intersect(low_output_markers, gene_universe)
  deg_set_universe       <- intersect(sig_low_DEGs, gene_universe)
  
  overlap_genes <- intersect(reference_set_universe, deg_set_universe)
  
  k <- length(overlap_genes)
  m <- length(reference_set_universe)
  q <- length(deg_set_universe)
  N <- length(gene_universe)
  
  p_hyper <- stats::phyper(
    q = k - 1,
    m = m,
    n = N - m,
    k = q,
    lower.tail = FALSE
  )
  
  message(glue::glue("   • Overlap size: {k}"))
  message(glue::glue("   • Hypergeometric P-value: {signif(p_hyper, 3)}"))
  
  venn_list_filtered <- list(
    low_output_markers = reference_set_universe,
    low_output_DEG     = deg_set_universe
  )
  
  p_venn_filtered <- ggvenn::ggvenn(
    venn_list_filtered,
    fill_color = c("#4E9AC7", "#F79A63"),
    text_size = 12,
    stroke_size = 0.6,
    show_percentage = FALSE
  ) +
    ggplot2::ggtitle(glue::glue("{sample_label}: Low Output — Universe-filtered")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5))
  
  ggplot2::ggsave(
    file.path(sample_dir, glue::glue("{sample_label}_LowOutput_venn_filtered.pdf")),
    p_venn_filtered,
    width = 6.5,
    height = 5.5
  )
  
  message("Done: Filtered Venn diagram saved.\n")
  
  # -----------------------------------------
  # 6. Contour UMAP
  # -----------------------------------------
  message("🔹 Step 6: Creating contour UMAP …")
  
  meta_df <- seurat_obj@meta.data %>%
    dplyr::select(celltype) %>%
    cbind(
      Seurat::Embeddings(seurat_obj, "umap") %>%
        as.data.frame() %>%
        stats::setNames(c("UMAP_1", "UMAP_2"))
    )
  
  p_contour <- ggplot2::ggplot(meta_df, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggrastr::rasterise(
      ggplot2::geom_point(
        ggplot2::aes(fill = celltype),
        shape = 21, color = "black",
        size = 1.2, stroke = 0.15, alpha = 0.8
      ),
      dpi = 300
    ) +
    ggplot2::geom_density_2d(color = "black", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = celltype_colors) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 18) +
    ggplot2::ggtitle(glue::glue("{sample_label} — Contour UMAP"))
  
  ggplot2::ggsave(
    filename = file.path(sample_dir, glue::glue("{sample_label}_contour_umap.pdf")),
    plot = p_contour,
    width = 9, height = 7, device = grDevices::cairo_pdf
  )
  
  message("Done: Contour UMAP saved.\n")
  
  # -----------------------------------------
  # 7. Save results object
  # -----------------------------------------
  results <- list(
    OA_by_clone = clone_freq,
    OA_by_cell = cell_OA_df,
    low_clones = low_clones,
    high_clones = high_clones,
    DEG_results = deg_out,
    folder = sample_dir,
    plot_dim = p_dim,
    plot_OA = p_oa,
    plot_DEG_volcano = deg_out$volcano_plot,
    plot_ref_volcano = if (exists("p_ref_volcano")) p_ref_volcano else NULL,
    plot_venn_filtered = p_venn_filtered,
    plot_contour = p_contour
  )
  
  saveRDS(results, file = file.path(sample_dir, "results_object.rds"))
  
  message(glue::glue("\n COMPLETE! All outputs saved to: {sample_dir}\n"))
  
  return(results)
}