#' Run clone-wise differential expression suite
#'
#' Runs clone-wise DEG between two clone sets using multiple gene universes
#' (all genes, top HVGs, filtered HVGs) and multiple DEG methods (wilcox, MAST).
#' Saves an Excel workbook (summary + per-result tabs) and optionally a volcano plot.
#'
#' @param seurat_obj A Seurat object.
#' @param clone_set1_ids Vector of clone IDs for Group1.
#' @param clone_set2_ids Vector of clone IDs for Group2.
#' @param clone_col Metadata column name storing clone IDs.
#' @param assay_name Assay to use.
#' @param top_hvgs Number of HVGs.
#' @param min_detect_prop Minimum detection proportion.
#' @param min_mean_norm Minimum mean normalized expression.
#' @param min_mean_raw Minimum mean raw counts.
#' @param excel_name Excel output file name.
#' @param volcano_name Volcano plot file name.
#' @param group1_label Label for Group1.
#' @param group2_label Label for Group2.
#' @param fc_cutoff log2FC cutoff.
#' @param fdr_cutoff Adjusted p-value cutoff.
#' @param output_dir Output directory.
#'
#' @return Invisible list containing summary, results, volcano plot, and folder path.
#' @export
run_clonewise_DEG_suite <- function(
    seurat_obj,
    clone_set1_ids,
    clone_set2_ids,
    clone_col   = "CloneID",
    assay_name  = "RNA",
    top_hvgs    = 2000,
    min_detect_prop = 0.10,
    min_mean_norm  = 0.10,
    min_mean_raw   = 2,
    excel_name  = "Clonewise_DEG.xlsx",
    volcano_name = "Volcano_MAST_TopHVGs.png",
    group1_label = "Low_OA",
    group2_label = "High_OA",
    fc_cutoff  = 0.25,
    fdr_cutoff = 0.05,
    output_dir = "Clonewise_DEG_Results"
) {

  # -----------------------------
  # Safety checks
  # -----------------------------
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }

  if (!(assay_name %in% names(seurat_obj@assays))) {
    stop("Assay not found in seurat_obj: ", assay_name)
  }

  if (!(clone_col %in% colnames(seurat_obj@meta.data))) {
    stop("clone_col not found in metadata: ", clone_col)
  }

  Seurat::DefaultAssay(seurat_obj) <- assay_name

  # -----------------------------
  # Create output folder
  # -----------------------------
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_folder <- file.path(
    output_dir,
    paste0(group1_label, "_vs_", group2_label, "_", timestamp)
  )

  dir.create(run_folder, recursive = TRUE, showWarnings = FALSE)

  excel_path   <- file.path(run_folder, excel_name)
  volcano_path <- file.path(run_folder, volcano_name)

  message("Output directory: ", run_folder)

  # -----------------------------
  # Assign groups
  # -----------------------------
  clone_ids <- seurat_obj@meta.data[[clone_col]]

  group_label <- rep(NA_character_, length(clone_ids))
  group_label[clone_ids %in% clone_set1_ids] <- "Group1"
  group_label[clone_ids %in% clone_set2_ids] <- "Group2"

  seurat_obj$group_label <- group_label

  keep_cells <- !is.na(seurat_obj$group_label)
  seurat_sub <- seurat_obj[, keep_cells]

  Seurat::Idents(seurat_sub) <- seurat_sub$group_label

  n1 <- sum(seurat_sub$group_label == "Group1")
  n2 <- sum(seurat_sub$group_label == "Group2")

  message("Cells: ", n1, " vs ", n2)

  # -----------------------------
  # Helper: adjust + avg expression
  # -----------------------------
  adjust_and_append_avg <- function(deg_df, obj) {

    if (is.null(deg_df) || nrow(deg_df) == 0) {
      return(data.frame())
    }

    deg_df$gene <- rownames(deg_df)

    deg_df$p_val_adj <- stats::p.adjust(deg_df$p_val, method = "BH")
    deg_df$FDR_BY_manual <- stats::p.adjust(deg_df$p_val, method = "BY")

    ae <- Seurat::AverageExpression(obj, group.by = "group_label")
    ae_df <- as.data.frame(ae[[assay_name]])
    ae_df$gene <- rownames(ae_df)

    colnames(ae_df)[colnames(ae_df) == "Group1"] <- "avg_Group1"
    colnames(ae_df)[colnames(ae_df) == "Group2"] <- "avg_Group2"

    merge(deg_df, ae_df, by = "gene", all.x = TRUE)
  }

  # -----------------------------
  # DEG runner
  # -----------------------------
  do_one_deg <- function(obj, features, method) {
    if (length(features) == 0) {
      return(data.frame())
    }

    Seurat::FindMarkers(
      obj,
      ident.1 = "Group1",
      ident.2 = "Group2",
      features = features,
      test.use = method,
      logfc.threshold = 0,
      min.pct = 0
    )
  }

  # -----------------------------
  # Gene universes
  # -----------------------------
  all_genes <- rownames(seurat_sub)

  seurat_sub <- Seurat::FindVariableFeatures(
    seurat_sub,
    selection.method = "vst",
    nfeatures = top_hvgs,
    verbose = FALSE
  )

  hvgs <- Seurat::VariableFeatures(seurat_sub)

  data_mat <- Seurat::GetAssayData(seurat_sub, slot = "data")
  raw_mat  <- tryCatch(
    Seurat::GetAssayData(seurat_sub, slot = "counts"),
    error = function(e) NULL
  )

  hvgs <- intersect(hvgs, rownames(data_mat))

  det_prop  <- Matrix::rowMeans(data_mat[hvgs, , drop = FALSE] > 0)
  mean_norm <- Matrix::rowMeans(data_mat[hvgs, , drop = FALSE])

  if (!is.null(raw_mat)) {
    mean_raw <- Matrix::rowMeans(raw_mat[hvgs, , drop = FALSE])
  } else {
    mean_raw <- rep(0, length(hvgs))
  }

  filt_genes <- hvgs[
    det_prop >= min_detect_prop &
      mean_norm >= min_mean_norm &
      mean_raw >= min_mean_raw
  ]

  # -----------------------------
  # DEG loop
  # -----------------------------
  gene_sets <- list(
    DEG1_AllGenes = all_genes,
    DEG2_TopHVGs = hvgs,
    DEG3_FilteredHVGs = filt_genes
  )

  methods <- c("wilcox", "MAST")

  deg_results <- list()
  summary_list <- list()

  for (setting in names(gene_sets)) {

    feats <- gene_sets[[setting]]

    for (m in methods) {

      message("Running ", m, " on ", setting)

      res <- do_one_deg(seurat_sub, feats, m)
      res2 <- adjust_and_append_avg(res, seurat_sub)

      key <- paste0(setting, "_", toupper(m))
      deg_results[[key]] <- res2

      if (nrow(res2) > 0) {
        n_sig <- sum(res2$p_val_adj < fdr_cutoff, na.rm = TRUE)
        n_up1 <- sum(res2$avg_log2FC >  fc_cutoff & res2$p_val_adj < fdr_cutoff, na.rm = TRUE)
        n_up2 <- sum(res2$avg_log2FC < -fc_cutoff & res2$p_val_adj < fdr_cutoff, na.rm = TRUE)
      } else {
        n_sig <- 0; n_up1 <- 0; n_up2 <- 0
      }

      summary_list[[key]] <- data.frame(
        Setting = setting,
        Method  = toupper(m),
        Cells_Group1 = n1,
        Cells_Group2 = n2,
        Genes_Test = length(feats),
        Sig_FDR = n_sig,
        Up_in_Group1 = n_up1,
        Up_in_Group2 = n_up2
      )
    }
  }

  summary_tab <- do.call(rbind, summary_list)

  # -----------------------------
  # Save Excel
  # -----------------------------
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Summary")
  openxlsx::writeData(wb, "Summary", summary_tab)

  for (nm in names(deg_results)) {
    openxlsx::addWorksheet(wb, nm)
    df <- deg_results[[nm]]
    if (nrow(df) == 0) {
      openxlsx::writeData(wb, nm, data.frame(note = "No results"))
    } else {
      openxlsx::writeData(wb, nm, df)
    }
  }

  openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)

  # -----------------------------
  # Volcano plot (MAST × Top HVGs)
  # -----------------------------
  mast_key <- "DEG2_TopHVGs_MAST"
  volc_plot <- NULL
  mast_df <- deg_results[[mast_key]]

  if (!is.null(mast_df) && nrow(mast_df) > 0) {

    df <- mast_df
    df$neg_log10_fdr <- -log10(pmax(df$p_val_adj, 1e-300))

    df$significance <- "Not significant"
    df$significance[
      df$p_val_adj < fdr_cutoff & df$avg_log2FC > fc_cutoff
    ] <- "↑Group1"

    df$significance[
      df$p_val_adj < fdr_cutoff & df$avg_log2FC < -fc_cutoff
    ] <- "↑Group2"

    # Top genes selection
    top_left <- df[
      df$significance == "↑Group2",
    ]
    top_left <- top_left[order(top_left$p_val_adj), ]
    top_left <- head(top_left, 20)

    top_right <- df[
      df$significance == "↑Group1",
    ]
    top_right <- top_right[order(top_right$p_val_adj), ]
    top_right <- head(top_right, 20)

    top_genes <- rbind(top_left, top_right)

    volc_plot <- ggplot2::ggplot(
      df,
      ggplot2::aes(avg_log2FC, neg_log10_fdr, color = significance)
    ) +
      ggplot2::geom_point(alpha = 0.8) +
      ggrepel::geom_text_repel(
        data = top_genes,
        ggplot2::aes(label = gene),
        size = 4,
        box.padding = 0.4,
        point.padding = 0.4,
        max.overlaps = Inf
      ) +
      ggplot2::scale_color_manual(values = c(
        "↑Group1" = "#D85B59",
        "↑Group2" = "#5271AE",
        "Not significant" = "grey80"
      )) +
      ggplot2::geom_vline(
        xintercept = c(-fc_cutoff, fc_cutoff),
        linetype = "dashed"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(fdr_cutoff),
        linetype = "dashed"
      ) +
      ggplot2::theme_classic(base_size = 16) +
      ggplot2::labs(
        title = paste0("Clone-wise DEG Volcano (MAST, Top ", top_hvgs, " HVGs)"),
        x = expression(log[2]("FC (Group1 / Group2)")),
        y = expression(-log[10]("FDR")),
        color = NULL
      )

    ggplot2::ggsave(
      volcano_path,
      volc_plot,
      width = 8.5,
      height = 7,
      dpi = 300
    )
  }

  invisible(list(
    summary = summary_tab,
    results = deg_results,
    volcano_plot = volc_plot,
    folder = run_folder
  ))
}
