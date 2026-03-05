#' Compute PCA-based clone heterogeneity (intra- and inter-clone distances)
#'
#' @importFrom magrittr %>%
#' @export
compute_pca_clone_heterogeneity <- function(
    seu,
    celltype_keep = c("LT-HSC", "ST-HSC"),
    clone_col = "CloneID_clean",
    sample_col = "sampleName",
    assay = "RNA",
    n_pcs = 30,
    min_pct = 0.1,
    min_mean = 0.1,
    min_cells = 2,
    group_recode = NULL,
    sample_order = NULL,
    prefix = "Clone Heterogeneity"
) {

  # ==========================================================
  # Step 1/5: Subset cells safely
  # ==========================================================
  message("Step 1/5: Subsetting cells...")

  if (!(clone_col %in% colnames(seu@meta.data)))
    stop("clone_col not found in metadata: ", clone_col)

  if (!(sample_col %in% colnames(seu@meta.data)))
    stop("sample_col not found in metadata: ", sample_col)

  if (!("celltype" %in% colnames(seu@meta.data)))
    stop("metadata column 'celltype' not found.")

  cells_keep <- rownames(seu@meta.data)[
    seu@meta.data$celltype %in% celltype_keep &
      seu@meta.data[[clone_col]] != "0"
  ]

  if (length(cells_keep) == 0)
    stop("No cells left after filtering.")

  seu <- subset(seu, cells = cells_keep)

  seu$CloneID <- seu@meta.data[[clone_col]]

  Seurat::DefaultAssay(seu) <- assay

  # ==========================================================
  # Step 2/5: HVG + PCA
  # ==========================================================
  message("Step 2/5: HVG filtering and PCA...")

  seu <- Seurat::FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = 2000
  )

  expr_data <- Seurat::GetAssayData(seu, layer = "data")

  hvgs <- Seurat::VariableFeatures(seu)

  gene_pct <- Matrix::rowMeans(expr_data[hvgs, , drop = FALSE] > 0)
  gene_avg <- Matrix::rowMeans(expr_data[hvgs, , drop = FALSE])

  genes_pass <- names(gene_pct)[
    gene_pct >= min_pct & gene_avg >= min_mean
  ]

  seu <- subset(seu, features = genes_pass)

  seu <- Seurat::ScaleData(seu, verbose = FALSE)

  seu <- Seurat::RunPCA(seu, npcs = n_pcs, verbose = FALSE)

  # ==========================================================
  # Step 3/5: Extract PCs
  # ==========================================================
  message("Step 3/5: Extracting PCs...")

  pcs <- Seurat::Embeddings(seu, "pca")[, 1:n_pcs, drop = FALSE]

  meta <- seu@meta.data %>%
    dplyr::mutate(cell_id = rownames(seu@meta.data))

  pc_df <- as.data.frame(pcs) %>%
    dplyr::mutate(cell_id = rownames(pcs)) %>%
    dplyr::inner_join(
      meta[, c("cell_id", "CloneID", sample_col)],
      by = "cell_id"
    ) %>%
    dplyr::rename(sampleName = dplyr::all_of(sample_col))

  # ==========================================================
  # Step 4/5: Intra-clone distance
  # ==========================================================
  message("Step 4/5: Computing intra-clone distances...")

  clone_ids_all <- unique(pc_df$CloneID)

  pb <- progress::progress_bar$new(
    format = "  intra [:bar] :percent eta: :eta",
    total = length(clone_ids_all),
    clear = FALSE
  )

  compute_clone_intra <- function(clone_id) {

    pb$tick()

    clone_df <- pc_df %>%
      dplyr::filter(CloneID == clone_id)

    if (nrow(clone_df) < 2) {
      return(dplyr::tibble(
        CloneID = clone_id,
        sampleName = clone_df$sampleName[1],
        n_cells = nrow(clone_df),
        intra_dist = NA_real_,
        centroid = list(NA)
      ))
    }

    pc_mat <- as.matrix(
      clone_df[, grep("^PC", colnames(clone_df)), drop = FALSE]
    )

    centroid <- colMeans(pc_mat)

    centered <- sweep(pc_mat, 2, centroid, "-")

    dists <- sqrt(rowSums(centered^2))

    dplyr::tibble(
      CloneID = clone_id,
      sampleName = clone_df$sampleName[1],
      n_cells = nrow(clone_df),
      intra_dist = mean(dists),
      centroid = list(centroid)
    )
  }

  intra_results <- purrr::map_dfr(
    clone_ids_all,
    compute_clone_intra
  )

  # ==========================================================
  # Step 5/5: Inter-clone distance
  # ==========================================================
  message("Step 5/5: Computing inter-clone distances...")

  intra_filtered <- intra_results %>%
    dplyr::filter(!is.na(intra_dist), n_cells >= min_cells)

  centroid_list <- intra_filtered %>%
    dplyr::select(CloneID, centroid) %>%
    tibble::deframe()

  clone_ids <- names(centroid_list)

  pairwise_df <- expand.grid(
    CloneA = clone_ids,
    CloneB = clone_ids,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(CloneA != CloneB)

  pb2 <- progress::progress_bar$new(
    format = "  inter [:bar] :percent eta: :eta",
    total = nrow(pairwise_df),
    clear = FALSE
  )

  calc_inter_dist <- function(clone_A, clone_B) {

    pb2$tick()

    df_A <- pc_df %>%
      dplyr::filter(CloneID == clone_A)

    pcs_A <- as.matrix(
      df_A[, grep("^PC", colnames(df_A)), drop = FALSE]
    )

    centroid_B <- centroid_list[[clone_B]]

    if (is.null(centroid_B) || nrow(pcs_A) == 0)
      return(NA_real_)

    dists <- sqrt(rowSums(
      (pcs_A -
         matrix(
           centroid_B,
           nrow = nrow(pcs_A),
           ncol = length(centroid_B),
           byrow = TRUE
         ))^2
    ))

    mean(dists)
  }

  pairwise_df$inter_dist <- purrr::map2_dbl(
    pairwise_df$CloneA,
    pairwise_df$CloneB,
    calc_inter_dist
  )

  avg_inter <- pairwise_df %>%
    dplyr::group_by(CloneA) %>%
    dplyr::summarise(
      mean_inter_dist = mean(inter_dist, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(CloneID = CloneA)

  summary_df <- intra_results %>%
    dplyr::left_join(avg_inter, by = "CloneID")

  # ==========================================================
  # Violin plot + statistics
  # ==========================================================
  clone_long <- summary_df %>%
    dplyr::select(sampleName, intra_dist, mean_inter_dist) %>%
    tidyr::pivot_longer(
      cols = c(intra_dist, mean_inter_dist),
      names_to = "DistanceType",
      values_to = "Distance"
    ) %>%
    dplyr::mutate(
      DistanceType = dplyr::recode(
        DistanceType,
        intra_dist = "Intra-clone",
        mean_inter_dist = "Inter-clone"
      )
    ) %>%
    tidyr::drop_na(Distance)

  if (!is.null(group_recode)) {
    clone_long$sampleName <-
      dplyr::recode(clone_long$sampleName, !!!group_recode)
  }
  if (!is.null(sample_order)) {

    clone_long$sampleName <- factor(
      clone_long$sampleName,
      levels = sample_order
    )

  } else {

    clone_long$sampleName <- factor(clone_long$sampleName)

  }
  clone_long$DistanceType <- factor(
    clone_long$DistanceType,
    levels = c("Intra-clone", "Inter-clone")
  )

  # Within-group tests
  within_tests <- clone_long %>%
    dplyr::group_by(sampleName) %>%
    dplyr::summarise(
      p_value =
        stats::wilcox.test(Distance ~ DistanceType,
                           exact = FALSE)$p.value,
      intra_mean =
        mean(Distance[DistanceType == "Intra-clone"],
             na.rm = TRUE),
      inter_mean =
        mean(Distance[DistanceType == "Inter-clone"],
             na.rm = TRUE),
      .groups = "drop"
    )

  # Between-group tests
  group_levels <- levels(clone_long$sampleName)

  between_tests <- purrr::map_dfr(
    combn(group_levels, 2, simplify = FALSE),
    function(pair) {
      clone_long %>%
        dplyr::filter(sampleName %in% pair) %>%
        dplyr::group_by(DistanceType) %>%
        dplyr::summarise(
          group1 = pair[1],
          group2 = pair[2],
          p_value =
            stats::wilcox.test(Distance ~ sampleName,
                               exact = FALSE)$p.value,
          .groups = "drop"
        )
    }
  )

  # ==========================================================
  # Prepare significance labels
  # ==========================================================

  p_to_star <- function(p) {
    dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  }

  within_tests$label <- p_to_star(within_tests$p_value)
  between_tests$label <- p_to_star(between_tests$p_value)

  y_max <- max(clone_long$Distance, na.rm = TRUE)

  # ==========================================================
  # Manual annotation table (same structure as your script)
  # ==========================================================

  sample_levels <- levels(clone_long$sampleName)

  manual_annot <- tibble::tibble(
    x_start = c(
      1 - 0.2,                 # Y: Intra
      2 - 0.2,                 # O: Intra
      1 - 0.2,                 # Intra Y vs O
      1 + 0.2                  # Inter Y vs O
    ),
    x_end = c(
      1 + 0.2,                 # Y: Inter
      2 + 0.2,                 # O: Inter
      2 - 0.2,                 # Intra Y vs O
      2 + 0.2                  # Inter Y vs O
    ),
    y_pos = c(
      y_max * 1.05,
      y_max * 1.05,
      y_max * 1.15,
      y_max * 1.25
    ),
    label = c(
      within_tests$label[1],
      within_tests$label[2],
      between_tests$label[between_tests$DistanceType == "Intra-clone"],
      between_tests$label[between_tests$DistanceType == "Inter-clone"]
    )
  )

  # ==========================================================
  # Plot
  # ==========================================================

  p <- ggplot2::ggplot(
    clone_long,
    ggplot2::aes(
      x = sampleName,
      y = Distance,
      fill = DistanceType
    )
  ) +
    ggplot2::geom_violin(
      trim = TRUE,
      alpha = 0.7,
      scale = "width",
      position = ggplot2::position_dodge(width = 0.8)
    ) +
    ggplot2::geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      position = ggplot2::position_dodge(width = 0.8)
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Intra-clone" = "cornflowerblue",
      "Inter-clone" = "orange3"
    )) +
    ggplot2::labs(
      x = "",
      y = "PCA distance",
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20, hjust = 0.5),
      axis.title.y = ggplot2::element_text(size = 21, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 18, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 16),
      legend.position = "top"
    ) +
    ggplot2::geom_segment(
      data = manual_annot,
      ggplot2::aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    ggplot2::geom_text(
      data = manual_annot,
      ggplot2::aes(
        x = (x_start + x_end) / 2,
        y = y_pos + y_max * 0.02,
        label = label
      ),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold"
    )

  list(
    summary = summary_df,
    pairwise = pairwise_df,
    long_format = clone_long,
    within_group_tests = within_tests,
    between_group_tests = between_tests,
    plot = p
  )
}
