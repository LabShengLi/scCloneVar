#' Plot variance, CV, and mean-adjusted variance summaries
#'
#' @importFrom magrittr %>%
#' @export
plot_variance_cv_summaries <- function(
    df,
    fdr_col = "fdr_brown_forsythe",
    fc_col = "log2FC_mean_adjusted_variance",
    fc_col_var = "log2FC_variance",
    fdr_cutoff = 0.05,
    fc_cutoff = 0.25,
    color_palette = c("#2F7BAA", "#E58B1C"),
    prefix = "Day0",
    group_order = NULL
) {
  
  # ---------------------------------------------------------
  # Detect group names from variance columns
  # ---------------------------------------------------------
  var_cols <- grep("^var_", colnames(df), value = TRUE)
  groups_detected <- gsub("^var_", "", var_cols)
  
  if (length(groups_detected) != 2)
    stop("Exactly two groups must be present.")
  
  if (is.null(group_order)) {
    g1 <- groups_detected[1]
    g2 <- groups_detected[2]
  } else {
    if (!all(group_order %in% groups_detected))
      stop("group_order must match detected groups: ",
           paste(groups_detected, collapse = ", "))
    g1 <- group_order[1]
    g2 <- group_order[2]
  }
  
  message("Plotting comparison: ", g2, " vs ", g1)
  
  get_wilcox_label <- function(p)
    if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
  
  # ==========================================================
  # Variance Boxplot
  # ==========================================================
  long_var <- df %>%
    dplyr::select(gene, dplyr::all_of(var_cols)) %>%
    tidyr::pivot_longer(
      -gene,
      names_to = "Group",
      values_to = "Variance"
    ) %>%
    dplyr::mutate(
      Group = gsub("^var_", "", Group),
      Group = factor(Group, levels = c(g1, g2))
    )
  
  p_wilcox <- stats::wilcox.test(Variance ~ Group, data = long_var)$p.value
  
  p_var <- ggplot2::ggplot(long_var,
                           ggplot2::aes(Group, Variance, color = Group)) +
    ggplot2::geom_boxplot(width = 0.5,
                          fill = "white",
                          outlier.shape = NA,
                          linewidth = 1) +
    ggplot2::geom_jitter(width = 0.15,
                         alpha = 0.5,
                         size = 1.6) +
    ggsignif::geom_signif(
      comparisons = list(c(g1, g2)),
      annotations = get_wilcox_label(p_wilcox),
      y_position = max(long_var$Variance, na.rm = TRUE) * 1.05
    ) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(
      title = sprintf("%s: Per-Gene Variance", prefix),
      subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes",
                         p_wilcox, nrow(df)),
      y = "Variance",
      x = NULL
    )
  
  # ==========================================================
  # CV Boxplot
  # ==========================================================
  cv_cols <- grep("^CV_", colnames(df), value = TRUE)
  
  long_cv <- df %>%
    dplyr::select(gene, dplyr::all_of(cv_cols)) %>%
    tidyr::pivot_longer(
      -gene,
      names_to = "Group",
      values_to = "CV"
    ) %>%
    dplyr::mutate(
      Group = gsub("^CV_", "", Group),
      Group = factor(Group, levels = c(g1, g2))
    )
  
  p_cv <- stats::wilcox.test(CV ~ Group, data = long_cv)$p.value
  
  p_cv_box <- ggplot2::ggplot(long_cv,
                              ggplot2::aes(Group, CV, color = Group)) +
    ggplot2::geom_boxplot(width = 0.5,
                          fill = "white",
                          outlier.shape = NA,
                          linewidth = 1) +
    ggplot2::geom_jitter(width = 0.15,
                         alpha = 0.5,
                         size = 1.6) +
    ggsignif::geom_signif(
      comparisons = list(c(g1, g2)),
      annotations = get_wilcox_label(p_cv),
      y_position = max(long_cv$CV, na.rm = TRUE) * 1.05
    ) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(
      title = sprintf("%s: Coefficient of Variation", prefix),
      subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes",
                         p_cv, nrow(df)),
      y = "Coefficient of Variation",
      x = NULL
    )
  
  # ==========================================================
  # Mean-adjusted Variance Boxplot
  # ==========================================================
  mvar_cols <- grep("^mean_adjusted_var_", colnames(df), value = TRUE)
  
  long_mvar <- df %>%
    dplyr::select(gene, dplyr::all_of(mvar_cols)) %>%
    tidyr::pivot_longer(
      -gene,
      names_to = "Group",
      values_to = "MeanAdjustedVar"
    ) %>%
    dplyr::mutate(
      Group = gsub("^mean_adjusted_var_", "", Group),
      Group = factor(Group, levels = c(g1, g2))
    )
  
  p_mv <- stats::wilcox.test(MeanAdjustedVar ~ Group,
                             data = long_mvar)$p.value
  
  p_mv_box <- ggplot2::ggplot(long_mvar,
                              ggplot2::aes(Group,
                                           MeanAdjustedVar,
                                           color = Group)) +
    ggplot2::geom_boxplot(width = 0.5,
                          fill = "white",
                          outlier.shape = NA,
                          linewidth = 1) +
    ggplot2::geom_jitter(width = 0.15,
                         alpha = 0.5,
                         size = 1.6) +
    ggsignif::geom_signif(
      comparisons = list(c(g1, g2)),
      annotations = get_wilcox_label(p_mv),
      y_position = max(long_mvar$MeanAdjustedVar,
                       na.rm = TRUE) * 1.05
    ) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(
      title = sprintf("%s: Mean-Adjusted Variance", prefix),
      subtitle = sprintf("Wilcoxon p = %.2e | n = %d genes",
                         p_mv, nrow(df)),
      y = "Mean-Adjusted Variance",
      x = NULL
    )
  
  # ==========================================================
  # Density Plot
  # ==========================================================
  p_density <- ggplot2::ggplot(
    long_mvar,
    ggplot2::aes(x = MeanAdjustedVar, color = Group)
  ) +
    ggplot2::geom_density(linewidth = 1.2, adjust = 1.1) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(
      title = sprintf("%s: Distribution of Mean-Adjusted Variance",
                      prefix),
      x = "Mean-Adjusted Variance",
      y = "Density"
    )
  
  # ==========================================================
  # Summary Counts
  # ==========================================================
  fc_thresholds <- c(0, 0.1, 0.25, 0.5)
  
  summary_counts_mvar <- purrr::map_dfr(
    fc_thresholds,
    function(fc_thr) {
      df %>%
        dplyr::filter(.data[[fdr_col]] < 0.5) %>%
        dplyr::summarise(
          log2FC_cutoff = fc_thr,
          up_in_g2 = sum(.data[[fc_col]] > fc_thr,
                         na.rm = TRUE),
          up_in_g1 = sum(.data[[fc_col]] < -fc_thr,
                         na.rm = TRUE)
        ) %>%
        dplyr::mutate(metric = "Mean-Adjusted Variance")
    }
  )
  
  summary_counts_var <- purrr::map_dfr(
    fc_thresholds,
    function(fc_thr) {
      df %>%
        dplyr::filter(.data[[fdr_col]] < 0.5) %>%
        dplyr::summarise(
          log2FC_cutoff = fc_thr,
          up_in_g2 = sum(.data[[fc_col_var]] > fc_thr,
                         na.rm = TRUE),
          up_in_g1 = sum(.data[[fc_col_var]] < -fc_thr,
                         na.rm = TRUE)
        ) %>%
        dplyr::mutate(metric = "Raw Variance")
    }
  )
  
  summary_counts <- dplyr::bind_rows(
    summary_counts_mvar,
    summary_counts_var
  )
  
  # ==========================================================
  # Volcano Plots
  # ==========================================================
  df_vol_mvar <- df %>%
    dplyr::mutate(
      neg_log10_fdr =
        -log10(pmax(.data[[fdr_col]], 1e-300)),
      significance = dplyr::case_when(
        .data[[fdr_col]] < fdr_cutoff &
          .data[[fc_col]] > fc_cutoff ~ paste0("↑", g2),
        .data[[fdr_col]] < fdr_cutoff &
          .data[[fc_col]] < -fc_cutoff ~ paste0("↑", g1),
        TRUE ~ "Not significant"
      )
    )
  
  top_up_g2 <- df_vol_mvar %>%
    dplyr::filter(significance == paste0("↑", g2)) %>%
    dplyr::arrange(desc(.data[[fc_col]])) %>%
    dplyr::slice_head(n = 10)
  
  top_up_g1 <- df_vol_mvar %>%
    dplyr::filter(significance == paste0("↑", g1)) %>%
    dplyr::arrange(.data[[fc_col]]) %>%
    dplyr::slice_head(n = 10)
  
  top_genes_mvar <- dplyr::bind_rows(top_up_g2, top_up_g1)
  
  volcano_colors <- c("#5271AE", "#D85B59", "grey80")
  names(volcano_colors) <- c(
    paste0("↑", g1),
    paste0("↑", g2),
    "Not significant"
  )
  
  p_volcano_mvar <- ggplot2::ggplot(
    df_vol_mvar,
    ggplot2::aes(x = .data[[fc_col]],
                 y = neg_log10_fdr)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = significance),
      alpha = 0.85,
      size = 2.4
    ) +
    ggrepel::geom_text_repel(
      data = top_genes_mvar,
      ggplot2::aes(label = gene),
      size = 5,
      color = "black",
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.4,
      segment.size = 0.3
    ) +
    ggplot2::scale_color_manual(values = volcano_colors) +
    ggplot2::geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed"
    ) +
    ggplot2::xlim(-3.5, 3.5) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = paste0("log2 Mean-Adjusted Variance (",
                 g2, "/", g1, ")"),
      y = expression(-log[10]("FDR (Brown–Forsythe)"))
    )
  
  # Raw variance volcano
  df_vol_var <- df %>%
    dplyr::mutate(
      neg_log10_fdr =
        -log10(pmax(.data[[fdr_col]], 1e-300)),
      significance_var = dplyr::case_when(
        .data[[fdr_col]] < fdr_cutoff &
          .data[[fc_col_var]] > fc_cutoff ~ paste0("↑", g2),
        .data[[fdr_col]] < fdr_cutoff &
          .data[[fc_col_var]] < -fc_cutoff ~ paste0("↑", g1),
        TRUE ~ "Not significant"
      )
    )
  
  top_up_g2_var <- df_vol_var %>%
    dplyr::filter(significance_var ==
                    paste0("↑", g2)) %>%
    dplyr::arrange(desc(.data[[fc_col_var]])) %>%
    dplyr::slice_head(n = 10)
  
  top_up_g1_var <- df_vol_var %>%
    dplyr::filter(significance_var ==
                    paste0("↑", g1)) %>%
    dplyr::arrange(.data[[fc_col_var]]) %>%
    dplyr::slice_head(n = 10)
  
  top_genes_var <- dplyr::bind_rows(
    top_up_g2_var,
    top_up_g1_var
  )
  
  p_volcano_var <- ggplot2::ggplot(
    df_vol_var,
    ggplot2::aes(x = .data[[fc_col_var]],
                 y = neg_log10_fdr)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = significance_var),
      alpha = 0.85,
      size = 2.4
    ) +
    ggrepel::geom_text_repel(
      data = top_genes_var,
      ggplot2::aes(label = gene),
      size = 5,
      color = "black",
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.4,
      segment.size = 0.3
    ) +
    ggplot2::scale_color_manual(values = volcano_colors) +
    ggplot2::geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed"
    ) +
    ggplot2::xlim(-3.5, 3.5) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = paste0("log2 Variance (",
                 g2, "/", g1, ")"),
      y = expression(-log[10]("FDR (Brown–Forsythe)"))
    )
  
  # ==========================================================
  # Return
  # ==========================================================
  list(
    variance_box = p_var,
    cv_box = p_cv_box,
    mean_adj_var_box = p_mv_box,
    density_mean_adj = p_density,
    volcano_mean_adj = p_volcano_mvar,
    volcano_variance = p_volcano_var,
    summary_counts = summary_counts
  )
}