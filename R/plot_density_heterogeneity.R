#' Plot gene expression density between two groups
#'
#' Visualizes gene expression distributions and optional mean-adjusted variance
#' between two groups of cells in a Seurat object.
#'
#' @param seurat_obj Seurat object
#' @param gene Gene name
#' @param group_col Metadata column defining groups
#' @param group1_labels Labels belonging to group1
#' @param group2_labels Labels belonging to group2
#' @param group1_name Display name of group1
#' @param group2_name Display name of group2
#' @param title_prefix Plot title prefix
#' @param dvgs_df Dataframe containing mean-adjusted variance
#' @param mean_adj_var_group1_col Column for group1 variance
#' @param mean_adj_var_group2_col Column for group2 variance
#' @param group_colors Color vector for the two groups
#'
#' @return ggplot object
#' @export
plot_density_gene <- function(
    seurat_obj,
    gene,
    group_col,
    group1_labels,
    group2_labels,
    group1_name = "Group1",
    group2_name = "Group2",
    title_prefix = "",
    dvgs_df = NULL,
    mean_adj_var_group1_col = NULL,
    mean_adj_var_group2_col = NULL,
    group_colors = c("#89AEEB", "#F5A36C")
){
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
  })
  # check if the gene exist
  if(!(gene %in% rownames(seurat_obj))){
    message("Skipping ", gene, " (not found in Seurat object)")
    return(NULL)
  }
  
  # Fetch expression + metadata
  expr_df <- Seurat::FetchData(seurat_obj, vars = c(gene, group_col)) %>%
    dplyr::mutate(GroupRaw = .data[[group_col]]) %>%
    dplyr::mutate(
      Group = dplyr::case_when(
        GroupRaw %in% group1_labels ~ group1_name,
        GroupRaw %in% group2_labels ~ group2_name,
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Group)) %>%
    dplyr::mutate(Group = factor(Group, levels = c(group1_name, group2_name)),
                  Expression = .data[[gene]]) %>%
    dplyr::filter(Expression > 0)
  color_palette <- c(
    setNames(group_colors[1], group1_name),
    setNames(group_colors[2], group2_name)
  )
  # Density plot
  p_density <- ggplot(expr_df, aes(x = Expression, color = Group, fill = Group)) +
    geom_density(alpha = 0.25, linewidth = 1.1) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    theme_classic(base_size = 15) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    ) +
    labs(
      title = paste0(title_prefix, " — ", gene),
      x = "Normalized Expression",
      y = "Density"
    )
  # Variance bar plot
  if(!is.null(dvgs_df) &&
     gene %in% dvgs_df$gene &&
     !is.null(mean_adj_var_group1_col) &&
     !is.null(mean_adj_var_group2_col)){
    var_vals <- dvgs_df %>%
      dplyr::filter(.data$gene == gene)
    if(nrow(var_vals) == 0){
      message("Gene ", gene, " not found in dvgs_df")
      return(p_density)
    }
    mean_var_df <- tibble::tibble(
      Group = factor(c(group1_name, group2_name),
                     levels = c(group1_name, group2_name)),
      MeanAdjVar = c(
        var_vals[[mean_adj_var_group1_col]][1],
        var_vals[[mean_adj_var_group2_col]][1]
      )
    )
    p_bar <- ggplot(mean_var_df,
                    aes(x = Group, y = MeanAdjVar, fill = Group)) +
      geom_bar(stat = "identity", width = 0.6, color = "black") +
      scale_fill_manual(values = color_palette) +
      theme_classic(base_size = 15) +
      theme(
        axis.title = element_text(face = "bold"),
        legend.position = "none"
      ) +
      labs(
        title = "Mean-Adjusted Variance",
        x = NULL,
        y = "Mean-Adj Var"
      ) +
      scale_y_continuous(expand = c(0,0), limits = c(0,NA))
    p_final <- p_density + p_bar + patchwork::plot_layout(widths = c(2.2,1))
  } else {
    p_final <- p_density
  }
  return(p_final)
}