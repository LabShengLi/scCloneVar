#' Wrapper function that produce distribution plots and statistical summary
#'
#' @export
run_clone_distribution_engine <- function(
    comparison_list,
    output_file = NULL,
    rep_col = "Rep",
    clone_col = "CloneID"
) {

  # --------------------------
  # Palette builder
  # --------------------------
  build_palette <- function(levels_vec, other_lab) {

    pal_raw <- c(
      RColorBrewer::brewer.pal(9, "Set1"),
      RColorBrewer::brewer.pal(8, "Dark2"),
      RColorBrewer::brewer.pal(8, "Accent"),
      ggsci::pal_npg("nrc")(10),
      ggsci::pal_d3("category20")(20)
    )

    hsv_vals <- colorspace::coords(
      methods::as(colorspace::hex2RGB(pal_raw), "HSV")
    )

    pal_raw <- pal_raw[hsv_vals[, "S"] > 0.4]

    pal <- rep(pal_raw, length.out = length(levels_vec))
    names(pal) <- levels_vec

    pal[other_lab] <- "grey80"

    pal
  }

  # --------------------------
  # Stats computation
  # --------------------------
  compute_stats <- function(meta_df) {

    meta_df |>
      dplyr::filter(.data[[clone_col]] != "0") |>
      dplyr::group_by(.data$sampleName, .data[[clone_col]]) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(.data$sampleName) |>
      dplyr::mutate(freq = n / sum(n)) |>
      dplyr::summarise(
        n_clones = dplyr::n(),
        min_freq = base::min(freq),
        max_freq = base::max(freq),
        mean_freq = base::mean(freq),
        .groups = "drop"
      )
  }

  plot_list <- list()
  stats_list <- list()

  # ======================================================
  # LOOP THROUGH USER-DEFINED COMPARISONS
  # ======================================================
  for (comp in comparison_list) {

    seurat_obj   <- comp$seurat_obj
    sample_names <- comp$samples
    label_counts <- comp$label_counts
    threshold    <- comp$threshold
    title        <- comp$title
    comp_name    <- comp$name

    other_lab <- paste0("Other (<", threshold * 100, "%)")

    meta_df <- seurat_obj@meta.data |>
      dplyr::filter(
        .data$sampleName %in% sample_names,
        .data[[rep_col]] != "UnMapped",
        .data[[clone_col]] != "0"
      ) |>
      dplyr::mutate(
        Sample = .data$sampleName,
        Replicate = .data[[rep_col]],
        CloneID = as.character(.data[[clone_col]])
      )

    clone_freq <- meta_df |>
      dplyr::count(Sample, Replicate, CloneID, name = "n_cells") |>
      dplyr::group_by(Sample, Replicate) |>
      dplyr::mutate(RelFreq = n_cells / sum(n_cells)) |>
      dplyr::ungroup() |>
      dplyr::mutate(CloneLabel = ifelse(RelFreq < threshold, other_lab, CloneID)) |>
      dplyr::group_by(Sample, Replicate, CloneLabel) |>
      dplyr::summarise(RelFreq = sum(RelFreq), .groups = "drop") |>
      dplyr::left_join(label_counts, by = c("Replicate","Sample")) |>
      dplyr::mutate(FacetOrder = paste(Sample, Replicate, sep = "_"))

    clone_order <- clone_freq |>
      dplyr::group_by(CloneLabel) |>
      dplyr::summarise(total = sum(RelFreq)) |>
      dplyr::arrange(dplyr::desc(total)) |>
      dplyr::pull(CloneLabel)

    clone_freq$CloneLabel <- factor(
      clone_freq$CloneLabel,
      levels = c(other_lab, setdiff(clone_order, other_lab))
    )

    pal <- build_palette(levels(clone_freq$CloneLabel), other_lab)

    p <- ggplot2::ggplot(
      clone_freq,
      ggplot2::aes(
        x = SampleLabel,
        y = RelFreq,
        fill = CloneLabel
      )
    ) +
      ggplot2::geom_bar(
        stat = "identity",
        width = 0.9,
        position = ggplot2::position_stack(reverse = TRUE)
      ) +
      ggplot2::facet_wrap(~FacetOrder, nrow = 1) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        title = title,
        y = "Relative clone size",
        x = NULL
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none")

    plot_list[[comp_name]] <- p

    stats_list[[comp_name]] <- compute_stats(seurat_obj@meta.data)
  }

  # --------------------------
  # Combine all panels
  # --------------------------
  combined_panel <- patchwork::wrap_plots(plot_list)

  if (!is.null(output_file)) {
    ggplot2::ggsave(
      output_file,
      combined_panel,
      width = 18,
      height = 10,
      device = grDevices::cairo_pdf
    )
  }

  list(
    plots = plot_list,
    combined_plot = combined_panel,
    descriptive_stats = stats_list
  )
}
