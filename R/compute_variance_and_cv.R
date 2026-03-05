#' Compute variance tests, CV, and mean-adjusted variance metrics
#'
#' @importFrom magrittr %>%
#' @export
compute_variance_and_cv <- function(
    seu,
    assay = "RNA",
    slot = "data",
    group_col = "donor_age",
    n_hvg = 2000,
    selection.method = "vst",
    n_cores = 1,
    chunk_size = 250,
    tests = c("brown_forsythe", "levene", "bartlett"),
    min_pct = 0.1,
    min_mean = 0.1
) {

  stopifnot(inherits(seu, "Seurat"))
  stopifnot(group_col %in% colnames(seu@meta.data))

  groups_raw <- seu[[group_col]][, 1]
  groups_raw <- groups_raw[!is.na(groups_raw)]

  if (is.factor(groups_raw)) {
    groups <- levels(droplevels(groups_raw))
  } else {
    groups <- unique(as.character(groups_raw))
  }

  if (length(groups) != 2) {
    stop("group_col must contain exactly TWO groups.")
  }

  g1 <- groups[1]
  g2 <- groups[2]

  message("Comparing groups: ", g1, " vs ", g2)

  # ---------------- HVG selection ----------------
  hvgs <- Seurat::VariableFeatures(seu)

  if (length(hvgs) == 0) {
    seu <- Seurat::FindVariableFeatures(
      seu,
      assay = assay,
      selection.method = selection.method,
      nfeatures = n_hvg
    )
    hvgs <- Seurat::VariableFeatures(seu)
  } else if (length(hvgs) > n_hvg) {
    hvgs <- hvgs[1:n_hvg]
  }

  expr <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = slot),
    error = function(e)
      Seurat::GetAssayData(seu, assay = assay, slot = slot)
  )

  expr <- expr[intersect(rownames(expr), hvgs), , drop = FALSE]

  group_factor <- factor(seu[[group_col]][, 1], levels = groups)

  detect_pct <- Matrix::rowMeans(expr > 0)
  mean_expr  <- Matrix::rowMeans(expr)

  keep_genes <- names(which(detect_pct >= min_pct & mean_expr >= min_mean))

  expr <- expr[keep_genes, , drop = FALSE]

  message("Keeping ", length(keep_genes), " genes after filtering.")

  # ---------------- Variance tests ----------------
  do_tests <- function(x, group_factor) {

    res <- list(bf = NA, lev = NA, bart = NA)

    df <- data.frame(expr = x, grp = group_factor)

    if (length(unique(df$grp)) < 2) return(res)

    # Brown–Forsythe
    med <- tapply(df$expr, df$grp, median)
    dev <- abs(df$expr - med[df$grp])
    fit <- try(stats::aov(dev ~ grp, data = df), silent = TRUE)

    if (!inherits(fit, "try-error"))
      res$bf <- summary(fit)[[1]]["grp", "Pr(>F)"]

    # Levene (mean-centered)
    mu <- tapply(df$expr, df$grp, mean)
    dev <- abs(df$expr - mu[df$grp])
    fit <- try(stats::aov(dev ~ grp, data = df), silent = TRUE)

    if (!inherits(fit, "try-error"))
      res$lev <- summary(fit)[[1]]["grp", "Pr(>F)"]

    bt <- try(stats::bartlett.test(expr ~ grp, data = df), silent = TRUE)

    if (!inherits(bt, "try-error"))
      res$bart <- bt$p.value

    res
  }

  message("Running variance tests ...")

  progressr::handlers("txtprogressbar")

  chunks <- split(
    rownames(expr),
    ceiling(seq_along(rownames(expr)) / chunk_size)
  )

  results_list <- list()

  progressr::with_progress({

    p <- progressr::progressor(steps = length(chunks))

    for (i in seq_along(chunks)) {

      p(sprintf("Chunk %d / %d", i, length(chunks)))

      genes <- chunks[[i]]
      sub_expr <- expr[genes, , drop = FALSE]

      mat <- apply(sub_expr, 1, do_tests, group_factor = group_factor)

      bf_p   <- sapply(mat, `[[`, "bf")
      lev_p  <- sapply(mat, `[[`, "lev")
      bart_p <- sapply(mat, `[[`, "bart")

      var_tbl <- apply(sub_expr, 1, function(x) {

        df <- data.frame(expr = x, grp = group_factor)

        var_g1 <- stats::var(df$expr[df$grp == g1], na.rm = TRUE)
        var_g2 <- stats::var(df$expr[df$grp == g2], na.rm = TRUE)

        log2FC_variance <- ifelse(
          var_g1 > 0,
          log2(var_g2 / var_g1),
          NA_real_
        )

        c(var_g1, var_g2, log2FC_variance)
      })

      var_tbl <- as.data.frame(t(var_tbl))

      colnames(var_tbl) <- c(
        paste0("var_", g1),
        paste0("var_", g2),
        "log2FC_variance"
      )

      results_list[[i]] <- data.frame(
        gene = genes,
        p_brown_forsythe = bf_p,
        p_levene = lev_p,
        p_bartlett = bart_p,
        var_tbl
      )
    }
  })

  var_results <- dplyr::bind_rows(results_list) %>%
    dplyr::mutate(
      fdr_brown_forsythe = stats::p.adjust(p_brown_forsythe, "BH"),
      fdr_levene = stats::p.adjust(p_levene, "BH"),
      fdr_bartlett = stats::p.adjust(p_bartlett, "BH")
    )

  # ---------------- Mean-adjusted variance ----------------
  message("Computing mean-adjusted variance and SD ...")

  mean_g1 <- apply(expr[, group_factor == g1], 1, mean, na.rm = TRUE)
  mean_g2 <- apply(expr[, group_factor == g2], 1, mean, na.rm = TRUE)

  mean_var_df <- data.frame(
    gene = rownames(expr),
    group = rep(c(g1, g2), each = nrow(expr)),
    mean = c(mean_g1, mean_g2),
    var = c(
      var_results[[paste0("var_", g1)]],
      var_results[[paste0("var_", g2)]]
    )
  )

  mean_var_df <- mean_var_df %>%
    dplyr::filter(!is.na(mean) & mean > 0 & !is.na(var) & var > 0)

  fit_loess <- stats::loess(
    log10(var) ~ log10(mean),
    data = mean_var_df,
    span = 0.75
  )

  mean_var_df$expected_logVar <-
    stats::predict(fit_loess, newdata = mean_var_df)

  mean_var_df$residual_logVar <-
    log10(mean_var_df$var) - mean_var_df$expected_logVar

  mean_var_df$mean_adjusted_var <- 10 ^ mean_var_df$residual_logVar
  mean_var_df$mean_adjusted_sd  <- sqrt(mean_var_df$mean_adjusted_var)

  adj_tbl <- mean_var_df %>%
    dplyr::select(gene, group, mean_adjusted_var, mean_adjusted_sd) %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = c(mean_adjusted_var, mean_adjusted_sd)
    ) %>%
    dplyr::mutate(
      log2FC_mean_adjusted_variance =
        log2((.data[[paste0("mean_adjusted_var_", g2)]] + 1e-8) /
               (.data[[paste0("mean_adjusted_var_", g1)]] + 1e-8)),
      log2FC_mean_adjusted_SD =
        log2((.data[[paste0("mean_adjusted_sd_", g2)]] + 1e-8) /
               (.data[[paste0("mean_adjusted_sd_", g1)]] + 1e-8))
    )

  # ---------------- CV ----------------
  expr_full <- Seurat::GetAssayData(seu, assay = assay, slot = slot)

  pooled_mean <- Matrix::rowMeans(expr_full, na.rm = TRUE)

  cv_list <- purrr::map(groups, function(g) {

    cells <- colnames(seu)[seu[[group_col]][, 1] == g]

    mat <- expr_full[, cells, drop = FALSE]

    sd_vals <- apply(mat, 1, stats::sd, na.rm = TRUE)

    sd_vals / (pooled_mean + 1e-8)
  })

  cv_tbl <- as.data.frame(cv_list)

  colnames(cv_tbl) <- paste0("CV_", groups)

  cv_tbl$gene <- rownames(expr_full)

  cv_tbl <- cv_tbl %>%
    dplyr::mutate(
      log2_CV_ratio =
        log2((.data[[paste0("CV_", g2)]] + 1e-8) /
               (.data[[paste0("CV_", g1)]] + 1e-8))
    )

  merged <- var_results %>%
    dplyr::inner_join(cv_tbl, by = "gene") %>%
    dplyr::left_join(adj_tbl, by = "gene")

  message("Done. Generated table with variance tests, CV, and mean-adjusted variance metrics.")

  return(merged)
}
