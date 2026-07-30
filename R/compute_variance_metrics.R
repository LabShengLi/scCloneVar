#' Compute variance tests and mean-adjusted variance metrics
#'
#' @importFrom magrittr %>%
#' @export
compute_variance_metrics <- function(
    seu,
    assay = "RNA",
    layer = "data",
    group_col = "donor_age",
    n_hvg = 2000,
    selection.method = "vst",
    n_cores = 1,
    chunk_size = 250,
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
    Seurat::GetAssayData(seu, assay = assay, layer = layer),
    error = function(e)
      Seurat::GetAssayData(seu, assay = assay, slot = layer)
  )

  expr <- expr[intersect(rownames(expr), hvgs), , drop = FALSE]

  group_factor <- factor(seu[[group_col]][, 1], levels = groups)

  detect_pct <- Matrix::rowMeans(expr > 0)
  mean_expr  <- Matrix::rowMeans(expr)

  keep_genes <- names(which(detect_pct >= min_pct & mean_expr >= min_mean))

  expr <- expr[keep_genes, , drop = FALSE]

  message("Keeping ", length(keep_genes), " genes after filtering.")

  idx_g1 <- which(group_factor == g1)
  idx_g2 <- which(group_factor == g2)

  # ---------------- Mean-variance trend ----------------
  mean_g1 <- Matrix::rowMeans(expr[, idx_g1, drop = FALSE])
  mean_g2 <- Matrix::rowMeans(expr[, idx_g2, drop = FALSE])
  var_g1_raw <- apply(expr[, idx_g1, drop = FALSE], 1, stats::var)
  var_g2_raw <- apply(expr[, idx_g2, drop = FALSE], 1, stats::var)

  trend_df <- data.frame(
    mean = c(mean_g1, mean_g2),
    var  = c(var_g1_raw, var_g2_raw)
  )
  trend_df <- trend_df[
    !is.na(trend_df$mean) & trend_df$mean > 0 &
      !is.na(trend_df$var)  & trend_df$var  > 0,
  ]

  fit_loess <- stats::loess(
    log10(var) ~ log10(mean),
    data = trend_df,
    span = 0.75,
    control = stats::loess.control(surface = "direct")
  )

  expected_var <- function(mean_vec) {
    out <- rep(NA_real_, length(mean_vec))
    ok <- !is.na(mean_vec) & mean_vec > 0
    out[ok] <- 10 ^ stats::predict(fit_loess, newdata = data.frame(mean = mean_vec[ok]))
    out
  }

  scale_g1 <- sqrt(expected_var(mean_g1))
  scale_g2 <- sqrt(expected_var(mean_g2))

  testable <- !is.na(scale_g1) & !is.na(scale_g2) & scale_g1 > 0 & scale_g2 > 0

  adj_expr <- expr
  adj_expr[, idx_g1] <- expr[, idx_g1, drop = FALSE] / scale_g1
  adj_expr[, idx_g2] <- expr[, idx_g2, drop = FALSE] / scale_g2

  message(sum(!testable), " gene(s) excluded: mean-variance trend not evaluable.")

  # ---------------- Variance tests ----------------
  do_tests <- function(x, group_factor) {

    res <- list(bf = NA, lev = NA)

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

    res
  }

  message("Running variance tests ...")

  progressr::handlers("txtprogressbar")

  test_genes <- rownames(expr)[testable]

  chunks <- split(
    test_genes,
    ceiling(seq_along(test_genes) / chunk_size)
  )

  results_list <- list()

  progressr::with_progress({

    p <- progressr::progressor(steps = length(chunks))

    for (i in seq_along(chunks)) {

      p(sprintf("Chunk %d / %d", i, length(chunks)))

      genes <- chunks[[i]]
      sub_adj <- adj_expr[genes, , drop = FALSE]
      sub_raw <- expr[genes, , drop = FALSE]

      mat <- apply(sub_adj, 1, do_tests, group_factor = group_factor)

      bf_p   <- sapply(mat, `[[`, "bf")
      lev_p  <- sapply(mat, `[[`, "lev")

      var_tbl <- apply(sub_raw, 1, function(x) {

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
        var_tbl
      )
    }
  })

  var_results <- dplyr::bind_rows(results_list) %>%
    dplyr::mutate(
      fdr_brown_forsythe = stats::p.adjust(p_brown_forsythe, "BH"),
      fdr_levene = stats::p.adjust(p_levene, "BH")
    )

  # ---------------- Mean-adjusted variance ----------------
  message("Computing mean-adjusted variance and SD ...")

  mean_adjusted_var_g1 <- apply(adj_expr[, idx_g1, drop = FALSE], 1, stats::var, na.rm = TRUE)
  mean_adjusted_var_g2 <- apply(adj_expr[, idx_g2, drop = FALSE], 1, stats::var, na.rm = TRUE)

  adj_tbl <- data.frame(gene = rownames(expr))
  adj_tbl[[paste0("mean_adjusted_var_", g1)]] <- mean_adjusted_var_g1
  adj_tbl[[paste0("mean_adjusted_var_", g2)]] <- mean_adjusted_var_g2
  adj_tbl[[paste0("mean_adjusted_sd_", g1)]]  <- sqrt(mean_adjusted_var_g1)
  adj_tbl[[paste0("mean_adjusted_sd_", g2)]]  <- sqrt(mean_adjusted_var_g2)

  adj_tbl <- adj_tbl %>%
    dplyr::mutate(
      log2FC_mean_adjusted_variance =
        log2((.data[[paste0("mean_adjusted_var_", g2)]] + 1e-8) /
               (.data[[paste0("mean_adjusted_var_", g1)]] + 1e-8)),
      log2FC_mean_adjusted_SD =
        log2((.data[[paste0("mean_adjusted_sd_", g2)]] + 1e-8) /
               (.data[[paste0("mean_adjusted_sd_", g1)]] + 1e-8))
    )

  merged <- var_results %>% dplyr::left_join(adj_tbl, by = "gene")

  message("Done. Generated table with variance tests and mean-adjusted variance metrics.")

  return(merged)
}
