#' Run MSigDB GSEA across multiple collections
#'
#' @importFrom magrittr %>%
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @export
run_msigdb_gsea_all_collections <- function(
    hsc_results,
    species = "Mus musculus",
    msig_collections = c(
      "H",
      "C2:CGP", "C2:BIOCARTA", "C2:KEGG_MEDICUS", "C2:PID", "C2:REACTOME",
      "C2:WIKIPATHWAYS", "C2:KEGG_LEGACY",
      "C3:MIRDB", "C3:MIR_LEGACY", "C3:GTRD", "C3:TFT_LEGACY",
      "C4:3CA", "C4:CGN", "C4:CM",
      "C5:BP", "C5:CC", "C5:MF", "C5:HPO",
      "C6",
      "C7:IMMUNESIGDB", "C7:VAX",
      "C8"
    ),
    rank_var = "log2FC_mean_adjusted_variance_YO",
    top_n = 6,
    title_prefix = "HSC"
) {

  # sanity checks
  stopifnot("gene" %in% colnames(hsc_results))
  stopifnot(rank_var %in% colnames(hsc_results))

  set.seed(42)

  # Select OrgDb based on species
  OrgDb <- switch(
    species,
    "Mus musculus" = org.Mm.eg.db::org.Mm.eg.db,
    "Homo sapiens" = org.Hs.eg.db::org.Hs.eg.db,
    stop("Unsupported species. Use 'Mus musculus' or 'Homo sapiens'.")
  )

  # SYMBOL → ENTREZ
  gene_map <- clusterProfiler::bitr(
    unique(hsc_results$gene),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = OrgDb
  )

  hsc_results <- hsc_results %>%
    dplyr::left_join(gene_map, by = c("gene" = "SYMBOL")) %>%
    dplyr::filter(!is.na(ENTREZID))

  # lookup table for later
  entrez2symbol <- gene_map %>%
    dplyr::select(ENTREZID, SYMBOL)

  # ranked gene list
  ranked_tbl <- hsc_results %>%
    dplyr::arrange(desc(.data[[rank_var]])) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)

  geneList <- ranked_tbl[[rank_var]]
  names(geneList) <- ranked_tbl$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)

  n_genes_used <- length(geneList)

  message("Input genes: ", length(unique(hsc_results$gene)))
  message("Mapped ENTREZ: ", length(unique(entrez2symbol$ENTREZID)))
  message("Final ranked genes: ", n_genes_used)

  # convert core_enrichment to SYMBOLs
  convert_core_to_symbol <- function(core_string) {
    if (is.na(core_string) || core_string == "") return(NA_character_)
    entrez_ids <- unlist(strsplit(core_string, "/"))
    symbols <- entrez2symbol %>%
      dplyr::filter(ENTREZID %in% entrez_ids) %>%
      dplyr::pull(SYMBOL) %>%
      unique()
    paste(symbols, collapse = "; ")
  }

  # run one MSigDB collection
  run_one_collection <- function(msig_collection) {

    message("\n MSigDB collection: ", msig_collection)

    if (grepl(":", msig_collection)) {
      parts <- strsplit(msig_collection, ":")[[1]]
      m_df <- msigdbr::msigdbr(
        species = species,
        collection = parts[1],
        subcollection = parts[2]
      )
    } else {
      m_df <- msigdbr::msigdbr(
        species = species,
        collection = msig_collection
      )
    }

    if (nrow(m_df) == 0) {
      message("No gene sets found.")
      return(NULL)
    }

    term2gene <- m_df %>%
      dplyr::select(gs_name, ncbi_gene)

    term2name <- m_df %>%
      dplyr::select(gs_name, gs_description)

    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList     = geneList,
        TERM2GENE    = term2gene,
        TERM2NAME    = term2name,
        pvalueCutoff = 1,
        verbose      = FALSE
      )
    }, error = function(e) NULL)

    if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
      return(list(
        gsea_result  = gsea_res,
        gsea_table   = NULL,
        panel_up     = NULL,
        panel_down   = NULL
      ))
    }

    # summary table with genes
    gsea_tbl <- gsea_res@result %>%
      dplyr::mutate(
        collection = msig_collection,
        rank_variable = rank_var,
        condition = title_prefix,
        n_genes_used = n_genes_used,
        leading_edge_genes = vapply(
          core_enrichment,
          convert_core_to_symbol,
          FUN.VALUE = character(1)
        )
      )

    # FDR-filtered plotting
    sig_res <- gsea_res@result %>%
      dplyr::filter(p.adjust < 0.05)

    top_up <- sig_res %>%
      dplyr::filter(NES > 0) %>%
      dplyr::arrange(desc(NES)) %>%
      dplyr::slice_head(n = top_n)

    top_down <- sig_res %>%
      dplyr::filter(NES < 0) %>%
      dplyr::arrange(NES) %>%
      dplyr::slice_head(n = top_n)

    make_panel <- function(df, label) {

      if (nrow(df) == 0) return(NULL)

      plots <- lapply(seq_len(nrow(df)), function(i) {
        tryCatch(
          ggplotify::as.ggplot(
            enrichplot::gseaplot2(
              gsea_res,
              geneSetID = df$ID[i],
              title = paste0(
                stringr::str_wrap(df$ID[i], 50),
                "\nFDR=", signif(df$p.adjust[i], 3),
                ", NES=", round(df$NES[i], 2)
              )
            )
          ),
          error = function(e) NULL
        )
      })

      plots <- plots[!sapply(plots, is.null)]
      if (length(plots) == 0) return(NULL)

      patchwork::wrap_plots(plots)
    }

    list(
      gsea_result = gsea_res,
      gsea_table  = gsea_tbl,
      panel_up    = make_panel(top_up, "Up"),
      panel_down  = make_panel(top_down, "Down")
    )
  }

  # run all collections
  results <- lapply(msig_collections, run_one_collection)
  names(results) <- msig_collections

  # combine plots
  up_panels <- lapply(results, `[[`, "panel_up")
  up_panels <- up_panels[!sapply(up_panels, is.null)]

  down_panels <- lapply(results, `[[`, "panel_down")
  down_panels <- down_panels[!sapply(down_panels, is.null)]

  combined_up_panel <- if (length(up_panels) > 0)
    patchwork::wrap_plots(up_panels) else NULL

  combined_down_panel <- if (length(down_panels) > 0)
    patchwork::wrap_plots(down_panels) else NULL

  # combine summary tables
  gsea_summary_table <- dplyr::bind_rows(
    lapply(results, `[[`, "gsea_table")
  ) %>%
    dplyr::arrange(collection, desc(NES))

  return(list(
    per_collection      = results,
    combined_up_panel   = combined_up_panel,
    combined_down_panel = combined_down_panel,
    gsea_summary_table  = gsea_summary_table
  ))
}
