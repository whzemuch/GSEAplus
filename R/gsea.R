#' GSEA core utilities
#'
#' Functions for building ranked gene lists, running GSEA across contrasts,
#' retrieving specific results, and generating GSEA enrichment plots with
#' embedded statistics.

make_geneList <- function(dt, contrast_name) {
  #' make a geneList from result_df from DESeq2
  #' @export
  #'
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  dt[
    contrast == contrast_name & is.finite(stat),
    .SD[!duplicated(gene)],
    by = gene
  ][
    order(-stat),
    setNames(stat, gene)
  ]
}

make_geneLists <- function(res_df) {
  #' make a list of geneList from combined result_df from DESeq2
  #' @export
  #'

  contrasts <- unique(res_df$contrast)
  lists <- lapply(contrasts, function(ct) make_geneList(res_df, ct))
  names(lists) <- contrasts
  lists
}

get_gsea <- function(gsea_lst, cell, contrast) {
  #' extract gsea results for plotting
  #' @export
  #'
  if (!cell %in% names(gsea_lst)) {
    stop(paste("Cell type not found:", cell))
  }
  if (!contrast %in% names(gsea_lst[[cell]])) {
    stop(paste("Contrast not found:", contrast))
  }
  gsea_lst[[cell]][[contrast]]
}

gseaplot_with_stats <- function(
  gsea_obj,
  geneSetID,
  title = NULL,
  text_size = 6
) {
  #' Run gseaplots and then add annotations
  #' @export
  #'
  if (is.null(title)) {
    title <- geneSetID
  }

  stats <- gsea_obj@result |>
    dplyr::filter(ID == geneSetID) |>
    dplyr::select(NES, qvalue)

  if (nrow(stats) == 0) {
    stop(paste("Gene set", geneSetID, "not found"))
  }

  p <- enrichplot::gseaplot2(gsea_obj, geneSetID, title = title)

  es_plot <- p[[1]] +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      label = sprintf("NES = %.2f\nqval = %.2e", stats$NES, stats$qvalue),
      hjust = 1.05,
      vjust = 1.1,
      size = text_size
    )

  p[[1]] <- es_plot
  p
}
