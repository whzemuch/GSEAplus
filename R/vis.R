#' GSEAplus visualization workflow
#'
#' This file contains utilities for converting GSEA results into matrices
#' and producing publication-quality heatmaps.
#'
#' A typical end-to-end workflow:
#'
#' \dontrun{
#' library(GSEAplus)
#' library(SummarizedExperiment)
#'
#' # Load pseudobulk count matrices (from Scanpy → H5AD → SummarizedExperiment)
#' myo_se <- readH5AD("h5d/age_adata_myo_counts_20260110s.h5ad")
#' av_se  <- readH5AD("h5d/age_adata_AV_counts_2026012.h5ad")
#' hs_se  <- readH5AD("h5d/age_adata_HS_counts_2026012.h5ad")
#'
#' # Run DESeq2 pseudobulk
#' se_lst  <- list(myo = myo_se, av = av_se, hs = hs_se)
#' out_lst <- lapply(se_lst, run_pseudobulk_deseq)
#'
#' # Run multi-contrast GSEA
#' gsea_lst <- lapply(out_lst, function(x) {
#'   compareCluster(
#'     gene | stat ~ contrast,
#'     data = x$res_df,
#'     fun = "GSEA",
#'     TERM2GENE = H2T,
#'     eps = 0,
#'     pvalueCutoff = 1,
#'     minGSSize = 10,
#'     maxGSSize = 5000
#'   )
#' })
#'
#' # Convert to NES and q-value matrices
#' gsea_mats_lst <- lapply(gsea_lst, make_gsea_matrices_dt)
#'
#' # Build SummarizedExperiment for plotting
#' se_plot_lst <- lapply(gsea_mats_lst, function(mat) {
#'   SummarizedExperiment(
#'     assays = list(
#'       NES  = mat$NES,
#'       qval = mat$qvalue
#'     ),
#'     colData = S4Vectors::DataFrame(
#'       comparison = colnames(mat$NES),
#'       row.names  = colnames(mat$NES)
#'     )
#'   )
#' })
#'
#' # Extract matrices for heatmap
#' plot_mats_lst <- lapply(se_plot_lst, extract_plot_mats)
#'
#' # Plot one cell type
#' cell_type <- "myo"
#' resizeFig(15, 10)
#' draw_heatmap(
#'   plot_mats_lst[[cell_type]]$NES,
#'   plot_mats_lst[[cell_type]]$qval,
#'   plot_mats_lst[[cell_type]]$col_df
#' )
#' }

resizeFig <- function(width = 10, height = 5) {
  #' Resize Jupyter plot output
  #' @export
  #'
  options(repr.plot.width = width, repr.plot.height = height)
}

extract_plot_mats <- function(se) {
  #' get NES, qval and col_df for plotting heatmap
  #' @export
  #'
  col_df <- as.data.frame(SummarizedExperiment::colData(se))
  col_df$comparison <- rownames(col_df)

  list(
    col_df = col_df,
    NES = SummarizedExperiment::assay(se, "NES"),
    qval = SummarizedExperiment::assay(se, "qval")
  )
}

make_gsea_matrices_dt <- function(
  gsea_df,
  value_cols = c("NES", "qvalue"),
  row_col = "ID",
  col_col = "Cluster"
) {
  #' make a list mats for constructing the plotting SE object
  #' @export
  #'
  dt <- data.table::as.data.table(gsea_df)

  mats <- lapply(value_cols, function(val) {
    wide <- data.table::dcast(
      dt,
      as.formula(paste(row_col, "~", col_col)),
      value.var = val
    )
    mat <- as.matrix(wide[, -1, with = FALSE])
    rownames(mat) <- wide[[row_col]]
    mat
  })

  names(mats) <- value_cols
  mats
}

create_star_adder <- function(q_mat, threshold, fontsize = 15) {
  #' add start for sig pval on the heapmap
  #' @export
  #'
  if (!is.matrix(q_mat)) {
    stop("q_mat must be a matrix")
  }
  if (!any(q_mat < threshold)) {
    return(NULL)
  }

  function(j, i, x, y, w, h, fill) {
    v <- ComplexHeatmap::pindex(q_mat, i, j)
    l <- v < threshold
    grid::grid.text("*", x[l], y[l], gp = grid::gpar(fontsize = fontsize))
  }
}

make_heatmap_color <- function(break_low, break_high) {
  #' Create colors for heatmap
  #' @export
  #'
  col <- grDevices::colorRampPalette(c("blue", "white", "red"))(50)
  breaks <- seq(break_low, break_high, length.out = 50)
  circlize::colorRamp2(breaks, col)
}

draw_heatmap <- function(
  nes_mat,
  q_mat,
  col_df = NULL,
  group,
  cellwidth = 15,
  cellheight = 12
) {
  #' plot a heatmap using ComplexHeatmap
  #' @export
  #'
  layer_fun <- if (any(q_mat < 0.05)) {
    create_star_adder(q_mat, 0.05, 12)
  } else {
    NULL
  }

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(df = col_df)

  ComplexHeatmap::Heatmap(
    nes_mat,
    layer_fun = layer_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = make_heatmap_color(min(nes_mat), max(nes_mat)),
    width = ncol(nes_mat) * grid::unit(cellwidth, "pt"),
    height = nrow(nes_mat) * grid::unit(cellheight, "pt"),
    row_names_gp = grid::gpar(fontsize = 6, fontface = "italic"),
    top_annotation = top_annotation,
    name = "NES"
  )
}
