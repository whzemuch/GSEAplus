#' DESeq2 pseudobulk utilities
#'
#' Functions for running pseudobulk DESeq2 on SummarizedExperiment objects
#' and preparing contrast-level differential expression tables for downstream
#' GSEA and pathway analysis.

run_pseudobulk_deseq <- function(
  se,
  design_formula = ~age,
  contrasts = list(
    O_vs_Y = c("age", "O", "Y"),
    O_vs_M = c("age", "O", "M"),
    M_vs_Y = c("age", "M", "Y")
  ),
  min_count = 10,
  min_samples = 3
) {
  #' run  DESeq2 on pseudobulk data from single cell analysi
  #' @export
  dds <- DESeq2::DESeqDataSet(se, design = design_formula)
  keep <- rowSums(DESeq2::counts(dds) >= min_count) > min_samples
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)

  res_df <- data.table::rbindlist(
    lapply(names(contrasts), function(contrast_name) {
      res <- DESeq2::results(
        dds,
        contrast = contrasts[[contrast_name]],
        format = "DataFrame"
      )
      dt <- data.table::as.data.table(res, keep.rownames = "gene")
      dt[, contrast := contrast_name]
      dt
    }),
    use.names = TRUE,
    fill = TRUE
  )

  list(dds = dds, res_df = res_df)
}
