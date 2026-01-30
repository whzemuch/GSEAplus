#' Extract the GSEA results table
#'
#' @param x An object with GSEA results.
#' @param ... Passed to methods.
#' @return A data.frame of GSEA results.
#' @export
gsea_result <- function(x, ...) {
  UseMethod("gsea_result")
}

#' @export
gsea_result.gsea_case <- function(x, ...) {
  x$gsea@result
}
