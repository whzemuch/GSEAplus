#' Create a gsea_case object
#'
#' @param deg Differential expression results.
#' @param gsea A \code{gseaResult} object.
#' @param cell_type Cell type label.
#' @param pathway Pathway (gene set) identifier.
#' @param contrast Contrast label.
#' @param cohort Cohort label.
#' @return An object of class \code{gsea_case}.
#' @export
new_gsea_case <- function(deg, gsea, cell_type, pathway,
                          contrast = "OvsY", cohort = "aging") {
  stopifnot(inherits(gsea, "gseaResult"))

  structure(
    list(
      deg = deg,
      gsea = gsea,
      meta = list(
        cell_type = cell_type,
        pathway = pathway,
        contrast = contrast,
        cohort = cohort
      )
    ),
    class = "gsea_case"
  )
}
