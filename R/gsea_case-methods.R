#' Print a gsea_case object
#'
#' @param x A \code{gsea_case}.
#' @param ... Unused.
#' @return The input object, invisibly.
#' @export
print.gsea_case <- function(x, ...) {
  cat(
    "<gsea_case>\n",
    " cohort   :", x$meta$cohort, "\n",
    " cell_type:", x$meta$cell_type, "\n",
    " contrast :", x$meta$contrast, "\n",
    " pathway  :", x$meta$pathway, "\n",
    sep = ""
  )
  invisible(x)
}

#' Plot a gsea_case object
#'
#' @param x A \code{gsea_case}.
#' @param pathway Optional pathway (gene set) identifier to override the stored
#'   one.
#' @param text_size Text size used for the NES/q-value annotation.
#' @param ... Passed to \code{gseaplot_with_stats()}.
#' @return A gseaplot object list returned by \code{gseaplot_with_stats()}.
#' @export
plot.gsea_case <- function(x, pathway = NULL, text_size = 6, ...) {
  geneSetID <- if (is.null(pathway)) x$meta$pathway else pathway

  if (!geneSetID %in% x$gsea@result$ID) {
    stop("Gene set not found: ", geneSetID)
  }

  title <- paste(
    x$meta$cohort,
    x$meta$cell_type,
    x$meta$contrast,
    geneSetID,
    sep = " | "
  )

  gseaplot_with_stats(
    x$gsea,
    geneSetID,
    title = title,
    text_size = text_size,
    ...
  )
}
