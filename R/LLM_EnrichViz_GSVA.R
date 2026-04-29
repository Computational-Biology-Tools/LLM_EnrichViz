#' Run Gene Set Variation Analysis (GSVA)
#'
#' Computes per-sample pathway activity scores using \code{GSVA::gsva()}.
#' Gene sets are retrieved via \code{msigdbr} for the requested databases.
#'
#' @param expression_data Numeric matrix (genes x samples). Row names must be
#'   gene SYMBOL identifiers.
#' @param gene_mapping Data frame with columns \code{ENTREZID} and \code{SYMBOL}.
#' @param organism \code{"human"} or \code{"mouse"}.
#' @param databases Character vector of databases to query.
#' @param output_dir Directory for result files.
#'
#' @return Named list; each element is a list with \code{summary} (data frame)
#'   and \code{scores_matrix} (matrix of GSVA scores).
#'
#' @importFrom GSVA gsva
#' @importFrom dplyr arrange desc slice mutate
#' @export
run_GSVA <- function(expression_data, gene_mapping, organism,
                     databases, output_dir) {

  results   <- list()
  gene_sets <- get_gene_sets(databases, organism)
  expr_mat  <- as.matrix(expression_data)

  for (db_name in names(gene_sets)) {
    message("    GSVA/", db_name, " ...")
    results[[db_name]] <- tryCatch({

      gsva_scores <- GSVA::gsva(
        expr          = expr_mat,
        gset.idx.list = gene_sets[[db_name]],
        method        = "gsva",
        verbose       = FALSE,
        parallel.sz   = 1L
      )

      mean_scores <- rowMeans(gsva_scores, na.rm = TRUE)
      sd_scores   <- apply(gsva_scores, 1, stats::sd, na.rm = TRUE)

      df <- data.frame(
        GeneSet   = names(mean_scores),
        MeanScore = mean_scores,
        SDScore   = sd_scores,
        stringsAsFactors = FALSE
      ) %>%
        dplyr::arrange(dplyr::desc(abs(.data[["MeanScore"]]))) %>%
        dplyr::slice(1:min(nrow(.), 100))

      df$Genes <- sapply(df$GeneSet, function(gs) {
        gns <- gene_sets[[db_name]][[gs]]
        paste(gns[gns %in% gene_mapping$SYMBOL], collapse = ",")
      })

      write.table(df,
                  file.path(output_dir, paste0("GSVA_", db_name, "_summary.txt")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      write.table(gsva_scores,
                  file.path(output_dir, paste0("GSVA_", db_name, "_scores_matrix.txt")),
                  sep = "\t", quote = FALSE)

      list(summary = df, scores_matrix = gsva_scores)
    }, error = function(e) {
      message("    GSVA/", db_name, " error: ", e$message)
      NULL
    })
  }

  return(results)
}
