#' Run Over-Representation Analysis (ORA)
#'
#' Performs ORA across the requested pathway databases using a Fisher /
#' hypergeometric test. Internally calls \code{clusterProfiler::enrichKEGG},
#' \code{clusterProfiler::enrichGO}, \code{ReactomePA::enrichPathway}, and
#' \code{clusterProfiler::enrichWP}.
#'
#' @param gene_ids Character vector of ENTREZID identifiers (significant genes).
#' @param gene_mapping Data frame with columns \code{ENTREZID} and \code{SYMBOL}.
#' @param organism Organism name. Supported: human, mouse, rat, chicken,
#'   zebrafish, fly, worm, yeast, arabidopsis, pig, dog, cow, ecoli.
#' @param databases Character vector of databases to query.
#' @param pvalue Numeric p-value cutoff. Default \code{0.05}.
#' @param ont GO ontology (\code{"BP"}, \code{"CC"}, or \code{"MF"}).
#' @param output_dir Directory for result files.
#'
#' @return Named list of filtered data frames, one per database.
#'
#' @importFrom clusterProfiler enrichKEGG enrichGO enrichWP
#' @importFrom ReactomePA enrichPathway
#' @importFrom dplyr filter arrange slice mutate left_join group_by summarize
#' @importFrom tidyr unnest
#' @export
run_ORA <- function(gene_ids, gene_mapping, organism,
                    databases, pvalue = 0.05, ont = "BP", output_dir) {

  org_info <- .get_org_info(organism)
  org_code <- org_info$kegg
  org_db   <- org_info$org_db
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop("Please install ", org_db, " via: BiocManager::install('", org_db, "')")
  }

  results <- list()

  if ("KEGG" %in% databases) {
    message("    ORA/KEGG ...")
    results$KEGG <- tryCatch({
      res <- clusterProfiler::enrichKEGG(gene = gene_ids, organism = org_code,
                                          pvalueCutoff = 1)
      .format_ora(res, gene_mapping, pvalue, output_dir, "ORA_KEGG")
    }, error = function(e) { message("    ORA/KEGG error: ", e$message); NULL })
  }

  if ("GO" %in% databases) {
    message("    ORA/GO ...")
    results$GO <- tryCatch({
      res <- clusterProfiler::enrichGO(gene = gene_ids, OrgDb = org_db,
                                        ont = ont, pvalueCutoff = 1)
      .format_ora(res, gene_mapping, pvalue, output_dir, "ORA_GO")
    }, error = function(e) { message("    ORA/GO error: ", e$message); NULL })
  }

  if ("Reactome" %in% databases) {
    message("    ORA/Reactome ...")
    if (!org_info$reactome) {
      message("    ORA/Reactome skipped: ReactomePA only supports human/mouse.")
    } else {
      results$Reactome <- tryCatch({
        res <- ReactomePA::enrichPathway(gene = gene_ids, organism = organism,
                                          pvalueCutoff = 1)
        .format_ora(res, gene_mapping, pvalue, output_dir, "ORA_Reactome")
      }, error = function(e) { message("    ORA/Reactome error: ", e$message); NULL })
    }
  }

  if ("WikiPathways" %in% databases) {
    message("    ORA/WikiPathways ...")
    wp_org <- org_info$wp
    results$WikiPathways <- tryCatch({
      res <- clusterProfiler::enrichWP(gene = gene_ids, organism = wp_org)
      .format_ora(res, gene_mapping, pvalue, output_dir, "ORA_WP")
    }, error = function(e) { message("    ORA/WP error: ", e$message); NULL })
  }

  return(results)
}

# Internal helper: format and save ORA results
.format_ora <- function(res, gene_mapping, pvalue, output_dir, prefix) {
  if (is.null(res) || nrow(res@result) == 0) return(NULL)

  df <- res@result %>%
    dplyr::filter(.data[["pvalue"]] < pvalue) %>%
    dplyr::arrange(.data[["pvalue"]]) %>%
    dplyr::slice(1:min(nrow(.), 100)) %>%
    dplyr::mutate(geneID = strsplit(.data[["geneID"]], "/")) %>%
    tidyr::unnest(.data[["geneID"]]) %>%
    dplyr::left_join(gene_mapping, by = c("geneID" = "ENTREZID")) %>%
    dplyr::mutate(Gene = ifelse(!is.na(.data[["SYMBOL"]]),
                                .data[["SYMBOL"]], .data[["geneID"]])) %>%
    dplyr::group_by(.data[["ID"]], .data[["Description"]], .data[["pvalue"]],
                    .data[["qvalue"]], .data[["Count"]]) %>%
    dplyr::summarize(Gene = paste(unique(.data[["Gene"]]), collapse = ","),
                     .groups = "drop")

  write.table(df, file.path(output_dir, paste0(prefix, "_results.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  return(df)
}
