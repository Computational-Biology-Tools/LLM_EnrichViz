#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' Performs GSEA using a pre-ranked gene list. Internally calls
#' \code{clusterProfiler::gseKEGG}, \code{clusterProfiler::gseGO},
#' \code{ReactomePA::gsePathway}, and WikiPathways via
#' \code{msigdbr::msigdbr} + \code{clusterProfiler::GSEA}.
#'
#' @param ranked_list Named numeric vector (gene SYMBOL or ENTREZID as names,
#'   ranking statistic as values). Must be sorted in decreasing order.
#' @param gene_mapping Data frame with columns \code{ENTREZID} and \code{SYMBOL}.
#' @param gene_type Gene identifier type used in ranked_list names
#'   (one of \code{"SYMBOL"}, \code{"ENSEMBL"}, \code{"ENTREZID"}).
#' @param organism Organism name. Supported: human, mouse, rat, chicken,
#'   zebrafish, fly, worm, yeast, arabidopsis, pig, dog, cow, ecoli.
#' @param databases Character vector of databases to query.
#' @param pvalue Numeric p-value cutoff. Default \code{0.05}.
#' @param ont GO ontology (\code{"BP"}, \code{"CC"}, or \code{"MF"}).
#' @param output_dir Directory for result files.
#'
#' @return Named list of filtered data frames, one per database.
#'
#' @importFrom clusterProfiler gseKEGG gseGO bitr GSEA
#' @importFrom ReactomePA gsePathway
#' @importFrom dplyr filter arrange slice mutate left_join group_by summarize
#' @importFrom tidyr unnest
#' @importFrom msigdbr msigdbr
#' @export
run_GSEA <- function(ranked_list, gene_mapping, gene_type, organism,
                     databases, pvalue = 0.05, ont = "BP", output_dir) {

  if (!gene_type %in% c("SYMBOL", "ENSEMBL", "ENTREZID"))
    stop("'gene_type' must be one of: SYMBOL, ENSEMBL, ENTREZID.")

  org_info <- .get_org_info(organism)
  org_code <- org_info$kegg
  org_db   <- org_info$org_db
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop("Please install ", org_db, " via: BiocManager::install('", org_db, "')")
  }

  # Convert ranked list to ENTREZID-named vector
  input_names <- names(ranked_list)
  if (is.null(input_names)) stop("ranked_list must be a NAMED numeric vector.")

  ranked_entrez <- .ranked_to_entrez(ranked_list, gene_type, org_db)

  results <- list()

  if ("KEGG" %in% databases) {
    message("    GSEA/KEGG ...")
    results$KEGG <- tryCatch({
      res <- clusterProfiler::gseKEGG(
        geneList     = ranked_entrez,
        organism     = org_code,
        pvalueCutoff = pvalue,
        verbose      = FALSE,
        eps          = 0
      )
      .format_gsea(res, gene_mapping, output_dir, "GSEA_KEGG")
    }, error = function(e) { message("    GSEA/KEGG error: ", e$message); NULL })
  }

  if ("GO" %in% databases) {
    message("    GSEA/GO ...")
    results$GO <- tryCatch({
      res <- clusterProfiler::gseGO(
        geneList     = ranked_entrez,
        OrgDb        = org_db,
        ont          = ont,
        pvalueCutoff = pvalue,
        verbose      = FALSE,
        eps          = 0
      )
      .format_gsea(res, gene_mapping, output_dir, "GSEA_GO")
    }, error = function(e) { message("    GSEA/GO error: ", e$message); NULL })
  }

  if ("Reactome" %in% databases) {
    message("    GSEA/Reactome ...")
    if (!org_info$reactome) {
      message("    GSEA/Reactome skipped: ReactomePA only supports human/mouse.")
    } else {
      results$Reactome <- tryCatch({
        res <- ReactomePA::gsePathway(
          geneList     = ranked_entrez,
          organism     = organism,
          pvalueCutoff = pvalue,
          verbose      = FALSE,
          eps          = 0
        )
        .format_gsea(res, gene_mapping, output_dir, "GSEA_Reactome")
      }, error = function(e) { message("    GSEA/Reactome error: ", e$message); NULL })
    }
  }

  if ("WikiPathways" %in% databases) {
    message("    GSEA/WikiPathways ...")
    results$WikiPathways <- tryCatch({
      wp_species <- org_info$wp
      wp_df <- msigdbr::msigdbr(
        species     = wp_species,
        category    = "C2",
        subcategory = "CP:WIKIPATHWAYS"
      )
      term2gene <- wp_df[, c("gs_name", "entrez_gene")]
      term2gene <- term2gene[!is.na(term2gene$entrez_gene), , drop = FALSE]
      res <- clusterProfiler::GSEA(
        geneList     = ranked_entrez,
        TERM2GENE    = term2gene,
        pvalueCutoff = pvalue,
        verbose      = FALSE
      )
      .format_gsea(res, gene_mapping, output_dir, "GSEA_WP")
    }, error = function(e) { message("    GSEA/WP error: ", e$message); NULL })
  }

  return(results)
}

# Internal: convert ranked list to ENTREZID-named
.ranked_to_entrez <- function(ranked_list, gene_type, org_db) {
  names_in <- names(ranked_list)
  if (gene_type == "ENSEMBL") {
    names_in <- sub("\\..*$", "", names_in)
  }

  if (gene_type == "ENTREZID") {
    ranked_out <- ranked_list
    names(ranked_out) <- names_in
    ranked_out <- ranked_out[!is.na(names(ranked_out))]
    return(sort(ranked_out, decreasing = TRUE))
  }

  gene_df <- suppressWarnings(
    clusterProfiler::bitr(names_in,
                           fromType = gene_type,
                           toType   = "ENTREZID",
                           OrgDb    = org_db)
  )
  matched     <- match(gene_df[[gene_type]], names_in)
  ranked_out  <- ranked_list[matched]
  names(ranked_out) <- gene_df$ENTREZID
  ranked_out  <- ranked_out[!is.na(ranked_out)]
  ranked_out  <- sort(ranked_out, decreasing = TRUE)
  return(ranked_out)
}

# Internal helper: format and save GSEA results
.format_gsea <- function(res, gene_mapping, output_dir, prefix) {
  if (is.null(res) || nrow(res@result) == 0) return(NULL)

  df <- res@result %>%
    dplyr::arrange(.data[["pvalue"]]) %>%
    dplyr::slice(1:min(nrow(.), 100)) %>%
    dplyr::mutate(core_enrichment = strsplit(.data[["core_enrichment"]], "/")) %>%
    tidyr::unnest(.data[["core_enrichment"]]) %>%
    dplyr::left_join(gene_mapping,
                     by = c("core_enrichment" = "ENTREZID")) %>%
    dplyr::mutate(Gene = ifelse(!is.na(.data[["SYMBOL"]]),
                                .data[["SYMBOL"]],
                                .data[["core_enrichment"]])) %>%
    dplyr::group_by(.data[["ID"]], .data[["Description"]],
                    .data[["pvalue"]], .data[["NES"]],
                    .data[["enrichmentScore"]]) %>%
    dplyr::summarize(Gene = paste(unique(.data[["Gene"]]), collapse = ","),
                     .groups = "drop")

  write.table(df, file.path(output_dir, paste0(prefix, "_results.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  attr(df, "raw_result") <- res
  return(df)
}
