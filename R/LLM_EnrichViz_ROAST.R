#' Run ROAST Gene Set Testing
#'
#' Performs rotation-based gene set testing using \code{limma::mroast()}.
#' Identifies directional (Up/Down) and mixed differential enrichment.
#'
#' @param expression_data Numeric matrix (genes x samples). Row names should
#'   be gene SYMBOL identifiers.
#' @param design_matrix Numeric design matrix (samples x coefficients) as
#'   returned by \code{model.matrix()}.
#' @param contrast_vector Numeric contrast vector (e.g., \code{c(0, 1)} for
#'   the second coefficient). If \code{NULL}, the last coefficient is tested.
#' @param gene_mapping Data frame with columns \code{ENTREZID} and \code{SYMBOL}.
#' @param organism \code{"human"} or \code{"mouse"}.
#' @param databases Character vector of databases to query.
#' @param pvalue Numeric p-value cutoff. Default \code{0.05}.
#' @param nrot Number of rotations. Default \code{999}.
#' @param output_dir Directory for result files.
#'
#' @return Named list of filtered data frames, one per database.
#'
#' @importFrom limma mroast
#' @importFrom dplyr filter arrange mutate
#' @importFrom tibble rownames_to_column
#' @export
run_ROAST <- function(expression_data, design_matrix, contrast_vector = NULL,
                       gene_mapping, organism, databases,
                       pvalue = 0.05, nrot = 999, output_dir) {

  results   <- list()
  row_genes <- rownames(expression_data)
  gene_sets <- get_gene_sets(databases, organism)

  if (is.null(contrast_vector)) {
    contrast_vector <- rep(0, ncol(design_matrix))
    contrast_vector[ncol(design_matrix)] <- 1
    message("    contrast_vector not provided — testing last coefficient.")
  }

  for (db_name in names(gene_sets)) {
    message("    ROAST/", db_name, " ...")
    results[[db_name]] <- tryCatch({

      # Build index list (gene positions in expression matrix)
      index_list <- lapply(gene_sets[[db_name]], function(gs_genes) {
        which(row_genes %in% gs_genes)
      })
      # Keep gene sets with at least 5 mapped genes
      index_list <- index_list[sapply(index_list, length) >= 5]

      if (length(index_list) == 0) {
        message("    ROAST/", db_name, ": no gene sets with >= 5 genes.")
        return(NULL)
      }

      roast_res <- limma::mroast(
        y        = expression_data,
        index    = index_list,
        design   = design_matrix,
        contrast = contrast_vector,
        nrot     = nrot
      )

      df <- roast_res %>%
        tibble::rownames_to_column("GeneSet") %>%
        dplyr::filter(.data[["PValue"]] < pvalue) %>%
        dplyr::arrange(.data[["PValue"]]) %>%
        dplyr::slice(1:min(nrow(.), 100))

      df$Genes <- sapply(df$GeneSet, function(gs) {
        if (gs %in% names(gene_sets[[db_name]])) {
          gns <- gene_sets[[db_name]][[gs]]
          paste(gns[gns %in% gene_mapping$SYMBOL], collapse = ",")
        } else { "" }
      })

      write.table(df,
                  file.path(output_dir, paste0("ROAST_", db_name, "_results.txt")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      df
    }, error = function(e) {
      message("    ROAST/", db_name, " error: ", e$message)
      NULL
    })
  }

  return(results)
}


#' Retrieve Gene Sets for GSVA and ROAST
#'
#' Downloads gene sets from MSigDB (via \code{msigdbr}) or KEGG for use
#' in GSVA and ROAST analyses.
#'
#' @param databases Character vector of databases (\code{"KEGG"}, \code{"GO"},
#'   \code{"Reactome"}, \code{"WikiPathways"}).
#' @param organism \code{"human"} or \code{"mouse"}.
#'
#' @return Named list of gene set lists (names = gene set names, values =
#'   character vectors of gene SYMBOL).
#'
#' @importFrom msigdbr msigdbr
#' @importFrom limma getGeneKEGGLinks
#' @export
get_gene_sets <- function(databases, organism) {
  organism <- tolower(trimws(organism))
  if (!organism %in% c("human", "mouse"))
    stop("GSVA/ROAST gene sets are currently available only for human or mouse.")

  org_info     <- .get_org_info(organism)
  org_db       <- org_info$org_db
  species_name <- ifelse(organism == "human", "Homo sapiens", "Mus musculus")
  gene_sets    <- list()

  if ("KEGG" %in% databases) {
    tryCatch({
      org_code <- ifelse(organism == "human", "hsa", "mmu")
      kegg_df  <- limma::getGeneKEGGLinks(species.KEGG = org_code)
      if (!requireNamespace(org_db, quietly = TRUE)) {
        stop("Please install ", org_db, " via: BiocManager::install('", org_db, "')")
      }
      gene_map <- suppressWarnings(
        clusterProfiler::bitr(unique(kegg_df$GeneID),
                               fromType = "ENTREZID",
                               toType   = "SYMBOL",
                               OrgDb    = org_db)
      )
      kegg_df$GeneSymbol <- gene_map$SYMBOL[match(kegg_df$GeneID, gene_map$ENTREZID)]
      kegg_df <- kegg_df[!is.na(kegg_df$GeneSymbol), ]
      gene_sets$KEGG <- split(kegg_df$GeneSymbol, kegg_df$PathwayID)
    }, error = function(e) message("    get_gene_sets KEGG error: ", e$message))
  }

  if ("GO" %in% databases) {
    tryCatch({
      go_df <- msigdbr::msigdbr(species = species_name,
                                 category = "C5", subcategory = "GO:BP")
      gene_sets$GO <- split(go_df$gene_symbol, go_df$gs_name)
    }, error = function(e) message("    get_gene_sets GO error: ", e$message))
  }

  if ("Reactome" %in% databases) {
    tryCatch({
      r_df <- msigdbr::msigdbr(species = species_name,
                                category = "C2", subcategory = "CP:REACTOME")
      gene_sets$Reactome <- split(r_df$gene_symbol, r_df$gs_name)
    }, error = function(e) message("    get_gene_sets Reactome error: ", e$message))
  }

  if ("WikiPathways" %in% databases) {
    tryCatch({
      wp_df <- msigdbr::msigdbr(species = species_name,
                                 category = "C2", subcategory = "CP:WIKIPATHWAYS")
      gene_sets$WikiPathways <- split(wp_df$gene_symbol, wp_df$gs_name)
    }, error = function(e) message("    get_gene_sets WikiPathways error: ", e$message))
  }

  return(gene_sets)
}
