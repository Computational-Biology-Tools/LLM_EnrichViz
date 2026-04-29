#' Resolve organism name to database codes
#'
#' @param organism Character. Common name (e.g. "human", "mouse", "rat").
#' @return Named list with \code{org_db}, \code{kegg_code},
#'   \code{wp_name}, and \code{reactome_supported}.
#' @noRd
.get_org_info <- function(organism) {
  organism <- tolower(trimws(organism))
  mapping <- list(
    "human"      = list(org_db = "org.Hs.eg.db", kegg = "hsa",
                        wp = "Homo sapiens",      reactome = TRUE),
    "mouse"      = list(org_db = "org.Mm.eg.db", kegg = "mmu",
                        wp = "Mus musculus",      reactome = TRUE),
    "rat"        = list(org_db = "org.Rn.eg.db", kegg = "rno",
                        wp = "Rattus norvegicus", reactome = FALSE),
    "chicken"    = list(org_db = "org.Gg.eg.db", kegg = "gga",
                        wp = "Gallus gallus",     reactome = FALSE),
    "zebrafish"  = list(org_db = "org.Dr.eg.db", kegg = "dre",
                        wp = "Danio rerio",       reactome = FALSE),
    "fly"        = list(org_db = "org.Dm.eg.db", kegg = "dme",
                        wp = "Drosophila melanogaster", reactome = FALSE),
    "worm"       = list(org_db = "org.Ce.eg.db", kegg = "cel",
                        wp = "Caenorhabditis elegans", reactome = FALSE),
    "yeast"      = list(org_db = "org.Sc.sgd.db", kegg = "sce",
                        wp = "Saccharomyces cerevisiae", reactome = FALSE),
    "arabidopsis"= list(org_db = "org.At.tair.db", kegg = "ath",
                        wp = "Arabidopsis thaliana", reactome = FALSE),
    "pig"        = list(org_db = "org.Ss.eg.db", kegg = "ssc",
                        wp = "Sus scrofa",        reactome = FALSE),
    "dog"        = list(org_db = "org.Cf.eg.db", kegg = "cfa",
                        wp = "Canis familiaris",  reactome = FALSE),
    "cow"        = list(org_db = "org.Bt.eg.db", kegg = "bta",
                        wp = "Bos taurus",        reactome = FALSE),
    "ecoli"      = list(org_db = "org.EcK12.eg.db", kegg = "eco",
                        wp = "Escherichia coli",  reactome = FALSE)
  )
  if (!organism %in% names(mapping)) {
    stop("Unsupported organism: '", organism, "'. Supported: ",
         paste(names(mapping), collapse = ", "), ".")
  }
  info <- mapping[[organism]]
  info$name <- organism
  return(info)
}

#' Map Gene Identifiers for LLM_EnrichViz
#'
#' Converts gene identifiers from the input data (SYMBOL, ENSEMBL, or ENTREZID)
#' to a unified mapping table containing both ENTREZID and SYMBOL columns.
#' Handles the three input modes of LLM_EnrichViz: a discrete gene list
#' (ORA), a ranked vector (GSEA), and an expression matrix (GSVA/ROAST).
#'
#' @param deg_list Character vector of gene identifiers for ORA (optional).
#' @param ranked_list Named numeric vector for GSEA; names are gene identifiers
#'   (optional).
#' @param expression_data Numeric matrix for GSVA/ROAST; row names are gene
#'   identifiers (optional).
#' @param gene_type Source identifier type: \code{"SYMBOL"}, \code{"ENSEMBL"},
#'   or \code{"ENTREZID"}.
#' @param organism Organism name. Supported: human, mouse, rat, chicken,
#'   zebrafish, fly, worm, yeast, arabidopsis, pig, dog, cow, ecoli.
#'
#' @return A data frame with columns \code{ENTREZID} and \code{SYMBOL},
#'   deduplicated and with NA rows removed.
#'
#' @importFrom clusterProfiler bitr
#' @importFrom dplyr distinct filter
#' @export
map_gene_ids_enrichviz <- function(deg_list        = NULL,
                                    ranked_list     = NULL,
                                    expression_data = NULL,
                                    gene_type       = "SYMBOL",
                                    organism        = "human") {

  if (!gene_type %in% c("SYMBOL", "ENSEMBL", "ENTREZID")) {
    stop("'gene_type' must be one of: SYMBOL, ENSEMBL, ENTREZID.")
  }

  org_info <- .get_org_info(organism)
  org_db   <- org_info$org_db
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop("Please install ", org_db, " via: BiocManager::install('", org_db, "')")
  }

  # Collect all unique gene identifiers from all provided inputs
  all_ids <- unique(c(
    if (!is.null(deg_list))        as.character(deg_list),
    if (!is.null(ranked_list))     names(ranked_list),
    if (!is.null(expression_data)) rownames(expression_data)
  ))

  if (length(all_ids) == 0) stop("No gene identifiers found in any input.")

  # Strip ENSEMBL version suffix (e.g., ENSG00000139618.12 -> ENSG00000139618)
  if (gene_type == "ENSEMBL") {
    all_ids <- sub("\\..*$", "", all_ids)
  }

  if (gene_type == "ENTREZID") {
    # Already ENTREZID — map to SYMBOL only
    gene_mapping <- tryCatch(
      suppressWarnings(
        clusterProfiler::bitr(all_ids,
                               fromType = "ENTREZID",
                               toType   = "SYMBOL",
                               OrgDb    = org_db)
      ),
      error = function(e) {
        warning("bitr failed: ", e$message)
        data.frame(ENTREZID = character(0), SYMBOL = character(0))
      }
    )
  } else {
    gene_mapping <- tryCatch(
      suppressWarnings(
        clusterProfiler::bitr(all_ids,
                               fromType = gene_type,
                               toType   = c("ENTREZID", "SYMBOL"),
                               OrgDb    = org_db)
      ),
      error = function(e) {
        warning("bitr failed: ", e$message)
        data.frame(ENTREZID = character(0), SYMBOL = character(0))
      }
    )
  }

  gene_mapping <- gene_mapping[, c("ENTREZID", "SYMBOL"), drop = FALSE]
  gene_mapping <- gene_mapping[!is.na(gene_mapping$ENTREZID) &
                                 !is.na(gene_mapping$SYMBOL), ]
  gene_mapping <- unique(gene_mapping)

  if (nrow(gene_mapping) == 0) {
    warning("Gene ID mapping returned 0 rows. Check gene_type and organism.")
  } else {
    message("  Gene mapping: ", nrow(gene_mapping), " genes mapped.")
  }

  return(gene_mapping)
}
