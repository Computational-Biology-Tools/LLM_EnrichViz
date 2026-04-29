#' Create an LLM Prompt for a Given Enrichment Analysis
#'
#' Builds a method-specific prompt for the LLM agent, combining a standard
#' header (gene list, experimental design) with database-specific result
#' formatting and tailored analysis instructions.
#'
#' @param analysis_type One of \code{"ORA"}, \code{"GSEA"}, \code{"GSVA"},
#'   \code{"ROAST"}.
#' @param database Database name (e.g., \code{"KEGG"}, \code{"GO"}).
#' @param enrichment_results Data frame (or list with \code{summary}) of
#'   enrichment results for this method/database combination.
#' @param gene_symbols Character vector of gene symbols analysed.
#' @param experimental_design Optional free-text description of the experimental
#'   design.
#'
#' @return A character string containing the complete LLM prompt.
#' @export
create_enrichviz_prompt <- function(analysis_type, database,
                                     enrichment_results,
                                     gene_symbols,
                                     experimental_design = NULL) {

  if (is.null(enrichment_results) ||
      (is.data.frame(enrichment_results) && nrow(enrichment_results) == 0)) {
    return(paste0("No significant results available for ",
                  analysis_type, " / ", database, "."))
  }

  header  <- .build_base_header(gene_symbols, experimental_design)
  context <- switch(analysis_type,
    ORA   = .build_ORA_context(enrichment_results,  database),
    GSEA  = .build_GSEA_context(enrichment_results, database),
    GSVA  = .build_GSVA_context(enrichment_results, database),
    ROAST = .build_ROAST_context(enrichment_results, database),
    stop("Unknown analysis_type: ", analysis_type)
  )
  instructions <- switch(analysis_type,
    ORA   = .get_ORA_instructions(),
    GSEA  = .get_GSEA_instructions(),
    GSVA  = .get_GSVA_instructions(),
    ROAST = .get_ROAST_instructions()
  )

  paste0(header, context, instructions,
         "\n\n**Note:** Keep your response under 6000 words.\n\n")
}

# ── Header ────────────────────────────────────────────────────────────────────

.build_base_header <- function(gene_symbols, experimental_design) {
  h <- paste0(
    "I have performed a transcriptomics experiment and identified ",
    "differentially expressed genes (DEGs).\n\n",
    "**DEG list:**\n",
    paste(gene_symbols, collapse = ", "), "\n\n"
  )
  if (!is.null(experimental_design) && nchar(trimws(experimental_design)) > 0)
    h <- paste0(h, "**Experimental design:** ", experimental_design, "\n\n")
  h
}

# ── Context builders ──────────────────────────────────────────────────────────

.build_ORA_context <- function(results, database) {
  paste0(
    "**Method: ORA (Over-Representation Analysis) on ", database, "**\n",
    "ORA applies a hypergeometric/Fisher test to determine whether genes of ",
    "interest are statistically over-represented in predefined gene sets.\n",
    "Results are sorted by p-value (significant at p < 0.05).\n\n",
    "**ORA Results (", database, "):**\n",
    .fmt_table_generic(results, c("Description", "pvalue", "Gene")), "\n\n"
  )
}

.build_GSEA_context <- function(results, database) {
  paste0(
    "**Method: GSEA (Gene Set Enrichment Analysis) on ", database, "**\n",
    "GSEA evaluates whether members of a gene set are enriched at the top ",
    "or bottom of the ranked gene list. NES > 0 = upregulated/activated; ",
    "NES < 0 = downregulated/repressed.\n\n",
    "**GSEA Results (", database, "):**\n",
    .fmt_table_generic(results, c("Description", "NES", "pvalue", "Gene")), "\n\n"
  )
}

.build_GSVA_context <- function(results, database) {
  df <- if (is.list(results) && !is.null(results$summary)) results$summary else results
  paste0(
    "**Method: GSVA (Gene Set Variation Analysis) on ", database, "**\n",
    "GSVA computes a per-sample enrichment score for each gene set, capturing ",
    "inter-sample variation in pathway activity. MeanScore > 0 indicates ",
    "overall activation; MeanScore < 0 indicates repression.\n\n",
    "**GSVA Summary (", database, ") — Top 100 gene sets by |MeanScore|:**\n",
    .fmt_table_generic(df, c("GeneSet", "MeanScore", "SDScore", "Genes")), "\n\n"
  )
}

.build_ROAST_context <- function(results, database) {
  paste0(
    "**Method: ROAST (Rotation-based Gene Set Testing) on ", database, "**\n",
    "ROAST tests for coordinated differential expression within gene sets ",
    "using random rotations. Direction: Up = predominantly upregulated; ",
    "Down = predominantly downregulated; Mixed = heterogeneous.\n\n",
    "**ROAST Results (", database, "):**\n",
    .fmt_table_generic(results, c("GeneSet", "PValue", "Direction", "Genes")), "\n\n"
  )
}

# ── Table formatter ───────────────────────────────────────────────────────────

.fmt_table_generic <- function(df, cols) {
  if (is.null(df) || nrow(df) == 0) return("No significant results.\n")
  available <- intersect(cols, names(df))
  df_sub    <- df[, available, drop = FALSE]
  paste(apply(df_sub, 1, function(row) {
    paste(paste(available, row, sep = ": "), collapse = " | ")
  }), collapse = "\n")
}

# ── Instructions ──────────────────────────────────────────────────────────────

.get_ORA_instructions <- function() {
  paste0(
    "**Analysis Instructions (ORA):**\n\n",
    "1. **Summary & Categorization:** Provide a 4-line paragraph summarising ",
    "the enriched gene sets. Categorise them into the top 4 functional groups.\n\n",
    "2. **Known Roles & Upstream Regulators:** Summarise known roles in the ",
    "experimental context. Predict potential upstream regulators ",
    "(TFs, kinases) with reasoning.\n\n",
    "3. **Hub Gene Identification (step-by-step):**\n",
    "   - Step 1: Genes enriched in multiple pathways.\n",
    "   - Step 2: Highest-frequency genes across all enriched sets.\n",
    "   - Step 3: Integrate to name hub genes with justification.\n\n",
    "4. **ORA Limitations:** Note that ORA ignores expression magnitude and ",
    "non-significant genes; discuss potential bias.\n\n",
    "5. **Drug Targets / Kinases / Ligands:** For each hub gene, note if it ",
    "is a known drug target, kinase, or ligand.\n\n",
    "6. **Novelty Exploration:** Highlight unexpected findings.\n\n",
    "7. **Hypothesis:** Propose a mechanistic hypothesis.\n\n",
    "8. **System Network (TSV):**\n",
    "   Hub gene (uppercase) → pathway/GO term (Title Case).\n",
    "```tsv\nNode1\tEdge\tNode2\tExplanation\n```\n"
  )
}

.get_GSEA_instructions <- function() {
  paste0(
    "**Analysis Instructions (GSEA):**\n\n",
    "1. **NES Summary:** Summarise enriched gene sets, distinguishing ",
    "activated (NES > 0) from repressed (NES < 0).\n\n",
    "2. **Leading-Edge Genes:** Identify leading-edge (core enrichment) genes ",
    "shared across multiple enriched sets — these are the main drivers.\n\n",
    "3. **Directional Biology:** Interpret the biological significance of ",
    "enrichment direction given the experimental context.\n\n",
    "4. **GSEA Advantage:** Explain what GSEA detects that ORA would miss ",
    "(subtle signals below significance threshold).\n\n",
    "5. **Upstream Regulators & Hypothesis:** Propose regulators and a ",
    "mechanistic hypothesis.\n\n",
    "6. **Drug Targets / Kinases:** Note pharmacological relevance.\n\n",
    "7. **Novelty Exploration:** Unexpected connections.\n\n",
    "8. **System Network (TSV):**\n",
    "   Leading-edge gene (uppercase) → pathway (Title Case). ",
    "Include NES in Explanation.\n",
    "```tsv\nNode1\tEdge\tNode2\tExplanation\n```\n"
  )
}

.get_GSVA_instructions <- function() {
  paste0(
    "**Analysis Instructions (GSVA):**\n\n",
    "1. **Score Summary:** Summarise gene sets with the most extreme ",
    "MeanScore values (activated and repressed).\n\n",
    "2. **Inter-sample Variability:** Discuss what high SDScore values reveal ",
    "about sample heterogeneity or sub-populations.\n\n",
    "3. **Stable vs Variable Pathways:** Distinguish uniformly dysregulated ",
    "from highly variable pathways — what does variability imply?\n\n",
    "4. **Biomarker Potential:** Identify gene sets whose GSVA scores could ",
    "serve as pathway activity biomarkers.\n\n",
    "5. **Upstream Regulators & Hypothesis:** Propose regulators and a ",
    "mechanistic hypothesis explaining the observed variability.\n\n",
    "6. **Novelty Exploration:** Unexpected inter-sample patterns.\n\n",
    "7. **System Network (TSV):**\n",
    "   Gene set (Title Case) → key gene (uppercase). ",
    "Include MeanScore in Explanation.\n",
    "```tsv\nNode1\tEdge\tNode2\tExplanation\n```\n"
  )
}

.get_ROAST_instructions <- function() {
  paste0(
    "**Analysis Instructions (ROAST):**\n\n",
    "1. **Directional Summary:** Summarise significant gene sets, ",
    "distinguishing Up, Down, and Mixed.\n\n",
    "2. **Coordinated Expression:** Identify pathways with strongly ",
    "directional enrichment — what does coordinated up/down-regulation imply?\n\n",
    "3. **Statistical Robustness:** Discuss how the rotation-based approach ",
    "improves confidence relative to asymptotic tests.\n\n",
    "4. **Driver Genes:** For each significant gene set, identify the genes ",
    "most likely driving the directional signal.\n\n",
    "5. **Cross-Method Comparison:** If other methods were run, discuss ",
    "concordance and discordance with ROAST results.\n\n",
    "6. **Upstream Regulators & Hypothesis:** Propose regulators and a ",
    "mechanistic hypothesis.\n\n",
    "7. **Novelty Exploration:** Unexpected directional findings.\n\n",
    "8. **System Network (TSV):**\n",
    "   Driver gene (uppercase) → pathway (Title Case). ",
    "Include Direction in Explanation.\n",
    "```tsv\nNode1\tEdge\tNode2\tExplanation\n```\n"
  )
}
