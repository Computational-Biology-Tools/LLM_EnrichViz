#' Create a Multi-Method Integration Prompt
#'
#' Builds a final synthesis prompt that aggregates all agent responses and
#' instructs the LLM to identify cross-method convergence, complementary
#' insights, high-confidence hub genes, and a mechanistic model.
#'
#' @param agent_results List of agent result records (from
#'   \code{EnrichVizEnvironment$run_agents()}).
#' @param enrichment_results Named list of enrichment result lists, keyed by
#'   method (\code{"ORA"}, \code{"GSEA"}, \code{"GSVA"}, \code{"ROAST"}).
#' @param gene_symbols Character vector of gene symbols.
#' @param experimental_design Optional experimental design string.
#' @param analyses Character vector of methods that were run.
#'
#' @return A character string containing the integration prompt.
#' @export
create_integration_prompt <- function(agent_results, enrichment_results,
                                       gene_symbols, experimental_design = NULL,
                                       analyses) {

  # Header
  prompt <- paste0(
    "# Multi-Method Enrichment Integration\n\n",
    "**Analysed genes (", length(gene_symbols), "):**\n",
    paste(gene_symbols, collapse = ", "), "\n\n"
  )
  if (!is.null(experimental_design) && nchar(trimws(experimental_design)) > 0)
    prompt <- paste0(prompt,
                     "**Experimental design:** ", experimental_design, "\n\n")

  prompt <- paste0(prompt,
                   "**Methods run:** ", paste(analyses, collapse = ", "), "\n\n",
                   "---\n\n## Agent Responses\n\n")

  # Append each agent's response
  for (r in agent_results) {
    prompt <- paste0(prompt,
      "### Agent ", r$agent_id, " — ", r$prompt_type, "\n")
    if (is.null(r$response)) {
      prompt <- paste0(prompt, "_Task failed (no response)._\n\n")
    } else if (is.list(r$response) && !is.null(r$response$error)) {
      prompt <- paste0(prompt, "_Error:_ ", r$response$error, "\n\n")
    } else if (is.character(r$response)) {
      prompt <- paste0(prompt, r$response, "\n\n")
    }
  }

  # Synthesis instructions
  prompt <- paste0(prompt,
    "---\n\n## Integration Instructions\n\n",

    "### Step 1 — Cross-Method Convergence\n",
    "Identify pathways and genes detected by **multiple** methods. ",
    "Rank them by the number of methods that detected them (confidence tier: ",
    "4 methods = highest, 1 method = lowest). Explain *why* each method ",
    "captures or misses a given signal.\n\n",

    "### Step 2 — Complementary Insights per Method\n",
    "For each method, describe the unique biological information it contributes:\n",
    "- **ORA**: static over-representation (which gene sets are hit?)\n",
    "- **GSEA**: ranked enrichment (direction and magnitude of gene set shifts)\n",
    "- **GSVA**: per-sample variation (heterogeneity, sub-populations)\n",
    "- **ROAST**: coordinated directional expression (pathway coherence)\n\n",

    "### Step 3 — High-Confidence Hub Genes\n",
    "Identify genes that appear in:\n",
    "- ORA enriched sets\n",
    "- GSEA leading edges\n",
    "- ROAST directional sets (Up or Down)\n",
    "- Highest-scoring GSVA gene sets\n",
    "List hub genes with a confidence score (number of methods supporting them). ",
    "For each hub gene, note known drug targets, kinases, or ligands.\n\n",

    "### Step 4 — Integrated Mechanistic Model\n",
    "Propose a unified mechanistic model that incorporates findings from all ",
    "methods. Include:\n",
    "- Key regulatory nodes (upstream TFs, kinases)\n",
    "- Core effector pathways (highest-confidence)\n",
    "- Speculative pathways (low-confidence, single-method)\n",
    "- A short testable hypothesis\n\n",

    "### Step 5 — Title\n",
    "One title line (max 10 words) capturing the key finding, including one ",
    "hub gene name.\n\n",

    "### Step 6 — Integrated System Network (TSV)\n",
    "Generate a network in TSV format (4 columns: Node1, Edge, Node2, Explanation):\n",
    "- Node1: gene (uppercase) or pathway (Title Case)\n",
    "- Edge: 'activates', 'represses', 'associated with', 'targets'\n",
    "- Node2: gene (uppercase), pathway (Title Case), or GO term (Title Case)\n",
    "- Explanation: method(s) supporting the interaction, confidence tier\n",
    "Only include interactions supported by at least 2 methods.\n",
    "```tsv\nNode1\tEdge\tNode2\tExplanation\n```\n\n",

    "**Note:** Keep response under 6000 words.\n\n"
  )

  return(prompt)
}
