#' LLM_EnrichViz — Main Orchestrator
#'
#' Performs comprehensive gene set enrichment analyses (ORA, GSEA, GSVA, ROAST)
#' and interprets each analysis independently via a selected Large Language Model
#' (Gemini, OpenAI/GPT, Claude, or DeepSeek). A multi-agent system aggregates
#' per-method insights into an integrated mechanistic model with network
#' visualization and an automated HTML report.
#'
#' @param expression_data A numeric matrix of expression values (genes x samples),
#'   required for GSVA and ROAST. Row names must be gene identifiers matching
#'   \code{gene_type}.
#' @param ranked_list A named numeric vector for GSEA, where names are gene
#'   identifiers and values are ranking statistics (e.g., log2FC, signed
#'   -log10(p-value), or Wald statistic). Must be sorted in decreasing order.
#' @param deg_list A character vector of significant gene identifiers for ORA.
#' @param gene_type Character string specifying the gene identifier type used
#'   in inputs. One of \code{"SYMBOL"}, \code{"ENSEMBL"}, or \code{"ENTREZID"}.
#' @param organism Organism. Supported: human, mouse, rat, chicken, zebrafish,
#'   fly, worm, yeast, arabidopsis, pig, dog, cow, ecoli. Default is \code{"human"}.
#' @param design_matrix A numeric design matrix (samples x coefficients) for
#'   ROAST, typically produced by \code{model.matrix()}. Required if ROAST is
#'   in \code{analyses}.
#' @param contrast_vector A numeric contrast vector specifying the comparison
#'   for ROAST (e.g., \code{c(0, 1)} for the second coefficient). Required if
#'   ROAST is in \code{analyses}.
#' @param analyses Character vector specifying which analyses to run. Any
#'   combination of \code{"ORA"}, \code{"GSEA"}, \code{"GSVA"}, \code{"ROAST"}.
#'   Default is \code{c("ORA", "GSEA", "GSVA", "ROAST")}.
#' @param databases Character vector of pathway databases to query. Any
#'   combination of \code{"KEGG"}, \code{"GO"}, \code{"Reactome"},
#'   \code{"WikiPathways"}. Default is \code{c("KEGG", "GO", "Reactome",
#'   "WikiPathways")}.
#' @param pvalue Numeric p-value cutoff for filtering enrichment results.
#'   Default is \code{0.05}.
#' @param ont GO ontology for ORA/GSEA. One of \code{"BP"} (Biological Process),
#'   \code{"CC"} (Cellular Component), \code{"MF"} (Molecular Function).
#'   Default is \code{"BP"}.
#' @param llm_provider LLM provider to use. One of \code{"gemini"},
#'   \code{"openai"}, \code{"claude"}, \code{"deepseek"}.
#'   Default is \code{"gemini"}.
#' @param llm_model Specific model name (e.g., \code{"gpt-4o"},
#'   \code{"claude-sonnet-4-6"}, \code{"gemini-1.5-flash-latest"},
#'   \code{"deepseek-chat"}). If \code{NULL}, a sensible default is used for
#'   each provider.
#' @param api_key_file Path to a plain-text file containing the API key
#'   (one key per line; only the first line is read).
#' @param temperature Numeric value controlling LLM output randomness (0–1).
#'   Default is \code{0} (deterministic).
#' @param experimental_design A free-text description of the experimental
#'   design, included in all LLM prompts for biological context (optional).
#' @param output_dir Path to the directory where all results, plots, and the
#'   HTML report are saved. Created automatically if it does not exist.
#'   Default is \code{"enrichviz_results"}.
#'
#' @return Invisibly returns a named list with slots \code{enrichment_results},
#'   \code{agent_results}, and \code{final_response}. Primary output is written
#'   to \code{output_dir}.
#'
#' @details
#' ## Pipeline overview
#' \enumerate{
#'   \item **Gene ID mapping** — converts input identifiers to ENTREZID + SYMBOL.
#'   \item **Enrichment analyses** — runs the requested methods independently.
#'   \item **LLM prompt creation** — one prompt per (method × database) pair.
#'   \item **Multi-agent LLM execution** — agents run in parallel via
#'         \code{future}/\code{furrr}.
#'   \item **Integration prompt** — a synthesising LLM call across all agents.
#'   \item **Visualizations** — per-method plots + interactive network.
#'   \item **Report rendering** — an RMarkdown HTML report.
#' }
#'
#' ## Method-specific input requirements
#' | Method | Required inputs |
#' |--------|----------------|
#' | ORA    | \code{deg_list} |
#' | GSEA   | \code{ranked_list} |
#' | GSVA   | \code{expression_data} |
#' | ROAST  | \code{expression_data}, \code{design_matrix}, \code{contrast_vector} |
#'
#' @examples
#' \dontrun{
#' # Minimal ORA-only example
#' LLM_EnrichViz(
#'   deg_list            = c("TP53", "MYC", "BRCA1", "EGFR", "VEGFA"),
#'   gene_type           = "SYMBOL",
#'   organism            = "human",
#'   analyses            = "ORA",
#'   databases           = c("KEGG", "GO"),
#'   llm_provider        = "openai",
#'   llm_model           = "gpt-4o",
#'   api_key_file        = "~/.openai_key",
#'   experimental_design = "Breast cancer vs normal tissue, bulk RNA-Seq",
#'   output_dir          = "ora_results"
#' )
#'
#' # Full pipeline
#' LLM_EnrichViz(
#'   expression_data     = expr_matrix,
#'   ranked_list         = gene_rank_vec,
#'   deg_list            = sig_genes,
#'   gene_type           = "SYMBOL",
#'   organism            = "human",
#'   design_matrix       = model.matrix(~ condition, data = sample_info),
#'   contrast_vector     = c(0, 1),
#'   analyses            = c("ORA", "GSEA", "GSVA", "ROAST"),
#'   databases           = c("KEGG", "GO", "Reactome"),
#'   llm_provider        = "claude",
#'   llm_model           = "claude-sonnet-4-6",
#'   api_key_file        = "~/.claude_key",
#'   experimental_design = "PBMC: lupus patients vs healthy controls, scRNA-Seq",
#'   output_dir          = "full_results"
#' )
#' }
#'
#' @importFrom parallel detectCores
#' @importFrom rmarkdown render
#' @export
LLM_EnrichViz <- function(
    expression_data     = NULL,
    ranked_list         = NULL,
    deg_list            = NULL,
    gene_type           = NULL,
    organism            = "human",
    design_matrix       = NULL,
    contrast_vector     = NULL,
    analyses            = c("ORA", "GSEA", "GSVA", "ROAST"),
    databases           = c("KEGG", "GO", "Reactome", "WikiPathways"),
    pvalue              = 0.05,
    ont                 = "BP",
    llm_provider        = "gemini",
    llm_model           = NULL,
    api_key_file        = NULL,
    temperature         = 0,
    experimental_design = NULL,
    output_dir          = "enrichviz_results"
) {

  # ── 0. Validate inputs ────────────────────────────────────────────────────
  organism     <- tolower(organism)
  llm_provider <- tolower(llm_provider)

  org_info <- .get_org_info(organism)
  if (!llm_provider %in% c("gemini", "openai", "claude", "deepseek"))
    stop("'llm_provider' must be one of: gemini, openai, claude, deepseek.")
  if (is.null(api_key_file) || !file.exists(api_key_file))
    stop("'api_key_file' must point to an existing file containing your API key.")
  if (is.null(gene_type) || !gene_type %in% c("SYMBOL", "ENSEMBL", "ENTREZID"))
    stop("'gene_type' must be one of: SYMBOL, ENSEMBL, ENTREZID.")

  valid_analyses  <- c("ORA", "GSEA", "GSVA", "ROAST")
  invalid_methods <- setdiff(analyses, valid_analyses)
  if (length(invalid_methods) > 0)
    stop("Unknown analyses: ", paste(invalid_methods, collapse = ", "))

  if ("ORA" %in% analyses && is.null(deg_list))
    stop("'deg_list' must be provided when analyses includes 'ORA'.")
  if ("GSEA" %in% analyses && is.null(ranked_list))
    stop("'ranked_list' must be provided when analyses includes 'GSEA'.")
  if ("GSVA" %in% analyses && is.null(expression_data))
    stop("'expression_data' must be provided when analyses includes 'GSVA'.")
  if ("ROAST" %in% analyses &&
      (is.null(expression_data) || is.null(design_matrix)))
    stop("'expression_data' and 'design_matrix' must be provided when analyses includes 'ROAST'.")

  if (("GSVA" %in% analyses || "ROAST" %in% analyses) && gene_type != "SYMBOL")
    stop("GSVA/ROAST currently require gene_type = 'SYMBOL' for expression_data row names.")
  if (("GSVA" %in% analyses || "ROAST" %in% analyses) &&
      !organism %in% c("human", "mouse"))
    stop("GSVA/ROAST gene sets are currently available only for human or mouse.")

  # ── 1. Setup ──────────────────────────────────────────────────────────────
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  max_output_tokens <- 8192
  delay_seconds     <- 5
  num_workers       <- max(1L, parallel::detectCores() - 1L)

  api_key <- tryCatch(readLines(api_key_file, warn = FALSE)[1],
                      error = function(e) stop("Cannot read API key: ", e$message))
  if (is.na(api_key) || nchar(trimws(api_key)) == 0)
    stop("API key file is empty or unreadable.")

  if (is.null(llm_model)) llm_model <- get_default_model(llm_provider)
  llm_func <- get_llm_function(llm_provider)

  message("╔══════════════════════════════════════════════╗")
  message("║        LLM_EnrichViz Pipeline Started        ║")
  message("╚══════════════════════════════════════════════╝")
  message("Provider : ", llm_provider, " | Model: ", llm_model)
  message("Analyses : ", paste(analyses, collapse = ", "))
  message("Databases: ", paste(databases, collapse = ", "))
  message("Output   : ", output_dir)

  # ── 2. Gene ID mapping ───────────────────────────────────────────────────
  message("\n[Step 1/7] Gene ID mapping ...")
  gene_mapping <- map_gene_ids_enrichviz(
    deg_list        = deg_list,
    ranked_list     = ranked_list,
    expression_data = expression_data,
    gene_type       = gene_type,
    organism        = organism
  )
  gene_symbols <- unique(gene_mapping$SYMBOL)
  gene_ids     <- unique(gene_mapping$ENTREZID)
  message("  Mapped ", length(gene_symbols), " genes.")

  # ── 3. Enrichment analyses ────────────────────────────────────────────────
  message("\n[Step 2/7] Running enrichment analyses ...")
  enrichment_results <- list()

  if ("ORA" %in% analyses) {
    message("  → ORA ...")
    enrichment_results$ORA <- run_ORA(
      gene_ids     = gene_ids,
      gene_mapping = gene_mapping,
      organism     = organism,
      databases    = databases,
      pvalue       = pvalue,
      ont          = ont,
      output_dir   = output_dir
    )
  }

  if ("GSEA" %in% analyses) {
    message("  → GSEA ...")
    enrichment_results$GSEA <- run_GSEA(
      ranked_list  = ranked_list,
      gene_mapping = gene_mapping,
      gene_type    = gene_type,
      organism     = organism,
      databases    = databases,
      pvalue       = pvalue,
      ont          = ont,
      output_dir   = output_dir
    )
  }

  if ("GSVA" %in% analyses) {
    message("  → GSVA ...")
    enrichment_results$GSVA <- run_GSVA(
      expression_data = expression_data,
      gene_mapping    = gene_mapping,
      organism        = organism,
      databases       = databases,
      output_dir      = output_dir
    )
  }

  if ("ROAST" %in% analyses) {
    message("  → ROAST ...")
    enrichment_results$ROAST <- run_ROAST(
      expression_data = expression_data,
      design_matrix   = design_matrix,
      contrast_vector = contrast_vector,
      gene_mapping    = gene_mapping,
      organism        = organism,
      databases       = databases,
      pvalue          = pvalue,
      output_dir      = output_dir
    )
  }

  # ── 4. Build LLM prompts ──────────────────────────────────────────────────
  message("\n[Step 3/7] Building LLM prompts ...")
  prompts      <- list()
  prompt_types <- list()

  for (analysis in names(enrichment_results)) {
    for (db in names(enrichment_results[[analysis]])) {
      key <- paste0(analysis, "_", db)
      prompts[[key]] <- create_enrichviz_prompt(
        analysis_type       = analysis,
        database            = db,
        enrichment_results  = enrichment_results[[analysis]][[db]],
        gene_symbols        = gene_symbols,
        experimental_design = experimental_design
      )
      prompt_types[[key]] <- paste(analysis, db, "analysis")
    }
  }
  message("  Created ", length(prompts), " prompt(s).")

  # ── 5. Multi-agent LLM execution ─────────────────────────────────────────
  message("\n[Step 4/7] Running LLM agents in parallel ...")
  env <- EnrichVizEnvironment$new(prompts, prompt_types)
  agent_results <- env$run_agents(
    llm_func          = llm_func,
    temperature       = temperature,
    max_output_tokens = max_output_tokens,
    api_key           = api_key,
    llm_model         = llm_model,
    delay_seconds     = delay_seconds,
    num_workers       = num_workers
  )
  save_agent_responses(agent_results, output_dir)

  # ── 6. Integration prompt ─────────────────────────────────────────────────
  message("\n[Step 5/7] Running integration LLM call ...")
  integration_prompt <- create_integration_prompt(
    agent_results       = agent_results,
    enrichment_results  = enrichment_results,
    gene_symbols        = gene_symbols,
    experimental_design = experimental_design,
    analyses            = names(enrichment_results)
  )
  final_response <- llm_func(
    prompt            = integration_prompt,
    temperature       = temperature,
    max_output_tokens = max_output_tokens,
    api_key           = api_key,
    llm_model         = llm_model,
    delay_seconds     = delay_seconds
  )
  writeLines(
    if (is.character(final_response)) final_response else "Error: no response.",
    file.path(output_dir, "final_response.txt")
  )

  # ── 7. Visualizations ────────────────────────────────────────────────────
  message("\n[Step 6/7] Generating visualizations ...")
  generate_all_visualizations(
    enrichment_results = enrichment_results,
    agent_results      = agent_results,
    final_response     = final_response,
    gene_symbols       = gene_symbols,
    output_dir         = output_dir
  )

  # ── 8. Report ─────────────────────────────────────────────────────────────
  message("\n[Step 7/7] Rendering HTML report ...")
  rmd_path <- system.file("rmd", "report_template.Rmd", package = "LLM_EnrichViz")
  if (file.exists(rmd_path)) {
    tryCatch(
      rmarkdown::render(rmd_path, output_dir = output_dir,
                        params = list(output_dir = output_dir)),
      error = function(e) message("Report rendering failed: ", e$message)
    )
  } else {
    message("  Report template not found — skipping.")
  }

  message("\n╔══════════════════════════════════════════════╗")
  message("║  Pipeline complete. Results in: ", output_dir)
  message("╚══════════════════════════════════════════════╝")

  invisible(list(
    enrichment_results = enrichment_results,
    agent_results      = agent_results,
    final_response     = final_response
  ))
}
