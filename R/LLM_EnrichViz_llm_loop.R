#' Run LLM Analysis for a Single Enrichment Method
#'
#' Executes the LLM prompt + multi-agent loop for a single analysis type
#' (ORA, GSEA, GSVA, or ROAST) using precomputed enrichment results.
#'
#' @param analysis_type One of \code{"ORA"}, \code{"GSEA"}, \code{"GSVA"},
#'   \code{"ROAST"}.
#' @param enrichment_results Named list of enrichment results by database
#'   (e.g., output of \code{run_ORA()} for one method).
#' @param gene_symbols Character vector of gene symbols used in the analysis.
#' @param experimental_design Optional free-text experimental design string.
#' @param llm_provider LLM provider (\code{"gemini"}, \code{"openai"},
#'   \code{"claude"}, \code{"deepseek"}).
#' @param llm_model Specific model name. If \code{NULL}, a default is used.
#' @param api_key_file Path to a plain-text file containing the API key.
#' @param temperature Numeric value controlling randomness (0-1).
#' @param max_output_tokens Integer maximum tokens for the LLM response.
#' @param delay_seconds Pause (seconds) after each request.
#' @param num_workers Number of parallel LLM agents.
#' @param output_dir Directory where agent responses are saved.
#'
#' @return Named list with \code{prompts} and \code{agent_results}.
#' @export
run_llm_for_analysis <- function(
    analysis_type,
    enrichment_results,
    gene_symbols,
    experimental_design = NULL,
    llm_provider = "gemini",
    llm_model = NULL,
    api_key_file = NULL,
    temperature = 0,
    max_output_tokens = 8192,
    delay_seconds = 5,
    num_workers = max(1L, parallel::detectCores() - 1L),
    output_dir = "llm_results"
) {
  analysis_type <- toupper(trimws(analysis_type))
  if (!analysis_type %in% c("ORA", "GSEA", "GSVA", "ROAST")) {
    stop("'analysis_type' must be one of: ORA, GSEA, GSVA, ROAST.")
  }
  if (!is.list(enrichment_results)) {
    stop("'enrichment_results' must be a named list of results by database.")
  }
  if (is.null(names(enrichment_results)) || any(names(enrichment_results) == "")) {
    stop("'enrichment_results' must have database names (e.g., KEGG, GO).")
  }

  llm_provider <- tolower(trimws(llm_provider))
  if (!llm_provider %in% c("gemini", "openai", "claude", "deepseek")) {
    stop("'llm_provider' must be one of: gemini, openai, claude, deepseek.")
  }
  if (is.null(api_key_file) || !file.exists(api_key_file)) {
    stop("'api_key_file' must point to an existing file containing your API key.")
  }

  api_key <- tryCatch(readLines(api_key_file, warn = FALSE)[1],
                      error = function(e) stop("Cannot read API key: ", e$message))
  if (is.na(api_key) || nchar(trimws(api_key)) == 0) {
    stop("API key file is empty or unreadable.")
  }

  if (is.null(llm_model)) llm_model <- get_default_model(llm_provider)
  llm_func <- get_llm_function(llm_provider)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  prompts      <- list()
  prompt_types <- list()
  for (db in names(enrichment_results)) {
    prompts[[db]] <- create_enrichviz_prompt(
      analysis_type       = analysis_type,
      database            = db,
      enrichment_results  = enrichment_results[[db]],
      gene_symbols        = gene_symbols,
      experimental_design = experimental_design
    )
    prompt_types[[db]] <- paste(analysis_type, db, "analysis")
  }

  env <- EnrichVizEnvironment$new(prompts, prompt_types)
  agent_results <- env$run_agents(
    llm_func          = llm_func,
    temperature       = temperature,
    max_output_tokens = max_output_tokens,
    api_key           = api_key,
    llm_model         = llm_model,
    delay_seconds     = delay_seconds,
    num_workers       = max(1L, num_workers)
  )

  save_agent_responses(agent_results, output_dir)

  list(prompts = prompts, agent_results = agent_results)
}

#' Run LLM Analysis for ORA
#'
#' Convenience wrapper for \code{run_llm_for_analysis(analysis_type = "ORA")}.
#'
#' @inheritParams run_llm_for_analysis
#' @export
LLM_ORA <- function(
    enrichment_results,
    gene_symbols,
    experimental_design = NULL,
    llm_provider = "gemini",
    llm_model = NULL,
    api_key_file = NULL,
    temperature = 0,
    max_output_tokens = 8192,
    delay_seconds = 5,
    num_workers = max(1L, parallel::detectCores() - 1L),
    output_dir = "llm_results"
) {
  run_llm_for_analysis(
    analysis_type       = "ORA",
    enrichment_results  = enrichment_results,
    gene_symbols        = gene_symbols,
    experimental_design = experimental_design,
    llm_provider        = llm_provider,
    llm_model           = llm_model,
    api_key_file        = api_key_file,
    temperature         = temperature,
    max_output_tokens   = max_output_tokens,
    delay_seconds       = delay_seconds,
    num_workers         = num_workers,
    output_dir          = output_dir
  )
}

#' Run LLM Analysis for GSEA
#'
#' Convenience wrapper for \code{run_llm_for_analysis(analysis_type = "GSEA")}.
#'
#' @inheritParams run_llm_for_analysis
#' @export
LLM_GSEA <- function(
    enrichment_results,
    gene_symbols,
    experimental_design = NULL,
    llm_provider = "gemini",
    llm_model = NULL,
    api_key_file = NULL,
    temperature = 0,
    max_output_tokens = 8192,
    delay_seconds = 5,
    num_workers = max(1L, parallel::detectCores() - 1L),
    output_dir = "llm_results"
) {
  run_llm_for_analysis(
    analysis_type       = "GSEA",
    enrichment_results  = enrichment_results,
    gene_symbols        = gene_symbols,
    experimental_design = experimental_design,
    llm_provider        = llm_provider,
    llm_model           = llm_model,
    api_key_file        = api_key_file,
    temperature         = temperature,
    max_output_tokens   = max_output_tokens,
    delay_seconds       = delay_seconds,
    num_workers         = num_workers,
    output_dir          = output_dir
  )
}

#' Run LLM Analysis for GSVA
#'
#' Convenience wrapper for \code{run_llm_for_analysis(analysis_type = "GSVA")}.
#'
#' @inheritParams run_llm_for_analysis
#' @export
LLM_GSVA <- function(
    enrichment_results,
    gene_symbols,
    experimental_design = NULL,
    llm_provider = "gemini",
    llm_model = NULL,
    api_key_file = NULL,
    temperature = 0,
    max_output_tokens = 8192,
    delay_seconds = 5,
    num_workers = max(1L, parallel::detectCores() - 1L),
    output_dir = "llm_results"
) {
  run_llm_for_analysis(
    analysis_type       = "GSVA",
    enrichment_results  = enrichment_results,
    gene_symbols        = gene_symbols,
    experimental_design = experimental_design,
    llm_provider        = llm_provider,
    llm_model           = llm_model,
    api_key_file        = api_key_file,
    temperature         = temperature,
    max_output_tokens   = max_output_tokens,
    delay_seconds       = delay_seconds,
    num_workers         = num_workers,
    output_dir          = output_dir
  )
}

#' Run LLM Analysis for ROAST
#'
#' Convenience wrapper for \code{run_llm_for_analysis(analysis_type = "ROAST")}.
#'
#' @inheritParams run_llm_for_analysis
#' @export
LLM_ROAST <- function(
    enrichment_results,
    gene_symbols,
    experimental_design = NULL,
    llm_provider = "gemini",
    llm_model = NULL,
    api_key_file = NULL,
    temperature = 0,
    max_output_tokens = 8192,
    delay_seconds = 5,
    num_workers = max(1L, parallel::detectCores() - 1L),
    output_dir = "llm_results"
) {
  run_llm_for_analysis(
    analysis_type       = "ROAST",
    enrichment_results  = enrichment_results,
    gene_symbols        = gene_symbols,
    experimental_design = experimental_design,
    llm_provider        = llm_provider,
    llm_model           = llm_model,
    api_key_file        = api_key_file,
    temperature         = temperature,
    max_output_tokens   = max_output_tokens,
    delay_seconds       = delay_seconds,
    num_workers         = num_workers,
    output_dir          = output_dir
  )
}
