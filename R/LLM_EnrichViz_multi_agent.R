#' EnrichViz Agent Class
#'
#' An R6 class representing a single LLM agent with a method-specific prompt.
#'
#' @importFrom R6 R6Class
#' @export
EnrichVizAgent <- R6::R6Class(
  "EnrichVizAgent",
  public = list(
    #' @field id Agent identifier (integer).
    id          = NULL,
    #' @field prompt The text prompt sent to the LLM.
    prompt      = NULL,
    #' @field prompt_type Descriptive label (e.g., "ORA KEGG analysis").
    prompt_type = NULL,
    #' @field response LLM response text (populated after execution).
    response    = NULL,

    #' @description Initialise a new agent.
    #' @param id Integer agent ID.
    #' @param prompt Character string prompt.
    #' @param prompt_type Character string label.
    initialize = function(id, prompt, prompt_type) {
      self$id          <- id
      self$prompt      <- prompt
      self$prompt_type <- prompt_type
    },

    #' @description Execute the agent — send the prompt to the LLM.
    #' @param llm_func LLM request function.
    #' @param temperature Numeric.
    #' @param max_output_tokens Integer.
    #' @param api_key Character string.
    #' @param llm_model Character string.
    #' @param delay_seconds Numeric.
    #' @return LLM response (character or error list).
    get_response = function(llm_func, temperature, max_output_tokens,
                             api_key, llm_model, delay_seconds) {
      self$response <- llm_func(
        prompt            = self$prompt,
        temperature       = temperature,
        max_output_tokens = max_output_tokens,
        api_key           = api_key,
        llm_model         = llm_model,
        delay_seconds     = delay_seconds
      )
      return(self$response)
    }
  )
)


#' EnrichViz Multi-Agent Environment
#'
#' An R6 class that manages a pool of \code{EnrichVizAgent} objects and runs
#' them in parallel using \code{future}/\code{furrr}.
#'
#' @importFrom R6 R6Class
#' @importFrom future plan multisession sequential
#' @importFrom furrr future_map
#' @importFrom progressr with_progress progressor
#' @export
EnrichVizEnvironment <- R6::R6Class(
  "EnrichVizEnvironment",
  public = list(
    #' @field agents List of \code{EnrichVizAgent} objects.
    agents = list(),

    #' @description Initialise the environment from a list of prompts.
    #' @param prompts Named list of prompt strings.
    #' @param prompt_types Named list of descriptive labels (same names).
    initialize = function(prompts, prompt_types) {
      for (i in seq_along(prompts)) {
        self$agents[[i]] <- EnrichVizAgent$new(
          id          = i,
          prompt      = prompts[[i]],
          prompt_type = prompt_types[[i]]
        )
      }
    },

    #' @description Run all agents in parallel with automatic retry.
    #' @param llm_func LLM request function.
    #' @param temperature Numeric.
    #' @param max_output_tokens Integer.
    #' @param api_key Character string.
    #' @param llm_model Character string.
    #' @param delay_seconds Numeric.
    #' @param num_workers Number of parallel workers.
    #' @param max_retries Maximum retry attempts on 503 errors. Default 3.
    #' @return List of result records with fields \code{agent_id},
    #'   \code{prompt_type}, and \code{response}.
    run_agents = function(llm_func, temperature, max_output_tokens,
                           api_key, llm_model, delay_seconds,
                           num_workers, max_retries = 3) {

      future::plan(future::multisession, workers = num_workers)
      on.exit(future::plan(future::sequential))

      progressr::with_progress({
        p <- progressr::progressor(along = self$agents)

        results <- furrr::future_map(self$agents, function(agent) {
          p(sprintf("Agent %d: %s", agent$id, agent$prompt_type))
          message("Agent ", agent$id, " starting: ", agent$prompt_type)

          retries  <- 0L
          response <- NULL

          while (retries <= max_retries) {
            response <- tryCatch(
              agent$get_response(llm_func, temperature, max_output_tokens,
                                 api_key, llm_model, delay_seconds),
              error = function(e) {
                message("Agent ", agent$id, " error: ", e$message)
                list(error = e$message)
              }
            )

            is_503 <- is.list(response) &&
              !is.null(response$error) &&
              grepl("503", response$error, fixed = TRUE)

            if (is_503) {
              retries <- retries + 1L
              message("Agent ", agent$id, " 503 — retry ", retries, "/", max_retries)
              Sys.sleep(delay_seconds * 2)
            } else {
              break
            }
          }

          if (retries > max_retries) {
            message("Agent ", agent$id, " failed after ", max_retries, " retries.")
            return(list(agent_id    = agent$id,
                        prompt_type = agent$prompt_type,
                        response    = list(error = "Max retries exceeded")))
          }

          message("Agent ", agent$id, " completed: ", agent$prompt_type)
          list(agent_id    = agent$id,
               prompt_type = agent$prompt_type,
               response    = response)
        }, .progress = FALSE)
      })

      return(results)
    }
  )
)


#' Save Agent Responses to Disk
#'
#' Writes all agent responses to a single text file in \code{output_dir}.
#'
#' @param agent_results List returned by \code{EnrichVizEnvironment$run_agents()}.
#' @param output_dir Directory where the file is saved.
#' @return Invisibly returns the file path.
#' @export
save_agent_responses <- function(agent_results, output_dir) {
  out_file <- file.path(output_dir, "agent_responses.txt")
  tryCatch({
    con <- file(out_file, "w")
    on.exit(close(con))
    for (r in agent_results) {
      cat(sprintf("=== Agent %d | %s ===\n", r$agent_id, r$prompt_type),
          file = con)
      if (is.null(r$response)) {
        cat("Task failed (NULL response).\n\n", file = con)
      } else if (is.list(r$response) && !is.null(r$response$error)) {
        cat("Error:", r$response$error, "\n\n", file = con)
      } else if (is.character(r$response)) {
        cat(r$response, "\n\n", file = con)
      } else {
        cat("Unknown response format.\n\n", file = con)
      }
    }
    message("Agent responses saved to: ", out_file)
  }, error = function(e) message("Could not save agent responses: ", e$message))
  invisible(out_file)
}
