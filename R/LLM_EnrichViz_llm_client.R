#' Get Default Model Name for a LLM Provider
#'
#' @param provider One of \code{"gemini"}, \code{"openai"}, \code{"claude"},
#'   \code{"deepseek"}.
#' @return Character string with the default model name.
#' @export
get_default_model <- function(provider) {
  defaults <- list(
    gemini   = "gemini-1.5-flash-latest",
    openai   = "gpt-4o",
    claude   = "claude-sonnet-4-6",
    deepseek = "deepseek-chat"
  )
  model <- defaults[[provider]]
  if (is.null(model)) stop("Unknown provider: ", provider)
  return(model)
}

#' Get LLM Request Function for a Provider
#'
#' Returns the appropriate request function for the selected LLM provider.
#'
#' @param provider One of \code{"gemini"}, \code{"openai"}, \code{"claude"},
#'   \code{"deepseek"}.
#' @return A function with signature
#'   \code{function(prompt, temperature, max_output_tokens, api_key, llm_model,
#'   delay_seconds)}.
#' @export
get_llm_function <- function(provider) {
  switch(provider,
    "gemini"   = make_gemini_request,
    "openai"   = make_openai_request,
    "claude"   = make_claude_request,
    "deepseek" = make_deepseek_request,
    stop("Unknown provider: ", provider)
  )
}

# ── GEMINI ────────────────────────────────────────────────────────────────────

#' Send a Request to the Google Gemini API
#'
#' @param prompt Character string prompt.
#' @param temperature Numeric 0–1 (randomness).
#' @param max_output_tokens Integer maximum tokens in response.
#' @param api_key Character string API key.
#' @param llm_model Model name (e.g., \code{"gemini-1.5-flash-latest"}).
#' @param delay_seconds Pause (seconds) after the request.
#'
#' @return Character string response, or a named list with \code{error}.
#' @importFrom httr POST content_type_json content status_code
#' @export
make_gemini_request <- function(prompt, temperature, max_output_tokens,
                                 api_key, llm_model, delay_seconds) {
  model_query <- paste0(llm_model, ":generateContent")
  req_body <- list(
    contents         = list(parts = list(list(text = prompt))),
    generationConfig = list(temperature     = temperature,
                            maxOutputTokens = max_output_tokens)
  )
  response <- httr::POST(
    url   = paste0("https://generativelanguage.googleapis.com/v1beta/models/",
                   model_query),
    query = list(key = api_key),
    httr::content_type_json(),
    encode = "json",
    body   = req_body
  )
  Sys.sleep(delay_seconds)
  sc <- httr::status_code(response)
  if (sc == 200L) {
    cands <- httr::content(response)$candidates
    if (length(cands) > 0 && !is.null(cands[[1]]$content$parts[[1]]$text))
      return(cands[[1]]$content$parts[[1]]$text)
    return("No response from Gemini")
  }
  msg <- paste0("Gemini error ", sc, ": ",
                httr::content(response)$error$message)
  message(msg); return(list(error = msg))
}

# ── OPENAI ────────────────────────────────────────────────────────────────────

#' Send a Request to the OpenAI Chat Completions API
#'
#' @inheritParams make_gemini_request
#' @importFrom httr POST add_headers content_type_json content status_code
#' @export
make_openai_request <- function(prompt, temperature, max_output_tokens,
                                 api_key, llm_model, delay_seconds) {
  req_body <- list(
    model       = llm_model,
    messages    = list(list(role = "user", content = prompt)),
    temperature = temperature,
    max_tokens  = max_output_tokens
  )
  response <- httr::POST(
    url    = "https://api.openai.com/v1/chat/completions",
    httr::add_headers(Authorization = paste("Bearer", api_key)),
    httr::content_type_json(),
    encode = "json",
    body   = req_body
  )
  Sys.sleep(delay_seconds)
  sc <- httr::status_code(response)
  if (sc == 200L) {
    ct <- httr::content(response)
    return(ct$choices[[1]]$message$content)
  }
  msg <- paste0("OpenAI error ", sc, ": ",
                httr::content(response)$error$message)
  message(msg); return(list(error = msg))
}

# ── CLAUDE (Anthropic) ────────────────────────────────────────────────────────

#' Send a Request to the Anthropic Claude API
#'
#' @inheritParams make_gemini_request
#' @importFrom httr POST add_headers content_type_json content status_code
#' @export
make_claude_request <- function(prompt, temperature, max_output_tokens,
                                 api_key, llm_model, delay_seconds) {
  req_body <- list(
    model      = llm_model,
    max_tokens = max_output_tokens,
    messages   = list(list(role = "user", content = prompt))
  )
  response <- httr::POST(
    url = "https://api.anthropic.com/v1/messages",
    httr::add_headers(
      "x-api-key"         = api_key,
      "anthropic-version" = "2023-06-01"
    ),
    httr::content_type_json(),
    encode = "json",
    body   = req_body
  )
  Sys.sleep(delay_seconds)
  sc <- httr::status_code(response)
  if (sc == 200L) {
    ct <- httr::content(response)
    return(ct$content[[1]]$text)
  }
  msg <- paste0("Claude error ", sc, ": ",
                httr::content(response)$error$message)
  message(msg); return(list(error = msg))
}

# ── DEEPSEEK ──────────────────────────────────────────────────────────────────

#' Send a Request to the DeepSeek API
#'
#' @inheritParams make_gemini_request
#' @importFrom httr POST add_headers content_type_json content status_code
#' @export
make_deepseek_request <- function(prompt, temperature, max_output_tokens,
                                   api_key, llm_model, delay_seconds) {
  req_body <- list(
    model       = llm_model,
    messages    = list(list(role = "user", content = prompt)),
    temperature = temperature,
    max_tokens  = max_output_tokens
  )
  response <- httr::POST(
    url    = "https://api.deepseek.com/v1/chat/completions",
    httr::add_headers(Authorization = paste("Bearer", api_key)),
    httr::content_type_json(),
    encode = "json",
    body   = req_body
  )
  Sys.sleep(delay_seconds)
  sc <- httr::status_code(response)
  if (sc == 200L) {
    ct <- httr::content(response)
    return(ct$choices[[1]]$message$content)
  }
  msg <- paste0("DeepSeek error ", sc, ": ",
                httr::content(response)$error$message)
  message(msg); return(list(error = msg))
}
