#' Save Analysis Parameters to File
#'
#' Writes the key analysis parameters to a plain-text file in `output_dir`.
#'
#' @param params Named list of parameter values.
#' @param output_dir Directory where the file is saved.
#' @return Invisibly the file path.
#' @export
save_analysis_parameters <- function(params, output_dir) {
  out_file <- file.path(output_dir, "analysis_parameters.txt")
  tryCatch({
    lines <- mapply(function(nm, val) {
      paste0(nm, " = ", paste(as.character(val), collapse = ", "))
    }, names(params), params)
    writeLines(lines, out_file)
    message("Parameters saved to: ", out_file)
  }, error = function(e)
    message("Could not save parameters: ", e$message))
  invisible(out_file)
}


#' Get Example DEG List
#'
#' Returns the path to the built-in example DEG list (human gene symbols).
#'
#' @return Character string path to the example file.
#' @export
example_deg_path <- function() {
  system.file("extdata", "example_deg_list.txt", package = "LLM_EnrichViz")
}


#' Load Example DEG List
#'
#' Reads and returns the built-in example list of human DEG symbols.
#'
#' @return Character vector of gene symbols.
#' @export
load_example_deg <- function() {
  path <- example_deg_path()
  if (!file.exists(path)) stop("Example file not found.")
  readLines(path, warn = FALSE)
}


#' Check Package Dependencies
#'
#' Verifies that all required Bioconductor and CRAN packages are installed and
#' prints a summary.
#'
#' @return Invisibly a named logical vector (TRUE = installed).
#' @export
check_dependencies <- function() {
  required <- c(
    "clusterProfiler", "ReactomePA", "GSVA", "limma", "msigdbr",
    "httr", "R6", "future", "furrr", "progressr",
    "dplyr", "tidyr", "ggplot2", "visNetwork",
    "stringr", "rmarkdown", "tibble"
  )
  status <- vapply(required, function(pkg)
    requireNamespace(pkg, quietly = TRUE), logical(1L))

  ok  <- names(status)[status]
  bad <- names(status)[!status]

  if (length(ok)  > 0) message("Installed   : ", paste(ok,  collapse = ", "))
  if (length(bad) > 0) message("NOT installed: ", paste(bad, collapse = ", "),
                                "\n  Install via: install.packages() or ",
                                "BiocManager::install()")
  invisible(status)
}
