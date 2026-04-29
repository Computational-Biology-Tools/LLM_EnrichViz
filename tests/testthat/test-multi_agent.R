test_that("EnrichVizAgent initialises correctly", {
  agent <- EnrichVizAgent$new(id = 1L, prompt = "Hello LLM",
                               prompt_type = "ORA KEGG analysis")
  expect_equal(agent$id, 1L)
  expect_equal(agent$prompt, "Hello LLM")
  expect_equal(agent$prompt_type, "ORA KEGG analysis")
  expect_null(agent$response)
})

test_that("EnrichVizEnvironment creates correct number of agents", {
  prompts      <- list(p1 = "Prompt one", p2 = "Prompt two")
  prompt_types <- list(p1 = "ORA KEGG", p2 = "GSEA GO")
  env <- EnrichVizEnvironment$new(prompts, prompt_types)
  expect_equal(length(env$agents), 2L)
  expect_equal(env$agents[[1]]$prompt, "Prompt one")
  expect_equal(env$agents[[2]]$prompt_type, "GSEA GO")
})

test_that("save_agent_responses writes a file", {
  tmp <- tempdir()
  agent_results <- list(
    list(agent_id = 1L, prompt_type = "ORA KEGG",
         response = "KEGG pathways: Cell cycle was enriched."),
    list(agent_id = 2L, prompt_type = "GSEA GO",
         response = list(error = "503 error"))
  )
  save_agent_responses(agent_results, tmp)
  out_file <- file.path(tmp, "agent_responses.txt")
  expect_true(file.exists(out_file))
  content <- readLines(out_file)
  expect_true(any(grepl("Agent 1", content)))
  expect_true(any(grepl("Cell cycle", content)))
  expect_true(any(grepl("Error", content)))
  unlink(out_file)
})
