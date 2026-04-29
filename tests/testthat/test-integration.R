test_that("create_integration_prompt contains gene symbols", {
  syms <- c("TP53", "MYC", "BRCA1")
  ar   <- list(list(agent_id = 1L, prompt_type = "ORA KEGG",
                    response = "Cell cycle is enriched."))
  result <- create_integration_prompt(
    agent_results       = ar,
    enrichment_results  = list(ORA = list(KEGG = NULL)),
    gene_symbols        = syms,
    experimental_design = "Cancer study",
    analyses            = c("ORA")
  )
  for (g in syms) expect_true(grepl(g, result))
})

test_that("create_integration_prompt includes all six steps", {
  result <- create_integration_prompt(
    agent_results       = list(),
    enrichment_results  = list(),
    gene_symbols        = c("TP53"),
    analyses            = c("ORA", "GSEA", "GSVA", "ROAST")
  )
  for (step in c("Step 1", "Step 2", "Step 3", "Step 4", "Step 5", "Step 6")) {
    expect_true(grepl(step, result))
  }
})

test_that("create_integration_prompt handles agent error responses", {
  ar <- list(
    list(agent_id = 1L, prompt_type = "ORA KEGG",
         response = list(error = "503 Service Unavailable")),
    list(agent_id = 2L, prompt_type = "GSEA GO",
         response = "Immune pathways enriched (NES=1.9).")
  )
  result <- create_integration_prompt(
    agent_results       = ar,
    enrichment_results  = list(),
    gene_symbols        = c("IRF1"),
    analyses            = c("ORA", "GSEA")
  )
  expect_true(grepl("Error", result))
  expect_true(grepl("Immune", result))
})
