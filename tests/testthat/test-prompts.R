test_that("create_enrichviz_prompt returns character string", {
  fake_df <- data.frame(
    ID          = "hsa04110",
    Description = "Cell cycle",
    pvalue      = 0.001,
    Gene        = "TP53,CDK2",
    stringsAsFactors = FALSE
  )
  result <- create_enrichviz_prompt(
    analysis_type       = "ORA",
    database            = "KEGG",
    enrichment_results  = fake_df,
    gene_symbols        = c("TP53", "CDK2", "MYC"),
    experimental_design = "Test experiment"
  )
  expect_type(result, "character")
  expect_true(grepl("ORA", result))
  expect_true(grepl("KEGG", result))
  expect_true(grepl("TP53", result))
})

test_that("create_enrichviz_prompt handles NULL results gracefully", {
  result <- create_enrichviz_prompt(
    analysis_type      = "GSEA",
    database           = "GO",
    enrichment_results = NULL,
    gene_symbols       = c("BRCA1")
  )
  expect_type(result, "character")
  expect_true(grepl("No significant results", result))
})

test_that("create_enrichviz_prompt handles all four methods", {
  fake_df_gsea <- data.frame(
    ID = "GO:0001", Description = "immune response",
    NES = 1.8, pvalue = 0.01,
    enrichmentScore = 0.6, Gene = "IRF1,STAT1",
    stringsAsFactors = FALSE
  )
  for (method in c("ORA", "GSEA", "GSVA", "ROAST")) {
    df <- if (method == "GSEA") fake_df_gsea else
      data.frame(GeneSet = "TestSet", MeanScore = 0.5, SDScore = 0.1,
                 Genes = "GENE1", PValue = 0.01, Direction = "Up",
                 Description = "test", pvalue = 0.01, Gene = "GENE1",
                 ID = "ID1", stringsAsFactors = FALSE)
    res <- create_enrichviz_prompt(method, "KEGG", df, c("GENE1", "GENE2"))
    expect_type(res, "character")
  }
})

test_that("create_integration_prompt includes all method names", {
  agent_results <- list(
    list(agent_id = 1L, prompt_type = "ORA KEGG analysis",
         response = "Pathway A is enriched.")
  )
  result <- create_integration_prompt(
    agent_results       = agent_results,
    enrichment_results  = list(ORA = list(KEGG = NULL)),
    gene_symbols        = c("TP53", "MYC"),
    experimental_design = "Cancer study",
    analyses            = c("ORA", "GSEA")
  )
  expect_type(result, "character")
  expect_true(grepl("ORA", result))
  expect_true(grepl("GSEA", result))
  expect_true(grepl("TP53", result))
})
