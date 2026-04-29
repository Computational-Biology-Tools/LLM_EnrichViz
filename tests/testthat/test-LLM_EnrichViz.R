library(testthat)
library(LLM_EnrichViz)

# ── Gene mapper ───────────────────────────────────────────────────────────────

test_that("map_gene_ids_enrichviz rejects bad organism", {
  expect_error(
    map_gene_ids_enrichviz(deg_list = "TP53", gene_type = "SYMBOL",
                           organism = "martian"),
    "organism"
  )
})

test_that("map_gene_ids_enrichviz rejects bad gene_type", {
  expect_error(
    map_gene_ids_enrichviz(deg_list = "TP53", gene_type = "HUGO",
                           organism = "human"),
    "gene_type"
  )
})

test_that("map_gene_ids_enrichviz errors when no genes provided", {
  expect_error(
    map_gene_ids_enrichviz(gene_type = "SYMBOL", organism = "human"),
    "No gene identifiers"
  )
})

# ── LLM client helpers ────────────────────────────────────────────────────────

test_that("get_default_model returns correct defaults", {
  expect_equal(get_default_model("gemini"),   "gemini-1.5-flash-latest")
  expect_equal(get_default_model("openai"),   "gpt-4o")
  expect_equal(get_default_model("claude"),   "claude-sonnet-4-6")
  expect_equal(get_default_model("deepseek"), "deepseek-chat")
})

test_that("get_default_model errors on unknown provider", {
  expect_error(get_default_model("bard"), "Unknown provider")
})

test_that("get_llm_function returns a function", {
  expect_true(is.function(get_llm_function("gemini")))
  expect_true(is.function(get_llm_function("openai")))
  expect_true(is.function(get_llm_function("claude")))
  expect_true(is.function(get_llm_function("deepseek")))
})

test_that("get_llm_function errors on unknown provider", {
  expect_error(get_llm_function("mistral"), "Unknown provider")
})

# ── Prompt builder ────────────────────────────────────────────────────────────

test_that("create_enrichviz_prompt returns string for NULL results", {
  out <- create_enrichviz_prompt("ORA", "KEGG", NULL, c("TP53", "MYC"))
  expect_true(is.character(out))
  expect_match(out, "No significant results")
})

test_that("create_enrichviz_prompt includes gene symbols", {
  fake_df <- data.frame(
    ID          = "hsa04110",
    Description = "Cell cycle",
    pvalue      = 0.001,
    Gene        = "TP53,CDK2",
    stringsAsFactors = FALSE
  )
  out <- create_enrichviz_prompt("ORA", "KEGG", fake_df, c("TP53", "CDK2", "MYC"))
  expect_match(out, "TP53")
  expect_match(out, "ORA")
  expect_match(out, "KEGG")
})

test_that("create_enrichviz_prompt includes experimental design when provided", {
  fake_df <- data.frame(
    ID = "GO:0006915", Description = "apoptotic process",
    pvalue = 0.01, Gene = "TP53",
    stringsAsFactors = FALSE
  )
  out <- create_enrichviz_prompt(
    "ORA", "GO", fake_df, "TP53",
    experimental_design = "Tumour vs normal, n=10"
  )
  expect_match(out, "Tumour vs normal")
})

test_that("GSEA prompt mentions NES", {
  fake_df <- data.frame(
    ID = "hsa04110", Description = "Cell cycle",
    pvalue = 0.01, NES = 1.8, enrichmentScore = 0.6, Gene = "TP53",
    stringsAsFactors = FALSE
  )
  out <- create_enrichviz_prompt("GSEA", "KEGG", fake_df, "TP53")
  expect_match(out, "NES")
})

test_that("ROAST prompt mentions direction", {
  fake_df <- data.frame(
    GeneSet = "KEGG_CELL_CYCLE", NGenes = 10,
    PropDown = 0.1, PropUp = 0.7,
    Direction = "Up", PValue = 0.005, FDR = 0.02,
    Genes = "TP53,CDK2",
    stringsAsFactors = FALSE
  )
  out <- create_enrichviz_prompt("ROAST", "KEGG", fake_df, c("TP53", "CDK2"))
  expect_match(out, "Direction|direction|Up|Down", perl = TRUE)
})

test_that("GSVA prompt mentions MeanScore", {
  fake_res <- list(
    summary = data.frame(
      GeneSet = "HALLMARK_MYC_TARGETS_V1",
      MeanScore = 0.45, SDScore = 0.12,
      Genes = "MYC,CDK4",
      stringsAsFactors = FALSE
    )
  )
  out <- create_enrichviz_prompt("GSVA", "GO", fake_res, c("MYC", "CDK4"))
  expect_match(out, "MeanScore|Mean")
})

# ── Integration prompt ────────────────────────────────────────────────────────

test_that("create_integration_prompt produces non-empty string", {
  fake_agents <- list(
    list(agent_id = 1L, prompt_type = "ORA KEGG analysis",
         response = "Pathway X is enriched.")
  )
  out <- create_integration_prompt(
    agent_results       = fake_agents,
    enrichment_results  = list(),
    gene_symbols        = c("TP53", "MYC"),
    experimental_design = NULL,
    analyses            = "ORA"
  )
  expect_true(is.character(out))
  expect_gt(nchar(out), 100L)
})

test_that("create_integration_prompt handles NULL agent response gracefully", {
  fake_agents <- list(
    list(agent_id = 1L, prompt_type = "ORA KEGG analysis", response = NULL)
  )
  out <- create_integration_prompt(
    agent_results      = fake_agents,
    enrichment_results = list(),
    gene_symbols       = "TP53",
    analyses           = "ORA"
  )
  expect_match(out, "failed")
})

# ── Orchestrator validation ───────────────────────────────────────────────────

test_that("LLM_EnrichViz stops without api_key_file", {
  expect_error(
    LLM_EnrichViz(deg_list = "TP53", organism = "human",
                  analyses = "ORA", api_key_file = NULL),
    "api_key_file"
  )
})

test_that("LLM_EnrichViz stops with unknown analysis", {
  tmp_key <- tempfile()
  writeLines("FAKE_KEY", tmp_key)
  expect_error(
    LLM_EnrichViz(deg_list = "TP53", analyses = "FAKE",
                  api_key_file = tmp_key),
    "Unknown analyses"
  )
  unlink(tmp_key)
})

test_that("LLM_EnrichViz stops when deg_list missing for ORA", {
  tmp_key <- tempfile()
  writeLines("FAKE_KEY", tmp_key)
  expect_error(
    LLM_EnrichViz(analyses = "ORA", api_key_file = tmp_key),
    "deg_list"
  )
  unlink(tmp_key)
})

test_that("LLM_EnrichViz stops when ranked_list missing for GSEA", {
  tmp_key <- tempfile()
  writeLines("FAKE_KEY", tmp_key)
  expect_error(
    LLM_EnrichViz(analyses = "GSEA", api_key_file = tmp_key),
    "ranked_list"
  )
  unlink(tmp_key)
})
