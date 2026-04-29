# These tests require org.Hs.eg.db; they are skipped if it is not installed.

test_that("map_gene_ids_enrichviz errors on unsupported organism", {
  expect_error(
    map_gene_ids_enrichviz(deg_list = "TP53", gene_type = "SYMBOL",
                            organism = "zebrafish"),
    "Unsupported organism"
  )
})

test_that("map_gene_ids_enrichviz errors when no inputs provided", {
  skip_if_not_installed("org.Hs.eg.db")
  expect_error(
    map_gene_ids_enrichviz(gene_type = "SYMBOL", organism = "human"),
    "No gene identifiers"
  )
})

test_that("map_gene_ids_enrichviz maps SYMBOL → ENTREZID correctly", {
  skip_if_not_installed("org.Hs.eg.db")
  result <- map_gene_ids_enrichviz(
    deg_list  = c("TP53", "MYC", "BRCA1"),
    gene_type = "SYMBOL",
    organism  = "human"
  )
  expect_true(is.data.frame(result))
  expect_true("ENTREZID" %in% names(result))
  expect_true("SYMBOL"   %in% names(result))
  expect_true(nrow(result) > 0)
})
