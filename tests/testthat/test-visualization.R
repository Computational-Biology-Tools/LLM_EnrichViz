test_that("generate_all_visualizations runs without error on empty inputs", {
  tmp <- tempdir()
  expect_silent(
    generate_all_visualizations(
      enrichment_results = list(),
      agent_results      = list(),
      final_response     = "",
      gene_symbols       = c("TP53", "MYC"),
      output_dir         = tmp
    )
  )
})

test_that("generate_all_visualizations handles ORA data frame", {
  tmp    <- tempdir()
  ora_df <- make_fake_ora_df()
  er     <- list(ORA = list(KEGG = ora_df))
  expect_silent(
    generate_all_visualizations(
      enrichment_results = er,
      agent_results      = list(),
      final_response     = "",
      gene_symbols       = c("TP53", "CDK2"),
      output_dir         = tmp
    )
  )
  pdf_file <- file.path(tmp, "ORA_KEGG_barplot.pdf")
  expect_true(file.exists(pdf_file))
  unlink(pdf_file)
})

test_that("network is built from TSV in final_response", {
  tmp <- tempdir()
  fake_response <- paste0(
    "Some integration text.\n",
    "```tsv\n",
    "Node1\tEdge\tNode2\tExplanation\n",
    "TP53\tassociated with\tCell Cycle\tGrounded in ORA and GSEA\n",
    "MYC\tactivates\tCCND1\tLeading edge gene (GSEA)\n",
    "```\n"
  )
  expect_silent(
    generate_all_visualizations(
      enrichment_results = list(),
      agent_results      = list(),
      final_response     = fake_response,
      gene_symbols       = c("TP53", "MYC", "CCND1"),
      output_dir         = tmp
    )
  )
  net_file <- file.path(tmp, "integrated_network.html")
  expect_true(file.exists(net_file))
  unlink(net_file)
})
