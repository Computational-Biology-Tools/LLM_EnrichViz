# LLM_EnrichViz 0.1.0

## New features

* Initial release of `LLM_EnrichViz`.
* `LLM_EnrichViz()`: main orchestrator supporting ORA, GSEA, GSVA, and ROAST.
* `run_ORA()`: over-representation analysis via clusterProfiler / ReactomePA.
* `run_GSEA()`: gene set enrichment analysis via clusterProfiler / ReactomePA.
* `run_GSVA()`: gene set variation analysis via GSVA.
* `run_ROAST()`: rotation-based gene set testing via limma::mroast().
* `map_gene_ids_enrichviz()`: unified gene ID mapper (ENSEMBL / ENTREZID / SYMBOL).
* `get_gene_sets()`: gene set retrieval from KEGG, GO, Reactome, WikiPathways via msigdbr.
* LLM client supporting Gemini, OpenAI, Anthropic Claude, and DeepSeek.
* `EnrichVizEnvironment` / `EnrichVizAgent` R6 multi-agent parallel system.
* `create_enrichviz_prompt()`: method-aware prompt builder.
* `create_integration_prompt()`: cross-method synthesis prompt.
* `generate_all_visualizations()`: dot plots, NES plots, GSVA score bars, ROAST direction bars.
* Interactive HTML network via visNetwork.
* Automated HTML report via RMarkdown.
