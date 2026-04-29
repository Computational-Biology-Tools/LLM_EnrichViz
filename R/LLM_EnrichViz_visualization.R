#' Generate All Visualizations for LLM_EnrichViz
#'
#' Orchestrates the creation of per-method enrichment plots (bar, dot, GSEA
#' enrichment plots, GSVA heatmap) and the final interactive network.
#'
#' @param enrichment_results Named list of enrichment results by method.
#' @param agent_results List of agent result records.
#' @param final_response Character string: the final LLM integration response.
#' @param gene_symbols Character vector of gene symbols.
#' @param output_dir Directory for output files.
#'
#' @return Invisibly \code{NULL}.
#' @importFrom ggplot2 ggplot aes geom_col geom_point coord_flip theme_minimal
#'   scale_fill_gradient2 scale_color_gradient2 labs ggsave theme element_text
#' @importFrom visNetwork visNetwork visEdges visNodes visOptions visSave
#' @importFrom stringr str_match_all str_split
#' @importFrom dplyr arrange desc slice mutate
#' @export
generate_all_visualizations <- function(enrichment_results, agent_results,
                                         final_response, gene_symbols,
                                         output_dir) {

  # Per-method plots
  for (method in names(enrichment_results)) {
    for (db in names(enrichment_results[[method]])) {
      res <- enrichment_results[[method]][[db]]
      if (is.null(res)) next
      tryCatch({
        switch(method,
          ORA   = .plot_ora_bar(res, db, output_dir),
          GSEA  = .plot_gsea_bubble(res, db, output_dir),
          GSVA  = .plot_gsva_heatmap(res, db, output_dir),
          ROAST = .plot_roast_bar(res, db, output_dir)
        )
      }, error = function(e)
        message("Plot error [", method, "/", db, "]: ", e$message))
    }
  }

  # Final integrated network from LLM response
  if (is.character(final_response) && nchar(final_response) > 10) {
    tryCatch(
      .build_network_html(final_response, gene_symbols,
                          "integrated_network.html", output_dir),
      error = function(e)
        message("Network visualization error: ", e$message)
    )
  }

  invisible(NULL)
}

# ── ORA bar chart ─────────────────────────────────────────────────────────────

.plot_ora_bar <- function(df, db, output_dir) {
  if (!is.data.frame(df) || nrow(df) == 0) return(invisible(NULL))
  plot_df <- df %>%
    dplyr::arrange(.data[["pvalue"]]) %>%
    dplyr::slice(1:min(nrow(.), 20)) %>%
    dplyr::mutate(
      neg_log10p = -log10(.data[["pvalue"]]),
      Description = factor(.data[["Description"]],
                           levels = rev(.data[["Description"]]))
    )
  p <- ggplot2::ggplot(plot_df,
         ggplot2::aes(x = .data[["Description"]],
                      y = .data[["neg_log10p"]],
                      fill = .data[["neg_log10p"]])) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf",
                                   high = "#d73027", midpoint = 1.3) +
    ggplot2::labs(title = paste("ORA —", db, "(Top 20)"),
                  x = NULL, y = "-log10(p-value)", fill = "-log10(p)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  ggplot2::ggsave(file.path(output_dir, paste0("ORA_", db, "_barplot.pdf")),
                  plot = p, width = 10, height = 7)
  invisible(p)
}

# ── GSEA bubble plot ──────────────────────────────────────────────────────────

.plot_gsea_bubble <- function(df, db, output_dir) {
  if (!is.data.frame(df) || nrow(df) == 0 || !"NES" %in% names(df))
    return(invisible(NULL))
  plot_df <- df %>%
    dplyr::arrange(.data[["pvalue"]]) %>%
    dplyr::slice(1:min(nrow(.), 20)) %>%
    dplyr::mutate(
      NES = as.numeric(.data[["NES"]]),
      neg_log10p = -log10(as.numeric(.data[["pvalue"]])),
      Description = factor(.data[["Description"]],
                           levels = rev(.data[["Description"]]))
    )
  p <- ggplot2::ggplot(plot_df,
         ggplot2::aes(x = .data[["NES"]],
                      y = .data[["Description"]],
                      size = .data[["neg_log10p"]],
                      color = .data[["NES"]])) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_gradient2(low = "#4575b4", mid = "grey80",
                                    high = "#d73027", midpoint = 0) +
    ggplot2::labs(title = paste("GSEA —", db, "(Top 20)"),
                  x = "NES", y = NULL,
                  color = "NES", size = "-log10(p)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  ggplot2::ggsave(file.path(output_dir, paste0("GSEA_", db, "_bubble.pdf")),
                  plot = p, width = 10, height = 7)
  invisible(p)
}

# ── GSVA heatmap (simplified using ggplot2) ───────────────────────────────────

.plot_gsva_heatmap <- function(res, db, output_dir) {
  df <- if (is.list(res) && !is.null(res$summary)) res$summary else res
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  plot_df <- df %>%
    dplyr::arrange(dplyr::desc(abs(as.numeric(.data[["MeanScore"]])))) %>%
    dplyr::slice(1:min(nrow(.), 20)) %>%
    dplyr::mutate(
      MeanScore = as.numeric(.data[["MeanScore"]]),
      GeneSet   = factor(.data[["GeneSet"]],
                         levels = rev(.data[["GeneSet"]]))
    )
  p <- ggplot2::ggplot(plot_df,
         ggplot2::aes(x = 1, y = .data[["GeneSet"]],
                      fill = .data[["MeanScore"]])) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "white",
                                   high = "#d73027", midpoint = 0) +
    ggplot2::labs(title = paste("GSVA —", db, "(Top 20 by |MeanScore|)"),
                  x = NULL, y = NULL, fill = "Mean\nScore") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 8))
  ggplot2::ggsave(file.path(output_dir, paste0("GSVA_", db, "_heatmap.pdf")),
                  plot = p, width = 10, height = 7)
  invisible(p)
}

# ── ROAST bar chart ───────────────────────────────────────────────────────────

.plot_roast_bar <- function(df, db, output_dir) {
  if (!is.data.frame(df) || nrow(df) == 0) return(invisible(NULL))
  pval_col <- if ("PValue" %in% names(df)) "PValue" else "pvalue"
  dir_col  <- if ("Direction" %in% names(df)) "Direction" else NULL
  plot_df  <- df %>%
    dplyr::arrange(.data[[pval_col]]) %>%
    dplyr::slice(1:min(nrow(.), 20)) %>%
    dplyr::mutate(
      neg_log10p = -log10(as.numeric(.data[[pval_col]])),
      GeneSet    = factor(.data[["GeneSet"]],
                          levels = rev(.data[["GeneSet"]]))
    )
  fill_aes <- if (!is.null(dir_col)) {
    ggplot2::aes(x = .data[["GeneSet"]], y = .data[["neg_log10p"]],
                 fill = .data[[dir_col]])
  } else {
    ggplot2::aes(x = .data[["GeneSet"]], y = .data[["neg_log10p"]],
                 fill = .data[["neg_log10p"]])
  }
  p <- ggplot2::ggplot(plot_df, fill_aes) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste("ROAST —", db, "(Top 20)"),
                  x = NULL, y = "-log10(PValue)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  ggplot2::ggsave(file.path(output_dir, paste0("ROAST_", db, "_barplot.pdf")),
                  plot = p, width = 10, height = 7)
  invisible(p)
}

# ── Interactive network ───────────────────────────────────────────────────────

.build_network_html <- function(llm_response, gene_symbols,
                                 html_file, output_dir) {
  pattern <- "\\s*```(?i)tsv\\s*([\\s\\S]*?)\\s*```\\s*"
  matches <- stringr::str_match_all(llm_response, pattern)[[1]]
  if (nrow(matches) == 0) {
    message("No TSV block in final response — network skipped.")
    return(invisible(NULL))
  }
  tsv_content <- matches[nrow(matches), 2]  # Use last TSV block

  net_df <- tryCatch(
    read.table(text = tsv_content, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, fill = TRUE, quote = ""),
    error = function(e) { message("TSV parse error: ", e$message); NULL }
  )
  if (is.null(net_df) || ncol(net_df) < 3) return(invisible(NULL))
  colnames(net_df)[1:3] <- c("source", "interaction", "target")

  nodes <- data.frame(
    id    = unique(c(net_df$source, net_df$target)),
    label = unique(c(net_df$source, net_df$target)),
    stringsAsFactors = FALSE
  )
  nodes$group <- ifelse(nodes$id %in% toupper(gene_symbols), "Gene", "Process")
  nodes$color <- ifelse(nodes$group == "Gene", "#E8A838", "#6BAED6")

  edges <- data.frame(
    from  = net_df$source,
    to    = net_df$target,
    label = net_df$interaction,
    stringsAsFactors = FALSE
  )
  edges <- unique(edges)

  vis <- visNetwork::visNetwork(nodes = nodes, edges = edges,
                                 width = "100%", height = "700px") %>%
    visNetwork::visEdges(arrows = "to",
                          smooth = list(type = "cubicBezier")) %>%
    visNetwork::visNodes(size = 22,
                          font = list(size = 14, bold = TRUE)) %>%
    visNetwork::visOptions(highlightNearest   = TRUE,
                            selectedBy         = "group",
                            nodesIdSelection   = TRUE)

  visNetwork::visSave(vis, file = file.path(output_dir, html_file))
  message("Network saved: ", file.path(output_dir, html_file))
  invisible(vis)
}
