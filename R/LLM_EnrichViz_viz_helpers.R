#' Visualize ORA Results
#'
#' Generates a bar plot for ORA results (Top 20 by p-value).
#'
#' @param ora_results Data frame of ORA results or named list by database.
#' @param database Optional database name when a list is supplied.
#' @param output_dir Output directory for the plot.
#'
#' @return Invisibly returns the ggplot object.
#' @export
Viz_ORA <- function(ora_results, database = NULL, output_dir = "ora_viz") {
  if (is.list(ora_results) && !is.data.frame(ora_results)) {
    if (is.null(database)) {
      if (length(ora_results) != 1) {
        stop("Provide 'database' when ora_results has multiple entries.")
      }
      database <- names(ora_results)[1]
    }
    ora_results <- ora_results[[database]]
  }
  if (is.null(database) || database == "") database <- "ORA"
  if (is.null(ora_results)) return(invisible(NULL))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  p <- .plot_ora_bar(ora_results, database, output_dir)
  invisible(p)
}

#' Visualize GSEA Results
#'
#' Creates GSEA plots from a gseaResult object or the output of
#' \code{run_GSEA()}. Two styles are supported: \code{"simple"} uses
#' \code{enrichplot::gseaplot2()}, while \code{"diy"} generates a
#' multi-panel plot (running score + ranks + fold change) similar to the
#' reference tutorial.
#'
#' @param gsea_results A \code{gseaResult} object, a data frame returned by
#'   \code{run_GSEA()}, or a named list of GSEA results by database.
#' @param gene_set_id Character vector of pathway IDs to plot.
#' @param database Optional database name when \code{gsea_results} is a list.
#' @param style Either \code{"simple"} or \code{"diy"}.
#' @param ranked_list Optional named numeric vector used for DIY plotting.
#'   When \code{NULL}, the \code{geneList} from the GSEA object is used.
#' @param selected_genes Optional vector of gene symbols to label in DIY plots.
#' @param output_dir Output directory for plots and CSV files.
#' @param file_prefix Prefix for saved files.
#' @param save_running_scores Logical; save running-score CSV per pathway.
#' @param group_label Suffix used for running-score CSV filenames.
#' @param rel_heights Relative heights for DIY multi-panel layout.
#'
#' @return A named list of ggplot objects.
#' @export
Viz_GSEA <- function(
    gsea_results,
    gene_set_id,
    database = NULL,
    style = c("simple", "diy"),
    ranked_list = NULL,
    selected_genes = NULL,
    output_dir = "gsea_viz",
    file_prefix = "GSEA",
    save_running_scores = TRUE,
    group_label = "group1",
    rel_heights = c(1.5, 0.5, 1.5)
) {
  style <- match.arg(style)
  if (is.null(gene_set_id) || length(gene_set_id) == 0) {
    stop("'gene_set_id' must contain at least one pathway ID.")
  }
  gsea_raw <- .resolve_gsea_raw(gsea_results, database)

  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required for Viz_GSEA. Please install it.")
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  plots <- list()
  if (style == "simple") {
    for (gs_id in gene_set_id) {
      p <- enrichplot::gseaplot2(gsea_raw, geneSetID = gs_id, title = gs_id)
      out_file <- file.path(output_dir, paste0(file_prefix, "_", gs_id, ".pdf"))
      ggplot2::ggsave(out_file, plot = p, width = 7, height = 5)
      plots[[gs_id]] <- p
    }
    return(plots)
  }

  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for style = 'diy'.")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for style = 'diy'.")
  }

  gene_list <- ranked_list
  if (is.null(gene_list)) {
    gene_list <- gsea_raw@geneList
  }
  if (is.null(gene_list) || is.null(names(gene_list))) {
    stop("ranked_list (named numeric) is required for style = 'diy'.")
  }
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_symbols <- names(gene_list)

  gsdata <- do.call(rbind, lapply(gene_set_id, enrichplot:::gsInfo, object = gsea_raw))
  gsdata$gsym <- gene_symbols[gsdata$x]

  if (save_running_scores) {
    for (gs_id in gene_set_id) {
      gs_info <- enrichplot:::gsInfo(gsea_raw, gs_id)
      out_csv <- file.path(output_dir, paste0("gsea_genelist_", gs_id, "_",
                                              group_label, ".csv"))
      utils::write.csv(gs_info, out_csv, row.names = FALSE)
    }
  }

  # Running score plot
  p.res <- ggplot2::ggplot(gsdata, ggplot2::aes(x = x, y = runningScore,
                                               color = Description)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
    ggplot2::labs(x = NULL, y = "Enrichment\nScore") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.title = ggplot2::element_blank())

  # Rank plot
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }

  p.rank <- ggplot2::ggplot(gsdata, ggplot2::aes(x = x, ymin = ymin, ymax = ymax,
                                                color = Description)) +
    ggplot2::geom_linerange() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.position = "none",
                   axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank())

  # Fold-change plot
  df_fc <- data.frame(
    x    = seq_along(gene_list),
    y    = as.numeric(gene_list),
    gsym = gene_symbols,
    stringsAsFactors = FALSE
  )

  p.fc <- ggplot2::ggplot(df_fc, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_col(fill = "grey70", width = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
    ggplot2::labs(x = "Rank in ordered dataset", y = "Ranked list\nmetric") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  if (!is.null(selected_genes)) {
    sel_df <- df_fc[df_fc$gsym %in% selected_genes, , drop = FALSE]
    if (nrow(sel_df) > 0) {
      p.fc <- p.fc +
        ggrepel::geom_text_repel(
          data = sel_df,
          ggplot2::aes(label = gsym),
          direction = "x",
          angle = 90,
          size = 2.5,
          box.padding = grid::unit(0.35, "lines"),
          point.padding = grid::unit(0.3, "lines")
        )
    }
  }

  p.diy <- cowplot::plot_grid(
    p.res, p.rank, p.fc,
    ncol = 1,
    align = "v",
    rel_heights = rel_heights
  )

  out_file <- file.path(output_dir, paste0(file_prefix, "_multi_pathways.pdf"))
  ggplot2::ggsave(out_file, plot = p.diy, width = 7, height = 6)
  plots[["diy"]] <- p.diy

  plots
}

#' Visualize GSVA Results
#'
#' Generates GSVA visualizations. Use 
#' 
#' * `style = "heatmap"` for a top-mean-score heatmap, or
#' * `style = "limma_bar"` to compute limma statistics from GSVA scores and plot
#'   t-values as grouped bars.
#'
#' @param gsva_results Data frame, matrix, or named list by database. If the
#'   list comes from `run_GSVA()`, it should contain `summary` and
#'   `scores_matrix`.
#' @param database Optional database name when a list is supplied.
#' @param style Visualization style: `"heatmap"` or `"limma_bar"`.
#' @param group Optional vector of group labels (length = #samples) to build a
#'   design matrix for `style = "limma_bar"`.
#' @param design Optional design matrix (samples x groups).
#' @param contrast Optional contrast (string or matrix) for limma.
#' @param cutoff Threshold for grouping bars by t-value.
#' @param top_n Number of pathways to plot (set <= 0 for all).
#' @param strip_prefix Optional prefix to remove from pathway IDs.
#' @param output_dir Output directory for the plot.
#'
#' @return Invisibly returns a ggplot object for `style = "heatmap"`. For
#'   `style = "limma_bar"`, returns a list with `plot` and `table`.
#' @export
Viz_GSVA <- function(
    gsva_results,
    database = NULL,
    style = c("heatmap", "limma_bar"),
    group = NULL,
    design = NULL,
    contrast = NULL,
    cutoff = 1,
    top_n = 30,
    strip_prefix = NULL,
    output_dir = "gsva_viz"
) {
  style <- match.arg(style)
  if (is.list(gsva_results) && !is.data.frame(gsva_results)) {
    if (is.null(database)) {
      if (length(gsva_results) != 1) {
        stop("Provide 'database' when gsva_results has multiple entries.")
      }
      database <- names(gsva_results)[1]
    }
    gsva_results <- gsva_results[[database]]
  }
  if (is.null(database) || database == "") database <- "GSVA"
  if (is.null(gsva_results)) return(invisible(NULL))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (style == "heatmap") {
    p <- .plot_gsva_heatmap(gsva_results, database, output_dir)
    return(invisible(p))
  }

  scores_matrix <- NULL
  if (is.list(gsva_results) && !is.null(gsva_results$scores_matrix)) {
    scores_matrix <- gsva_results$scores_matrix
  } else if (is.matrix(gsva_results)) {
    scores_matrix <- gsva_results
  } else if (is.data.frame(gsva_results)) {
    scores_matrix <- as.matrix(gsva_results)
  }

  if (is.null(scores_matrix)) {
    stop("GSVA scores matrix is required for style = 'limma_bar'.")
  }

  res <- .plot_gsva_limma_bar(
    scores_matrix = scores_matrix,
    db = database,
    output_dir = output_dir,
    group = group,
    design = design,
    contrast = contrast,
    cutoff = cutoff,
    top_n = top_n,
    strip_prefix = strip_prefix
  )
  invisible(res)
}

.plot_gsva_limma_bar <- function(scores_matrix,
                                 db,
                                 output_dir,
                                 group = NULL,
                                 design = NULL,
                                 contrast = NULL,
                                 cutoff = 1,
                                 top_n = 30,
                                 strip_prefix = NULL) {
  if (is.null(design)) {
    if (is.null(group)) {
      stop("Provide 'group' or 'design' for style = 'limma_bar'.")
    }
    if (length(group) != ncol(scores_matrix)) {
      stop("'group' must have the same length as ncol(scores_matrix).")
    }
    group <- factor(group)
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    rownames(design) <- colnames(scores_matrix)
  } else if (nrow(design) != ncol(scores_matrix)) {
    stop("'design' must have nrow equal to ncol(scores_matrix).")
  }

  if (is.null(contrast)) {
    if (ncol(design) == 2) {
      contrast <- paste0(colnames(design)[2], "-", colnames(design)[1])
    } else {
      stop("Provide 'contrast' when design has more than 2 groups.")
    }
  }

  contrast_matrix <- if (is.character(contrast)) {
    limma::makeContrasts(contrasts = contrast, levels = design)
  } else {
    contrast
  }

  fit <- limma::lmFit(scores_matrix, design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  tt <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = "BH",
                        sort.by = "P")
  utils::write.csv(tt,
                   file.path(output_dir, paste0("GSVA_", db, "_limma.csv")),
                   quote = FALSE)

  if (nrow(tt) == 0) return(invisible(NULL))
  if (top_n > 0 && nrow(tt) > top_n) {
    tt <- tt[seq_len(top_n), , drop = FALSE]
  }

  df <- data.frame(ID = rownames(tt), score = tt$t, stringsAsFactors = FALSE)
  if (!is.null(strip_prefix) && nzchar(strip_prefix)) {
    df$ID <- sub(paste0("^", strip_prefix), "", df$ID)
  }
  df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),
                  labels = c("1", "2", "3"))
  df <- df[order(df$score), , drop = FALSE]
  df$ID <- factor(df$ID, levels = df$ID)
  utils::write.csv(df,
                   file.path(output_dir, paste0("GSVA_", db, "_t_values.csv")),
                   quote = FALSE, row.names = FALSE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = ID, y = score, fill = group)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("palegreen3", "snow3",
                                          "dodgerblue4"), guide = "none") +
    ggplot2::geom_hline(yintercept = c(-cutoff, cutoff),
                        color = "white", linetype = 2, linewidth = 0.3) +
    ggplot2::labs(title = paste("GSVA (limma) —", db),
                  x = NULL, y = "t value of GSVA score") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(linewidth = 0.6),
                   axis.text.y = ggplot2::element_text(size = 8))

  ggplot2::ggsave(
    file.path(output_dir, paste0("GSVA_", db, "_limma_barplot.pdf")),
    plot = p, width = 6, height = 8
  )

  list(plot = p, table = tt)
}

#' Visualize ROAST Results
#'
#' Generates a bar plot for ROAST results (Top 20 by PValue).
#'
#' @param roast_results Data frame or named list by database.
#' @param database Optional database name when a list is supplied.
#' @param output_dir Output directory for the plot.
#'
#' @return Invisibly returns the ggplot object.
#' @export
Viz_ROAST <- function(roast_results, database = NULL, output_dir = "roast_viz") {
  if (is.list(roast_results) && !is.data.frame(roast_results)) {
    if (is.null(database)) {
      if (length(roast_results) != 1) {
        stop("Provide 'database' when roast_results has multiple entries.")
      }
      database <- names(roast_results)[1]
    }
    roast_results <- roast_results[[database]]
  }
  if (is.null(database) || database == "") database <- "ROAST"
  if (is.null(roast_results)) return(invisible(NULL))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  p <- .plot_roast_bar(roast_results, database, output_dir)
  invisible(p)
}

.resolve_gsea_raw <- function(gsea_results, database = NULL) {
  if (inherits(gsea_results, "gseaResult")) return(gsea_results)

  if (is.list(gsea_results) && !is.data.frame(gsea_results)) {
    if (is.null(database)) {
      if (length(gsea_results) != 1) {
        stop("Provide 'database' when gsea_results has multiple entries.")
      }
      database <- names(gsea_results)[1]
    }
    gsea_results <- gsea_results[[database]]
  }

  if (is.null(gsea_results)) stop("No GSEA results found for visualization.")
  raw <- attr(gsea_results, "raw_result")
  if (is.null(raw)) {
    stop("Raw gseaResult object missing. Re-run run_GSEA() with updated version.")
  }
  raw
}
