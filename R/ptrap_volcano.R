#' Single volcano plot for PhosphoTRAP / TRAP-seq differential expression results
#'
#' A volcano plot puts log2 fold change (logFC, x-axis) against statistical
#' significance (-log p-value, y-axis). Genes that are both strongly enriched
#' *and* statistically reliable appear in the upper corners — the most
#' biologically meaningful candidates. Genes that fail either threshold are
#' shown in grey; those that pass both are coloured by direction (`"UP"` for
#' enriched in IP, `"DOWN"` for depleted).
#'
#' Takes a single tibble from [pTRAPPING::ptrap_de()] (for
#' `test_method = "paired.ttest"`, pass the `$results` component). Gene
#' labels are drawn only for genes listed in `genes.annot`, via
#' [ggrepel::geom_text_repel()] (static mode only). The DE classification is
#' recomputed inside this function from the supplied thresholds, so you can
#' explore different cutoffs without re-running `ptrap_de()`.
#'
#' @param de_result A tibble returned by [pTRAPPING::ptrap_de()]. For
#'   `test_method = "paired.ttest"`, pass the `$results` component.
#' @param fdr Logical. If `TRUE` (default), significance is assessed using
#'   `FDR` (false discovery rate — multiple-testing adjusted p-value); if
#'   `FALSE`, uses raw `PValue`.
#' @param lfc_threshold Minimum absolute log2 fold change for a gene to be
#'   coloured as DE. A value of `1` means at least a 2× change. Default `1`.
#' @param fdr_threshold Significance cutoff applied to the column selected by
#'   `fdr`. Default `0.05`.
#' @param log_base Numeric. Base of the logarithm for the p-value axis.
#'   Base 10 (default) is conventional; base 2 matches Tan et al. (2016).
#'   Must be a positive number other than `1`.
#' @param gene_col Column name for gene identifiers. Default `"Gene"`.
#' @param treatment_col Column name for treatment label. Default `"treatment"`.
#' @param region_col Column name for brain region label. Default `"BrainRegion"`.
#' @param colors Named vector mapping `"UP"` and `"DOWN"` to colours.
#'   Default is colourblind-friendly: `"#0072b2"` (UP), `"#E69f00"` (DOWN).
#' @param point_size Size of points. Default `3.5`.
#' @param point_alpha Opacity of highlighted points. Default `0.7`.
#' @param genes.annot Character vector of gene names to label. In static mode
#'   labels are drawn via [ggrepel::geom_text_repel()]; in interactive mode
#'   they are omitted (use the hover tooltip instead). Default `NULL`.
#' @param max_overlaps Passed to [ggrepel::geom_text_repel()]. Default `20`.
#' @param title Plot title. Auto-generated from region and treatment if `NULL`.
#' @param interactive Logical. If `TRUE`, returns an interactive
#'   [plotly::ggplotly()] object with hover tooltips showing gene name, logFC,
#'   p-value, and FDR. Gene labels (`genes.annot`) are omitted in this mode.
#'   Default `FALSE`.
#'
#' @return A [ggplot2::ggplot()] object, or a plotly object when
#'   `interactive = TRUE`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Default: -log10(FDR), static
#' ptrap_volcano(res)
#'
#' # Raw p-values, log2 scale -- as in Tan et al. 2016
#' ptrap_volcano(res, fdr = FALSE, log_base = 2)
#'
#' # Interactive with hover tooltips
#' ptrap_volcano(res, interactive = TRUE)
#' }
#'
#' @importFrom dplyr mutate case_when
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline scale_fill_manual theme_classic labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom plotly ggplotly layout toWebGL
#' @importFrom cli cli_abort

ptrap_volcano <- function(
  de_result,
  fdr = TRUE,
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  log_base = 10,
  gene_col = "Gene",
  treatment_col = "treatment",
  region_col = "BrainRegion",
  colors = NULL,
  point_size = 3.5,
  point_alpha = 0.7,
  genes.annot = NULL,
  max_overlaps = 20,
  title = NULL,
  interactive = FALSE
) {
  # --- validate log_base -------------------------------------------------
  if (
    !is.numeric(log_base) ||
      length(log_base) != 1L ||
      is.na(log_base) ||
      log_base <= 0 ||
      log_base == 1
  ) {
    cli::cli_abort(c(
      "{.arg log_base} must be a single positive number other than 1.",
      "x" = "You supplied {.val {log_base}}."
    ))
  }

  # --- extract labels ----------------------------------------------------
  treatment <- unique(de_result[[treatment_col]])
  region <- unique(de_result[[region_col]])

  # --- default colours ---------------------------------------------------
  if (is.null(colors)) {
    colors <- c("UP" = "#0072b2", "DOWN" = "#E69f00")
  }

  # --- choose p-value column ---------------------------------------------
  p_col <- if (fdr) "FDR" else "PValue"

  # --- auto-generate y-axis label ----------------------------------------
  base_int <- if (log_base == as.integer(log_base)) {
    as.integer(log_base)
  } else {
    log_base
  }
  # plotly cannot parse bquote() expressions; use plain text for interactive
  y_label <- if (interactive) {
    if (fdr) paste0("-log", base_int, "(FDR)") else paste0("-log", base_int, "(p-value)")
  } else if (fdr) {
    bquote(-log[.(base_int)](FDR))
  } else {
    bquote(-log[.(base_int)]("p-value"))
  }

  # --- recompute DE classification and build hover text ------------------
  plot_data <- de_result |>
    mutate(
      DE = case_when(
        .data$logFC > lfc_threshold &
          !is.na(.data[[p_col]]) &
          .data[[p_col]] < fdr_threshold ~ "UP",
        .data$logFC < -lfc_threshold &
          !is.na(.data[[p_col]]) &
          .data[[p_col]] < fdr_threshold ~ "DOWN",
        TRUE ~ "Not DE"
      ),
      y_val = .data[[p_col]],
      text_hover = paste0(
        "<b>", .data[[gene_col]], "</b><br>",
        "logFC: ", round(.data$logFC, 3), "<br>",
        "PValue: ", signif(.data$PValue, 3), "<br>",
        "FDR: ", signif(.data$FDR, 3)
      )
    )

  # --- annotation data ---------------------------------------------------
  if (!is.null(genes.annot)) {
    bad <- setdiff(genes.annot, plot_data[[gene_col]])
    if (length(bad) > 0) {
      stop(
        "The following names in 'genes.annot' were not found in the '",
        gene_col,
        "' column: ",
        paste(bad, collapse = ", ")
      )
    }
    annot_data <- plot_data[plot_data[[gene_col]] %in% genes.annot, ]
  } else {
    annot_data <- plot_data[integer(0), ]
  }

  # --- default title -----------------------------------------------------
  if (is.null(title)) {
    title <- paste0(region, " (", treatment, ")")
  }

  # --- build plot --------------------------------------------------------
  p <- ggplot(
    plot_data,
    aes(
      x = .data$logFC,
      y = -log(.data$y_val, base = log_base),
      text = .data$text_hover
    )
  ) +
    geom_point(
      data = plot_data[plot_data$DE == "Not DE", ],
      color = "grey80",
      alpha = 0.5,
      size = point_size
    ) +
    geom_point(
      data = plot_data[plot_data$DE != "Not DE", ],
      aes(fill = .data$DE),
      alpha = point_alpha,
      size = point_size,
      shape = 21
    ) +
    geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log(fdr_threshold, base = log_base),
      linetype = "dashed"
    ) +
    scale_fill_manual(values = colors) +
    geom_point(
      data = annot_data,
      aes(fill = .data$DE),
      size = point_size + 0.3,
      shape = 21,
      stroke = 0.5,
      color = "black"
    ) +
    theme_classic() +
    labs(
      x = if (interactive) "log2 Fold Change (IP / Input)" else expression(log[2] ~ "Fold Change (IP / Input)"),
      y = y_label,
      title = title,
      fill = "DE"
    )

  if (!interactive) {
    p <- p + geom_text_repel(
      data = annot_data,
      aes(label = .data[[gene_col]]),
      size = 4.3,
      color = "black",
      max.overlaps = max_overlaps,
      min.segment.length = 0
    )
  }

  if (interactive) {
    plt <- plotly::ggplotly(p, tooltip = "text", width = NULL, height = 500) |>
      plotly::layout(dragmode = "zoom", autosize = TRUE)
    # scattergl (used by toWebGL) does not support 'hoveron'; remove it first
    plt$x$data <- lapply(plt$x$data, function(tr) { tr$hoveron <- NULL; tr })
    plotly::toWebGL(plt)
  } else {
    p
  }
}
