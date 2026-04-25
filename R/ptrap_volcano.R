#' Single volcano plot for TRAP-seq or RNA-seq differential expression results
#'
#' Takes a single tibble produced by [pTRAPPING::ptrap_de()]
#' and returns a classic volcano plot (logFC on the x-axis,
#' \eqn{-\log_{10}(\text{FDR})} or \eqn{-\log_{10}(p\text{-value})} on the
#' y-axis). Non-significant genes are shown in grey; genes classified as
#' `"UP"` or `"DOWN"` are highlighted with distinct fill colours. Gene labels
#' are added via [ggrepel::geom_text_repel()] **only for genes supplied in
#' `genes.annot`**; by default (`genes.annot = NULL`) no labels are drawn,
#' keeping the plot uncluttered. Dashed threshold lines are drawn at
#' ±`lfc_threshold` (vertical) and at the chosen p-value cutoff (horizontal).
#'
#' The DE classification used for point colouring is computed **inside this
#' function** from `logFC`, the chosen p-value column (`FDR` or `PValue`), and
#' the supplied thresholds — so the plot always reflects the thresholds you
#' pass here, regardless of what was used in [pTRAPPING::ptrap_de()].
#'
#' @param de_result A tibble returned by [pTRAPPING::ptrap_de()], containing at
#'   minimum the columns `Gene` (or `gene_col`), `logFC`, `FDR`, and `PValue`.
#'   For `test_method = "paired.ttest"`, pass the `$results` component of the
#'   returned list.
#' @param fdr Logical. If `TRUE` (default), the y-axis shows
#'   \eqn{-\log_{10}(\text{FDR})} and significance is assessed against the
#'   BH-adjusted p-value (`FDR` column). If `FALSE`, the y-axis shows
#'   \eqn{-\log_{10}(p\text{-value})} and significance is assessed against the
#'   raw p-value (`PValue` column).
#' @param lfc_threshold Minimum absolute log2 fold change required to classify
#'   a gene as differentially expressed (also sets the vertical threshold lines).
#'   Default is `1`.
#' @param fdr_threshold P-value cutoff used to draw the horizontal threshold
#'   line and to classify genes as DE. Applied to `FDR` when `fdr = TRUE` and
#'   to `PValue` when `fdr = FALSE`. Default is `0.05`.
#' @param gene_col Name of the column containing gene identifiers. Default is
#'   `"Gene"`.
#' @param treatment_col Name of the column containing the treatment label. Used
#'   to build the default plot title. Default is `"Treatment"`.
#' @param region_col Name of the column containing the brain region label. Used
#'   to build the default plot title. Default is `"BrainRegion"`.
#' @param colors A named character vector mapping `"UP"` and `"DOWN"` to
#'   colours. If `NULL` (default), a colourblind-friendly palette is used
#'   (`"#D55E00"` for `"UP"`, `"#0072B2"` for `"DOWN"`).
#' @param point_size Size of the points. Default is `3.5`.
#' @param point_alpha Opacity of the highlighted (significant) points.
#'   Default is `0.7`.
#' @param genes.annot Character vector of gene names to label on the plot via
#'   [ggrepel::geom_text_repel()]. Names must match values in the column
#'   specified by `gene_col`; an error is raised if any names are not found.
#'   Default is `NULL` (no labels).
#' @param max_overlaps Passed to [ggrepel::geom_text_repel()]. Controls how
#'   many overlapping labels are suppressed. Default is `20`.
#' @param title Plot title. If `NULL` (default), a title is auto-generated from
#'   the brain region and treatment labels found in `de_result`.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb"
#' )
#'
#' # Default plot — y-axis uses FDR-adjusted p-values
#' ptrap_volcano(res)
#'
#' # Use raw (unadjusted) p-values on the y-axis
#' ptrap_volcano(res, fdr = FALSE)
#'
#' # Custom thresholds and colours
#' ptrap_volcano(res,
#'               lfc_threshold = 0.5,
#'               fdr_threshold = 0.1,
#'               colors        = c("UP" = "#E69F00", "DOWN" = "#56B4E9"))
#' }
#'
#' @importFrom dplyr mutate case_when
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline
#'   scale_fill_manual theme_classic labs
#' @importFrom ggrepel geom_text_repel

ptrap_volcano <- function(
  de_result,
  fdr = TRUE,
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  gene_col = "Gene",
  treatment_col = "Treatment",
  region_col = "BrainRegion",
  colors = NULL,
  point_size = 3.5,
  point_alpha = 0.7,
  genes.annot = NULL,
  max_overlaps = 20,
  title = NULL
) {
  # --- extract labels from the result table ------------------------------------
  treatment <- unique(de_result[[treatment_col]])
  region    <- unique(de_result[[region_col]])

  # --- default colourblind-friendly palette ------------------------------------
  if (is.null(colors)) {
    colors <- c("UP" = "#D55E00", "DOWN" = "#0072B2")
  }

  # --- choose p-value column and y-axis label based on fdr argument -----------
  p_col   <- if (fdr) "FDR" else "PValue"
  y_label <- if (fdr) {
    expression(-log[10](FDR))
  } else {
    expression(-log[10]("p-value"))
  }

  # --- recompute DE classification using the chosen p-value column ------------
  # This ensures the coloring always reflects the thresholds supplied HERE,
  # independent of what was used inside ptrap_de().
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
      y_val = .data[[p_col]]
    )

  # --- validate and prepare annotation data ------------------------------------
  if (!is.null(genes.annot)) {
    bad <- setdiff(genes.annot, plot_data[[gene_col]])
    if (length(bad) > 0) {
      stop(
        "The following names in 'genes.annot' were not found in the '",
        gene_col, "' column: ",
        paste(bad, collapse = ", ")
      )
    }
    annot_data <- plot_data[plot_data[[gene_col]] %in% genes.annot, ]
  } else {
    annot_data <- plot_data[integer(0), ]   # zero-row tibble, correct columns
  }

  # --- default title -----------------------------------------------------------
  if (is.null(title)) {
    title <- paste0(region, " (", treatment, ")")
  }

  # --- build plot --------------------------------------------------------------
  ggplot(plot_data, aes(x = .data$logFC, y = -log10(.data$y_val))) +
    geom_point(
      data  = plot_data[plot_data$DE == "Not DE", ],
      color = "grey80",
      alpha = 0.5,
      size  = point_size
    ) +
    geom_point(
      data  = plot_data[plot_data$DE != "Not DE", ],
      aes(fill = .data$DE),
      alpha = point_alpha,
      size  = point_size,
      shape = 21
    ) +
    geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype   = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype   = "dashed"
    ) +
    scale_fill_manual(values = colors) +
    geom_point(
      data   = annot_data,
      aes(fill = .data$DE),
      size   = point_size + 0.3,
      shape  = 21,
      stroke = 0.5,
      color  = "black"
    ) +
    geom_text_repel(
      data             = annot_data,
      aes(label        = .data[[gene_col]]),
      size             = 4.3,
      color            = "black",
      max.overlaps     = max_overlaps,
      min.segment.length = 0
    ) +
    theme_classic() +
    labs(
      x     = expression(log[2] ~ "Fold Change (IP / Input)"),
      y     = y_label,
      title = title,
      fill  = "DE"
    )
}
