#' Single volcano plot for TRAP-seq or RNA-seq differential expression results
#'
#' Takes a single tibble produced by [pTRAPPING::ptrap_de()]
#' and returns a classic volcano plot (logFC on the x-axis,
#' -log_base(p) on the y-axis). Non-significant genes are shown in grey;
#' genes classified as `"UP"` or `"DOWN"` are highlighted with distinct
#' fill colours. Gene labels are added via [ggrepel::geom_text_repel()]
#' **only for genes supplied in `genes.annot`**.
#'
#' The DE classification is recomputed inside this function from `logFC`,
#' the chosen p-value column, and the supplied thresholds.
#'
#' @param de_result A tibble returned by [pTRAPPING::ptrap_de()]. For
#'   `test_method = "paired.ttest"`, pass the `$results` component.
#' @param fdr Logical. If `TRUE` (default), uses `FDR`; if `FALSE`, uses `PValue`.
#' @param lfc_threshold Minimum absolute log2 fold change. Default `1`.
#' @param fdr_threshold P-value cutoff. Default `0.05`.
#' @param log_base Numeric. Base of the logarithm for the p-value axis.
#'   Default is `10` (-log10, current behaviour). Use `2` for -log2(p)
#'   as in Tan et al. (2016). Must be a positive number other than `1`.
#' @param gene_col Column name for gene identifiers. Default `"Gene"`.
#' @param treatment_col Column name for treatment label. Default `"Treatment"`.
#' @param region_col Column name for brain region label. Default `"BrainRegion"`.
#' @param colors Named vector mapping `"UP"` and `"DOWN"` to colours.
#'   Default is colourblind-friendly: `"#D55E00"` (UP), `"#0072B2"` (DOWN).
#' @param point_size Size of points. Default `3.5`.
#' @param point_alpha Opacity of highlighted points. Default `0.7`.
#' @param genes.annot Character vector of gene names to label. Default `NULL`.
#' @param max_overlaps Passed to [ggrepel::geom_text_repel()]. Default `20`.
#' @param title Plot title. Auto-generated from region and treatment if `NULL`.
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Default: -log10(FDR)
#' ptrap_volcano(res)
#'
#' # Raw p-values, log2 scale -- as in Tan et al. 2016
#' ptrap_volcano(res, fdr = FALSE, log_base = 2)
#' }
#'
#' @importFrom dplyr mutate case_when
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline scale_fill_manual theme_classic labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom cli cli_abort

ptrap_volcano <- function(
  de_result,
  fdr = TRUE,
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  log_base = 10,
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
    colors <- c("UP" = "#D55E00", "DOWN" = "#0072B2")
  }

  # --- choose p-value column ---------------------------------------------
  p_col <- if (fdr) "FDR" else "PValue"

  # --- auto-generate y-axis label ----------------------------------------
  base_int <- if (log_base == as.integer(log_base)) {
    as.integer(log_base)
  } else {
    log_base
  }
  y_label <- if (fdr) {
    bquote(-log[.(base_int)](FDR))
  } else {
    bquote(-log[.(base_int)]("p-value"))
  }

  # --- recompute DE classification ---------------------------------------
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
  ggplot(
    plot_data,
    aes(x = .data$logFC, y = -log(.data$y_val, base = log_base))
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
    geom_text_repel(
      data = annot_data,
      aes(label = .data[[gene_col]]),
      size = 4.3,
      color = "black",
      max.overlaps = max_overlaps,
      min.segment.length = 0
    ) +
    theme_classic() +
    labs(
      x = expression(log[2] ~ "Fold Change (IP / Input)"),
      y = y_label,
      title = title,
      fill = "DE"
    )
}
