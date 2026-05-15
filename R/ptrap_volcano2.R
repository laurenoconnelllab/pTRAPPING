#' Compare IP enrichment across two treatment conditions in a single scatter plot
#'
#' When you have run [pTRAPPING::ptrap_de()] for two treatment groups (e.g.,
#' pair-bonded and non-bonded animals) in the same brain region, looking at
#' two separate volcano plots makes it hard to see which genes behave similarly
#' or differently across conditions. This function places both results on a
#' single 2D scatter: the logFC (IP vs. input) of condition 1 on the y-axis
#' and condition 2 on the x-axis. Each point is one gene.
#'
#' Reading the plot:
#' * Genes near the **diagonal** (y = x dotted line) are equally enriched in
#'   both conditions.
#' * Genes **above** the diagonal are more enriched in condition 1; genes
#'   **below** are more enriched in condition 2.
#' * Genes that pass the fold-change and p-value thresholds in **both**
#'   conditions, or in only **one**, are highlighted in distinct colours.
#'   Genes that fail in both are shown in grey.
#'
#' Significance uses FDR-adjusted p-values by default (`fdr = TRUE`); set
#' `fdr = FALSE` to use raw p-values instead. The DE classification is
#' recomputed inside this function from the supplied thresholds, so you can
#' explore different cutoffs without re-running `ptrap_de()`.
#'
#' @param de_result_1 A tibble returned by [pTRAPPING::ptrap_de()] for the
#'   **first** treatment condition (plotted on the **y-axis**).
#' @param de_result_2 A tibble returned by [pTRAPPING::ptrap_de()] for the
#'   **second** treatment condition (plotted on the **x-axis**).
#' @param fdr Logical. If `TRUE` (default), significance is assessed using the
#'   `FDR` column (BH-adjusted p-values). If `FALSE`, the raw `PValue` column
#'   is used instead.
#' @param lfc_threshold Minimum absolute log2 fold change required to classify
#'   a gene as differentially expressed in a given condition (also sets the
#'   vertical and horizontal threshold lines). Default is `1`.
#' @param fdr_threshold P-value cutoff used to classify genes as DE. Applied
#'   to `FDR` when `fdr = TRUE` and to `PValue` when `fdr = FALSE`. Default
#'   is `0.05`.
#' @param gene_col Name of the column containing gene identifiers in both
#'   result tibbles. Default is `"Gene"`.
#' @param treatment_col Name of the column containing the treatment label in
#'   both result tibbles. Used to auto-generate axis labels and DE class names.
#'   Default is `"treatment"`.
#' @param region_col Name of the column containing the brain region label.
#'   Used to set the default plot title. Default is `"BrainRegion"`.
#' @param colors A named character vector mapping each DE class to a colour.
#'   Names must be `"DE in both"`, `"DE only <t1>"`, and `"DE only <t2>"`,
#'   where `<t1>` / `<t2>` are the treatment names found in the data (e.g.,
#'   `"DE only pb"` and `"DE only sol"`). See the examples for a template.
#'   If `NULL` (default), a colourblind-friendly palette is used.
#' @param point_size Size of the points. Default is `3.5`.
#' @param point_alpha Opacity of the coloured (significant) points.
#'   Default is `0.7`.
#' @param genes.annot Character vector of gene names to label on the plot via
#'   [ggrepel::geom_text_repel()] (static mode only). Names must match values
#'   in the column specified by `gene_col`; an error is raised if any names
#'   are not found. Default is `NULL` (no labels).
#' @param max_overlaps Passed to [ggrepel::geom_text_repel()]. Controls how
#'   many overlapping labels are allowed before they are hidden. Default
#'   is `20`.
#' @param title Plot title. If `NULL` (default), the brain region name
#'   extracted from `de_result_1` is used.
#' @param interactive Logical. If `TRUE`, returns an interactive
#'   [plotly::ggplotly()] object with hover tooltips showing gene name, logFC
#'   for each condition, p-value, and FDR. Gene labels (`genes.annot`) are
#'   omitted in this mode. Default `FALSE`.
#'
#' @return A [ggplot2::ggplot()] object, or a plotly object when
#'   `interactive = TRUE`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res_pb  <- ptrap_de(counts_mat, sample_df, gene_ids,
#'                     region_name = "POA", treatment_name = "pb")
#' res_sol <- ptrap_de(counts_mat, sample_df, gene_ids,
#'                     region_name = "POA", treatment_name = "sol")
#'
#' # Default plot — FDR on both axes
#' ptrap_volcano2(res_pb, res_sol)
#'
#' # Use raw p-values for the DE classification
#' ptrap_volcano2(res_pb, res_sol, fdr = FALSE)
#'
#' # Interactive with hover tooltips
#' ptrap_volcano2(res_pb, res_sol, interactive = TRUE)
#'
#' # Custom thresholds, title and colours
#' ptrap_volcano2(res_pb, res_sol,
#'                lfc_threshold = 0.5,
#'                fdr_threshold = 0.1,
#'                title         = "POA: pb vs sol",
#'                colors        = c("DE in both"  = "#D55E00",
#'                                  "DE only pb"  = "#0072B2",
#'                                  "DE only sol" = "#49a15f"))
#' }
#'
#' @importFrom dplyr select left_join mutate case_when arrange filter all_of
#' @importFrom rlang .data :=
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_vline geom_hline scale_fill_manual theme_classic labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom plotly ggplotly layout toWebGL

ptrap_volcano2 <- function(
  de_result_1,
  de_result_2,
  fdr = TRUE,
  lfc_threshold = 1,
  fdr_threshold = 0.05,
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
  # --- extract treatment and region names from the result tables --------------
  t1 <- unique(de_result_1[[treatment_col]])
  t2 <- unique(de_result_2[[treatment_col]])
  region <- unique(de_result_1[[region_col]])

  # --- DE class labels (derived from treatment names) -------------------------
  class_both <- "DE in both"
  class_only1 <- paste0("DE only ", t1)
  class_only2 <- paste0("DE only ", t2)
  class_none <- "Not DE"

  # --- default colourblind-friendly palette -----------------------------------
  if (is.null(colors)) {
    colors <- c("#0072b2", "#E69f00", "#009e73")
    names(colors) <- c(class_both, class_only1, class_only2)
  }

  # --- which p-value column drives DE classification --------------------------
  stat_col_1 <- if (fdr) "fdr_1" else "pval_1"
  stat_col_2 <- if (fdr) "fdr_2" else "pval_2"

  # --- join both conditions, keeping PValue and FDR for hover -----------------
  combined <- de_result_1 |>
    select(
      all_of(gene_col),
      logFC_1 = "logFC",
      pval_1 = "PValue",
      fdr_1 = "FDR"
    ) |>
    left_join(
      de_result_2 |>
        select(
          all_of(gene_col),
          logFC_2 = "logFC",
          pval_2 = "PValue",
          fdr_2 = "FDR"
        ),
      by = gene_col
    )

  # --- drop genes missing a logFC in either condition ------------------------
  # This happens when filterByExpr removes a gene in one condition but not the
  # other, leaving NAs after the left_join. Those genes cannot be placed on a
  # 2D scatter and are excluded from the plot with an informative message.
  n_missing <- sum(is.na(combined$logFC_1) | is.na(combined$logFC_2))
  if (n_missing > 0) {
    message(
      n_missing,
      " gene(s) were dropped because they were only tested in one ",
      "condition (likely filtered out by filterByExpr in the other). ",
      "They cannot be placed on a 2D scatter plot."
    )
    combined <- combined |>
      filter(!is.na(.data$logFC_1), !is.na(.data$logFC_2))
  }

  # --- classify genes and build hover text ------------------------------------
  combined <- combined |>
    mutate(
      sig_1 = !is.na(.data[[stat_col_1]]) &
        .data[[stat_col_1]] < fdr_threshold &
        abs(.data$logFC_1) > lfc_threshold,
      sig_2 = !is.na(.data[[stat_col_2]]) &
        .data[[stat_col_2]] < fdr_threshold &
        abs(.data$logFC_2) > lfc_threshold
    ) |>
    mutate(
      DE_class = case_when(
        .data$sig_1 & .data$sig_2 ~ class_both,
        .data$sig_1 & !.data$sig_2 ~ class_only1,
        !.data$sig_1 & .data$sig_2 ~ class_only2,
        TRUE ~ class_none
      ),
      text_hover = paste0(
        "<b>", .data[[gene_col]], "</b><br>",
        "logFC (", t1, "): ", round(.data$logFC_1, 3), "<br>",
        "logFC (", t2, "): ", round(.data$logFC_2, 3), "<br>",
        "PValue (", t1, "): ", signif(.data$pval_1, 3), "<br>",
        "FDR (", t1, "): ", signif(.data$fdr_1, 3), "<br>",
        "PValue (", t2, "): ", signif(.data$pval_2, 3), "<br>",
        "FDR (", t2, "): ", signif(.data$fdr_2, 3)
      )
    ) |>
    arrange(.data$DE_class)

  # --- validate and prepare annotation data ------------------------------------
  if (!is.null(genes.annot)) {
    bad <- setdiff(genes.annot, combined[[gene_col]])
    if (length(bad) > 0) {
      stop(
        "The following names in 'genes.annot' were not found in the '",
        gene_col,
        "' column: ",
        paste(bad, collapse = ", ")
      )
    }
    annot_data <- combined[combined[[gene_col]] %in% genes.annot, ]
  } else {
    annot_data <- combined[integer(0), ] # zero-row data frame, correct columns
  }

  # --- default title ----------------------------------------------------------
  if (is.null(title)) {
    title <- region
  }

  # --- build plot -------------------------------------------------------------
  p <- ggplot(
    combined,
    aes(
      x = .data$logFC_2,
      y = .data$logFC_1,
      text = .data$text_hover
    )
  ) +
    geom_point(
      data = combined[combined$DE_class == class_none, ],
      color = "grey80",
      alpha = 0.5,
      size = point_size
    ) +
    geom_point(
      data = combined[combined$DE_class != class_none, ],
      aes(fill = .data$DE_class),
      alpha = point_alpha,
      size = point_size,
      shape = 21
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = c(-lfc_threshold, lfc_threshold),
      linetype = "dashed"
    ) +
    scale_fill_manual(values = colors) +
    geom_point(
      data = annot_data,
      aes(fill = .data$DE_class),
      size = point_size + 0.3,
      shape = 21,
      stroke = 0.5,
      color = "black"
    ) +
    theme_classic() +
    labs(
      x = paste0("logFC IP/Input (", t2, ")"),
      y = paste0("logFC IP/Input (", t1, ")"),
      title = title,
      fill = "DE class"
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
