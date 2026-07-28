#' Bubble plot of differential expression across treatments
#'
#' Displays genes of interest (y-axis) against treatment conditions (x-axis) as
#' a bubble plot. Circle **fill** encodes the log2 fold change (logFC) via a
#' continuous colour gradient, and circle **size** encodes statistical
#' significance as \eqn{-\log_{10}(\mathrm{FDR})} (default) or
#' \eqn{-\log_{10}(p\text{-value})}.
#'
#' This is particularly useful after running [pTRAPPING::ptrap_de()] for
#' multiple treatment conditions (e.g., PACAP and BDNF) on the same set of
#' genes: bind the `$results` tibbles together with [dplyr::bind_rows()] and
#' pass the combined data frame to `ptrap_bubble()`.
#'
#' @param de_result A data frame of combined results from one or more calls to
#'   [pTRAPPING::ptrap_de()]. For `test_method = "paired.ttest"`, use the
#'   `$results` component of each and combine them with [dplyr::bind_rows()].
#'   Must contain columns for gene names, treatment, `logFC`, `FDR`, and
#'   `PValue`.
#' @param sig_size Character. Which column to use for the size aesthetic
#'   (displayed as \eqn{-\log_{10}}). Either `"FDR"` (default) or `"PValue"`.
#' @param colors_lfc Character vector of at least three colours defining the
#'   fill gradient mapped to logFC values. Passed to
#'   [ggplot2::scale_fill_gradientn()]. Default:
#'   `c("#25599b", "#ABD0DD", "#F2F9FE", "#F88705", "#B8351F")`.
#' @param gene_col Name of the column containing gene identifiers.
#'   Default `"Gene"`.
#' @param treatment_col Name of the column containing the treatment label.
#'   Default `"treatment"`.
#' @param size_range Numeric vector of length two controlling the minimum and
#'   maximum point sizes. Passed to [ggplot2::scale_size_continuous()].
#'   Default `c(2, 12)`.
#' @param point_alpha Opacity of the bubbles. Default `0.9`.
#' @param title Optional plot title. If `NULL` (default), no title is added.
#' @param interactive Logical. If `TRUE`, returns an interactive
#'   [plotly::ggplotly()] object with hover tooltips showing gene name, logFC,
#'   p-value, and FDR. Default `FALSE`.
#'
#' @return A [ggplot2::ggplot()] object, or a plotly object when
#'   `interactive = TRUE`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run DE for two treatments on the same gene set
#' pacap_de <- ptrap_de(counts, ..., treatment_name = "PACAP",
#'                      genes.filter = my_genes)
#' bdnf_de  <- ptrap_de(counts, ..., treatment_name = "BDNF",
#'                      genes.filter = my_genes)
#'
#' # Combine results and plot
#' library(dplyr)
#' bind_rows(pacap_de$results, bdnf_de$results) |>
#'   ptrap_bubble()
#'
#' # Use raw p-values for size and a custom colour palette
#' bind_rows(pacap_de$results, bdnf_de$results) |>
#'   ptrap_bubble(sig_size = "PValue",
#'                colors_lfc = c("blue", "white", "red"))
#'
#' # Interactive with hover tooltips
#' bind_rows(pacap_de$results, bdnf_de$results) |>
#'   ptrap_bubble(interactive = TRUE)
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous scale_fill_gradientn theme_bw theme element_text element_line element_blank labs
#' @importFrom plotly ggplotly layout toWebGL
#' @importFrom cli cli_abort

ptrap_bubble <- function(
  de_result,
  sig_size = c("FDR", "PValue"),
  colors_lfc = c("#25599b", "#ABD0DD", "#F2F9FE", "#F88705", "#B8351F"),
  gene_col = "Gene",
  treatment_col = "treatment",
  size_range = c(2, 12),
  point_alpha = 0.9,
  title = NULL,

  interactive = FALSE
) {
  # --- validate sig_size ----------------------------------------------------
  sig_size <- match.arg(sig_size)

  # --- validate colors_lfc -------------------------------------------------
  if (!is.character(colors_lfc) || length(colors_lfc) < 3L) {
    cli::cli_abort(c(
      "{.arg colors_lfc} must be a character vector with at least 3 colours.",
      "x" = "You supplied a vector of length {.val {length(colors_lfc)}}."
    ))
  }

  # --- validate required columns -------------------------------------------
  required_cols <- c(gene_col, treatment_col, "logFC", "FDR", "PValue")
  missing_cols <- setdiff(required_cols, names(de_result))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(c(
      "Column(s) not found in {.arg de_result}:",
      "x" = "{.val {missing_cols}}"
    ))
  }

  # --- size label -----------------------------------------------------------
  size_label <- if (interactive) {
    paste0("-log10(", sig_size, ")")
  } else {
    if (sig_size == "FDR") {
      expression(-log[10](FDR))
    } else {
      expression(-log[10]("p-value"))
    }
  }

  # --- build hover text -----------------------------------------------------
  plot_data <- de_result |>
    mutate(
      size_val = -log10(.data[[sig_size]]),
      text_hover = paste0(
        "<b>", .data[[gene_col]], "</b><br>",
        "logFC: ", round(.data$logFC, 3), "<br>",
        "PValue: ", signif(.data$PValue, 3), "<br>",
        "FDR: ", signif(.data$FDR, 3)
      )
    )

  # --- build plot -----------------------------------------------------------
  p <- ggplot(
    plot_data,
    aes(
      x = .data[[treatment_col]],
      y = .data[[gene_col]],
      size = .data$size_val,
      fill = .data$logFC,
      text = .data$text_hover
    )
  ) +
    geom_point(
      shape = 21,
      color = "black",
      alpha = point_alpha
    ) +
    scale_size_continuous(
      name = size_label,
      range = size_range
    ) +
    scale_fill_gradientn(
      colours = colors_lfc,
      name = "logFC"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = "Treatment",
      y = gene_col
    )

  if (!is.null(title)) {
    p <- p + labs(title = title)
  }

  # --- interactive mode -----------------------------------------------------
  if (interactive) {
    plt <- plotly::ggplotly(p, tooltip = "text") |>
      plotly::layout(dragmode = "zoom", autosize = TRUE)
    plt$x$data <- lapply(plt$x$data, function(tr) {
      tr$hoveron <- NULL
      tr
    })
    plotly::toWebGL(plt)
  } else {
    p
  }
}
