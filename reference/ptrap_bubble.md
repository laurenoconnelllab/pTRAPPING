# Bubble plot of differential expression across treatments

Displays genes of interest (y-axis) against treatment conditions
(x-axis) as a bubble plot. Circle **fill** encodes the log2 fold change
(logFC) via a continuous colour gradient, and circle **size** encodes
statistical significance as \\-\log\_{10}(\mathrm{FDR})\\ (default) or
\\-\log\_{10}(p\text{-value})\\.

## Usage

``` r
ptrap_bubble(
  de_result,
  sig_size = c("FDR", "PValue"),
  colors_lfc = c("#25599b", "#ABD0DD", "#F2F9FE", "#F88705", "#B8351F"),
  gene_col = "Gene",
  treatment_col = "treatment",
  size_range = c(2, 12),
  point_alpha = 0.9,
  title = NULL,
  interactive = FALSE
)
```

## Arguments

- de_result:

  A data frame of combined results from one or more calls to
  [`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md).
  For `test_method = "paired.ttest"`, use the `$results` component of
  each and combine them with
  [`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html).
  Must contain columns for gene names, treatment, `logFC`, `FDR`, and
  `PValue`.

- sig_size:

  Character. Which column to use for the size aesthetic (displayed as
  \\-\log\_{10}\\). Either `"FDR"` (default) or `"PValue"`.

- colors_lfc:

  Character vector of at least three colours defining the fill gradient
  mapped to logFC values. Passed to
  [`ggplot2::scale_fill_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html).
  Default: `c("#25599b", "#ABD0DD", "#F2F9FE", "#F88705", "#B8351F")`.

- gene_col:

  Name of the column containing gene identifiers. Default `"Gene"`.

- treatment_col:

  Name of the column containing the treatment label. Default
  `"treatment"`.

- size_range:

  Numeric vector of length two controlling the minimum and maximum point
  sizes. Passed to
  [`ggplot2::scale_size_continuous()`](https://ggplot2.tidyverse.org/reference/scale_size.html).
  Default `c(2, 12)`.

- point_alpha:

  Opacity of the bubbles. Default `0.9`.

- title:

  Optional plot title. If `NULL` (default), no title is added.

- interactive:

  Logical. If `TRUE`, returns an interactive
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)
  object with hover tooltips showing gene name, logFC, p-value, and FDR.
  Default `FALSE`.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object, or a plotly object when `interactive = TRUE`.

## Details

This is particularly useful after running
[`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md)
for multiple treatment conditions (e.g., PACAP and BDNF) on the same set
of genes: bind the `$results` tibbles together with
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
and pass the combined data frame to `ptrap_bubble()`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run DE for two treatments on the same gene set
pacap_de <- ptrap_de(counts, ..., treatment_name = "PACAP",
                     genes.filter = my_genes)
bdnf_de  <- ptrap_de(counts, ..., treatment_name = "BDNF",
                     genes.filter = my_genes)

# Combine results and plot
library(dplyr)
bind_rows(pacap_de$results, bdnf_de$results) |>
  ptrap_bubble()

# Use raw p-values for size and a custom colour palette
bind_rows(pacap_de$results, bdnf_de$results) |>
  ptrap_bubble(sig_size = "PValue",
               colors_lfc = c("blue", "white", "red"))

# Interactive with hover tooltips
bind_rows(pacap_de$results, bdnf_de$results) |>
  ptrap_bubble(interactive = TRUE)
} # }
```
