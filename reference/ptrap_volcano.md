# Single volcano plot for TRAP-seq or RNA-seq differential expression results

Takes a single tibble produced by
[`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md),
[`limma_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/limma_de.md)
or similar and returns a classic volcano plot (logFC on the x-axis,
-log10(FDR) on the y-axis). Non-significant genes are shown in grey;
genes classified as `"UP"` or `"DOWN"` in the `diffexpressed` column are
highlighted with distinct fill colours and labelled with their gene
identifier. Dashed threshold lines are drawn at ±`lfc_threshold`
(vertical) and `-log10(fdr_threshold)` (horizontal).

## Usage

``` r
ptrap_volcano(
  de_result,
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  gene_col = "Gene",
  treatment_col = "Treatment",
  region_col = "BrainRegion",
  colors = NULL,
  point_size = 3.5,
  point_alpha = 0.7,
  max_overlaps = 20,
  title = NULL
)
```

## Arguments

- de_result:

  A tibble returned by
  [`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md),
  containing at minimum columns `Gene` (or `gene_col`), `logFC`, `FDR`,
  and `diffexpressed`.

- lfc_threshold:

  Minimum absolute log2 fold change used to draw the vertical threshold
  lines. Should match the value used in
  [`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md).
  Default is `1`.

- fdr_threshold:

  FDR cutoff used to draw the horizontal threshold line (plotted as
  `-log10(fdr_threshold)`). Should match the value used in
  [`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md).
  Default is `0.05`.

- gene_col:

  Name of the column containing gene identifiers. Default is `"Gene"`.

- treatment_col:

  Name of the column containing the treatment label. Used to build the
  default plot title. Default is `"Treatment"`.

- region_col:

  Name of the column containing the brain region label. Used to build
  the default plot title. Default is `"BrainRegion"`.

- colors:

  A named character vector mapping `"UP"` and `"DOWN"` to colours. If
  `NULL` (default), a colourblind-friendly palette is used (`"#D55E00"`
  for `"UP"`, `"#0072B2"` for `"DOWN"`).

- point_size:

  Size of the points. Default is `3.5`.

- point_alpha:

  Opacity of the highlighted (significant) points. Default is `0.7`.

- max_overlaps:

  Passed to
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Controls how many overlapping labels are suppressed. Default is `20`.

- title:

  Plot title. If `NULL` (default), a title is auto-generated from the
  brain region and treatment labels found in `de_result`.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  region_name    = "POA",
  treatment_name = "pb"
)

# Default plot
ptrap_volcano(res)

# Custom thresholds and colours
ptrap_volcano(res,
              lfc_threshold = 0.5,
              fdr_threshold = 0.1,
              colors        = c("UP" = "#E69F00", "DOWN" = "#56B4E9"))
} # }
```
