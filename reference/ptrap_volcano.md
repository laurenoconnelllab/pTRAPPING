# Single volcano plot for TRAP-seq or RNA-seq differential expression results

Takes a single tibble produced by
[`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md)
and returns a classic volcano plot (logFC on the x-axis,
\\-\log\_{10}(\text{FDR})\\ or \\-\log\_{10}(p\text{-value})\\ on the
y-axis). Non-significant genes are shown in grey; genes classified as
`"UP"` or `"DOWN"` are highlighted with distinct fill colours. Gene
labels are added via
[`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
**only for genes supplied in `genes.annot`**; by default
(`genes.annot = NULL`) no labels are drawn, keeping the plot
uncluttered. Dashed threshold lines are drawn at ±`lfc_threshold`
(vertical) and at the chosen p-value cutoff (horizontal).

## Usage

``` r
ptrap_volcano(
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
)
```

## Arguments

- de_result:

  A tibble returned by
  [`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md),
  containing at minimum the columns `Gene` (or `gene_col`), `logFC`,
  `FDR`, and `PValue`. For `test_method = "paired.ttest"`, pass the
  `$results` component of the returned list.

- fdr:

  Logical. If `TRUE` (default), the y-axis shows
  \\-\log\_{10}(\text{FDR})\\ and significance is assessed against the
  BH-adjusted p-value (`FDR` column). If `FALSE`, the y-axis shows
  \\-\log\_{10}(p\text{-value})\\ and significance is assessed against
  the raw p-value (`PValue` column).

- lfc_threshold:

  Minimum absolute log2 fold change required to classify a gene as
  differentially expressed (also sets the vertical threshold lines).
  Default is `1`.

- fdr_threshold:

  P-value cutoff used to draw the horizontal threshold line and to
  classify genes as DE. Applied to `FDR` when `fdr = TRUE` and to
  `PValue` when `fdr = FALSE`. Default is `0.05`.

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

- genes.annot:

  Character vector of gene names to label on the plot via
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Names must match values in the column specified by `gene_col`; an
  error is raised if any names are not found. Default is `NULL` (no
  labels).

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

## Details

The DE classification used for point colouring is computed **inside this
function** from `logFC`, the chosen p-value column (`FDR` or `PValue`),
and the supplied thresholds — so the plot always reflects the thresholds
you pass here, regardless of what was used in
[`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md).

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

# Default plot — y-axis uses FDR-adjusted p-values
ptrap_volcano(res)

# Use raw (unadjusted) p-values on the y-axis
ptrap_volcano(res, fdr = FALSE)

# Custom thresholds and colours
ptrap_volcano(res,
              lfc_threshold = 0.5,
              fdr_threshold = 0.1,
              colors        = c("UP" = "#E69F00", "DOWN" = "#56B4E9"))
} # }
```
