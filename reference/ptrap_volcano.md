# Single volcano plot for TRAP-seq or RNA-seq differential expression results

Takes a single tibble produced by
[`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md)
and returns a classic volcano plot (logFC on the x-axis, -log_base(p) on
the y-axis). Non-significant genes are shown in grey; genes classified
as `"UP"` or `"DOWN"` are highlighted with distinct fill colours. Gene
labels are added via
[`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
**only for genes supplied in `genes.annot`**.

## Usage

``` r
ptrap_volcano(
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
)
```

## Arguments

- de_result:

  A tibble returned by
  [`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md).
  For `test_method = "paired.ttest"`, pass the `$results` component.

- fdr:

  Logical. If `TRUE` (default), uses `FDR`; if `FALSE`, uses `PValue`.

- lfc_threshold:

  Minimum absolute log2 fold change. Default `1`.

- fdr_threshold:

  P-value cutoff. Default `0.05`.

- log_base:

  Numeric. Base of the logarithm for the p-value axis. Default is `10`
  (-log10, current behaviour). Use `2` for -log2(p) as in Tan et al.
  (2016). Must be a positive number other than `1`.

- gene_col:

  Column name for gene identifiers. Default `"Gene"`.

- treatment_col:

  Column name for treatment label. Default `"Treatment"`.

- region_col:

  Column name for brain region label. Default `"BrainRegion"`.

- colors:

  Named vector mapping `"UP"` and `"DOWN"` to colours. Default is
  colourblind-friendly: `"#D55E00"` (UP), `"#0072B2"` (DOWN).

- point_size:

  Size of points. Default `3.5`.

- point_alpha:

  Opacity of highlighted points. Default `0.7`.

- genes.annot:

  Character vector of gene names to label. Default `NULL`.

- max_overlaps:

  Passed to
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Default `20`.

- title:

  Plot title. Auto-generated from region and treatment if `NULL`.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

The DE classification is recomputed inside this function from `logFC`,
the chosen p-value column, and the supplied thresholds.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: -log10(FDR)
ptrap_volcano(res)

# Raw p-values, log2 scale -- as in Tan et al. 2016
ptrap_volcano(res, fdr = FALSE, log_base = 2)
} # }
```
