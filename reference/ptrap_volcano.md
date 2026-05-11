# Single volcano plot for PhosphoTRAP / TRAP-seq differential expression results

A volcano plot puts log2 fold change (logFC, x-axis) against statistical
significance (-log p-value, y-axis). Genes that are both strongly
enriched *and* statistically reliable appear in the upper corners — the
most biologically meaningful candidates. Genes that fail either
threshold are shown in grey; those that pass both are coloured by
direction (`"UP"` for enriched in IP, `"DOWN"` for depleted).

## Usage

``` r
ptrap_volcano(
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
  title = NULL
)
```

## Arguments

- de_result:

  A tibble returned by
  [`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md).
  For `test_method = "paired.ttest"`, pass the `$results` component.

- fdr:

  Logical. If `TRUE` (default), significance is assessed using `FDR`
  (false discovery rate — multiple-testing adjusted p-value); if
  `FALSE`, uses raw `PValue`.

- lfc_threshold:

  Minimum absolute log2 fold change for a gene to be coloured as DE. A
  value of `1` means at least a 2× change. Default `1`.

- fdr_threshold:

  Significance cutoff applied to the column selected by `fdr`. Default
  `0.05`.

- log_base:

  Numeric. Base of the logarithm for the p-value axis. Base 10 (default)
  is conventional; base 2 matches Tan et al. (2016). Must be a positive
  number other than `1`.

- gene_col:

  Column name for gene identifiers. Default `"Gene"`.

- treatment_col:

  Column name for treatment label. Default `"treatment"`.

- region_col:

  Column name for brain region label. Default `"BrainRegion"`.

- colors:

  Named vector mapping `"UP"` and `"DOWN"` to colours. Default is
  colourblind-friendly: `"#0072b2"` (UP), `"#E69f00"` (DOWN).

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

Takes a single tibble from
[`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md)
(for `test_method = "paired.ttest"`, pass the `$results` component).
Gene labels are drawn only for genes listed in `genes.annot`, via
[`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
The DE classification is recomputed inside this function from the
supplied thresholds, so you can explore different cutoffs without
re-running
[`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: -log10(FDR)
ptrap_volcano(res)

# Raw p-values, log2 scale -- as in Tan et al. 2016
ptrap_volcano(res, fdr = FALSE, log_base = 2)
} # }
```
