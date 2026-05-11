# pTRAPPING

**pTRAPPING** takes the raw count data from a TRAP-seq or PhosphoTRAP
experiment and finds which genes are specifically enriched in your cell
population of interest — all in a few lines of R. You give it a counts
table where every sample column is labelled with its treatment group,
replicate number, and whether it is IP (immunoprecipitated) or INPUT
(total RNA); the package figures out the rest. Three functions cover the
full workflow: differential expression with your choice of statistical
engine
([`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md)),
a classic volcano plot for one condition
([`ptrap_volcano()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_volcano.md)),
and a paired scatter plot to compare two conditions head-to-head
([`ptrap_volcano2()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_volcano2.md)).

## Installation

``` r

# install.packages("devtools")
devtools::install_github("laurenoconnelllab/pTRAPPING")
```

Dependencies from Bioconductor (edgeR, limma, DESeq2) are installed
automatically.

## Functions at a glance

| Function | What it does |
|----|----|
| [`ptrap_de()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_de.md) | Differential expression: IP vs INPUT. Six statistical methods. |
| [`ptrap_volcano()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_volcano.md) | Volcano plot for a single condition. |
| [`ptrap_volcano2()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_volcano2.md) | Scatter plot comparing two conditions side-by-side. |

## Quick start

``` r

library(pTRAPPING)

# Load counts matrix
counts.mat <- read.delim(
  system.file("extdata", "TAN_etal_2016_raw.txt", package = "pTRAPPING")
)

# Table of differentially expressed genes of interest for one condition using the default method "LRT" from edgeR
ptrap_de(
  counts_mat = counts.mat,
  treatment_name = "PACAP",
  kable.out = TRUE,
  genes.filter = c(
    "Adcyap1",
    "Bdnf",
    "Ucn3",
    "Gng8",
    "Fosl2",
    "Junb",
    "Trappc12",
    "Gfap"
  )
)
```

| Gene     |  logFC | logCPM |     LR | PValue  | FDR     | treatment | diffexpressed |
|:---------|-------:|-------:|-------:|:--------|:--------|:----------|:--------------|
| Gng8     |  4.056 |  5.578 | 41.232 | \<0.001 | \<0.001 | PACAP     | UP            |
| Ucn3     |  3.948 |  5.218 | 25.888 | \<0.001 | \<0.001 | PACAP     | UP            |
| Bdnf     |  4.336 |  3.090 | 23.351 | \<0.001 | \<0.001 | PACAP     | UP            |
| Gfap     | -2.979 |  5.086 | 10.697 | 0.001   | 0.028   | PACAP     | DOWN          |
| Trappc12 | -3.350 |  1.809 | 10.479 | 0.001   | 0.030   | PACAP     | DOWN          |
| Fosl2    |  0.781 |  6.990 |  3.988 | 0.046   | 0.266   | PACAP     | NO            |
| Adcyap1  |  1.516 |  3.665 |  2.691 | 0.101   | 0.387   | PACAP     | NO            |
| Junb     |  0.757 |  4.568 |  1.343 | 0.247   | 0.597   | PACAP     | NO            |

``` r

# Volcano plot — label a few genes of interest
ptrap_de(
  counts_mat = counts.mat,
  treatment_name = "PACAP"
) |>
  ptrap_volcano(
    genes.annot = c(
      "Adcyap1",
      "Bdnf",
      "Ucn3",
      "Gng8",
      "Fosl2",
      "Junb",
      "Trappc12",
      "Gfap"
    ),
    title = "PACAP"
  )
```

![Volcano plot of IP vs INPUT enrichment in
PACAP](reference/figures/README-quick-volcano-1.png)

## Learn more

- **Full walkthrough:**
  [`vignette("getting-started", package = "pTRAPPING")`](https://laurenoconnelllab.github.io/pTRAPPING/articles/getting-started.md)

- The example dataset is from:

> Tan, C.L., Cooke, E.K., Leib, D.E., Lin, Y.-C., Daly, G.E., Zimmerman,
> C.A., and Knight, Z.A. (2016). Warm-sensitive neurons that control
> body temperature. *Cell* 167, 47–59.
> <https://doi.org/10.1016/j.cell.2016.08.028>

- The PhosphoTRAP method was introduced in:

> Knight, Z.A., Tan, K., Birsoy, K., Schmidt, S., Garrison, J.L.,
> Wysocki, R.W., Emiliano, A., Ekstrand, M.I., and Friedman, J.M.
> (2012). Molecular profiling of activated neurons by phosphorylated
> ribosome capture. *Cell* 151, 1126–1137.
> <https://doi.org/10.1016/j.cell.2012.10.039>
