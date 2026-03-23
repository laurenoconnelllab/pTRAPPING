# Perform differential expression analysis for TRAP-seq data using edgeR or a paired t-test

Compares IP vs INPUT fractions for a specified treatment condition using
one of three statistical approaches (see `test_method`). Both
`sample_df` and `gene_ids` are optional: the function can derive sample
metadata automatically from the column names of `counts_mat`, and gene
identifiers from its first column.

## Usage

``` r
ptrap_de(
  counts_mat,
  sample_df = NULL,
  gene_ids = NULL,
  region_name = NULL,
  treatment_name = NULL,
  sample_col = "sample",
  fraction_col = "fraction",
  block_col = "tube",
  region_col = NULL,
  treatment_col = "Treatment",
  ip_level = "IP",
  input_level = "INPUT",
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  test_method = c("LRT", "QLF", "paired.ttest"),
  pseudocount = 1,
  return_long = FALSE,
  ngenes.out = 20,
  kable.out = FALSE
)
```

## Arguments

- counts_mat:

  A counts matrix in one of two formats:

  - **Matrix** — numeric, genes × samples; column names are sample IDs;
    gene IDs are in `rownames` or supplied via `gene_ids`.

  - **Data frame / tibble** — first column is a character vector of gene
    IDs; remaining columns are numeric counts with sample names as
    column names. When `sample_df = NULL`, column names must follow the
    naming convention described in the *Automatic column-name parsing*
    section.

- sample_df:

  Optional data frame of sample metadata. When `NULL` (default),
  metadata is parsed automatically from the column names of
  `counts_mat`. When provided, the arguments `sample_col`,
  `fraction_col`, `block_col`, `region_col`, and `treatment_col` specify
  which columns to use (falling back to their defaults if the names
  match).

- gene_ids:

  Optional character vector of gene identifiers corresponding to the
  rows of `counts_mat`. When `NULL` (default), gene IDs are extracted
  from the first column of `counts_mat` (if it is a data frame) or from
  `rownames(counts_mat)` (if it is a matrix).

- region_name:

  The brain region to subset and analyze (e.g., `"POA"`). Only used when
  `region_col` is not `NULL`. Can be `NULL` when `region_col = NULL`
  (default) or when the data contain a single region.

- treatment_name:

  The treatment condition whose samples will be **subsetted** for the IP
  vs INPUT comparison (e.g., `"pb"` to analyse only the pair-bonded
  samples). Note that `treatment_name` is **not** the DE contrast
  variable — the contrast is always IP vs INPUT. To compare two
  treatments you call `ptrap_de()` once per treatment and then pass both
  results to
  [`ptrap_volcano2()`](https://laurenoconnelllab.github.io/pTRAPPING/reference/ptrap_volcano2.md).

  When `NULL` (default), the value is derived automatically from
  `sample_df` (or from the parsed column names when `sample_df = NULL`):
  if a single treatment is found the function proceeds silently; if
  multiple treatments are found an error asks the user to specify which
  one to analyse.

- sample_col:

  Name of the column in `sample_df` whose values match the column names
  of `counts_mat`. Ignored when `sample_df = NULL` (auto-set
  internally). Default is `"sample"`.

- fraction_col:

  Name of the column in `sample_df` that distinguishes IP from INPUT
  fractions. Default is `"fraction"`.

- block_col:

  Name of the column in `sample_df` used as the blocking / pairing
  variable (e.g., individual animal or tube). For `"LRT"` / `"QLF"` it
  enters the design matrix; for `"paired.ttest"` it aligns each animal's
  IP with its own INPUT. Default is `"tube"`.

- region_col:

  Name of the column in `sample_df` containing brain region labels. Set
  to `NULL` (default) to skip region filtering — recommended when the
  data already contain only one brain region. Only set this when
  `counts_mat` contains multiple regions and you want to analyze one at
  a time. Default is `NULL`.

- treatment_col:

  Name of the column in `sample_df` containing treatment labels. Default
  is `"Treatment"`.

- ip_level:

  The value in `fraction_col` that identifies the IP fraction. Default
  is `"IP"`.

- input_level:

  The value in `fraction_col` that identifies the INPUT fraction
  (reference level). Default is `"INPUT"`.

- lfc_threshold:

  Minimum absolute log2 fold change required to classify a gene as
  differentially expressed. Default is `1`.

- fdr_threshold:

  Maximum FDR (or adjusted p-value for `"paired.ttest"`) allowed to
  classify a gene as differentially expressed. Default is `0.05`.

- test_method:

  Statistical method to use. One of:

  - `"LRT"` — likelihood ratio test via `glmFit` + `glmLRT` (default).
    Suitable for ≥ 4 replicates and raw count data.

  - `"QLF"` — quasi-likelihood F-test via `glmQLFit` + `glmQLFTest`.
    More conservative than LRT; recommended for small sample sizes with
    raw count data.

  - `"paired.ttest"` — per-gene paired t-test between IP and INPUT
    across replicates, following Tan et al. (2016)
    [doi:10.1016/j.cell.2016.08.028](https://doi.org/10.1016/j.cell.2016.08.028)
    . Best suited for n = 3 replicates typical of PhosphoTRAP
    experiments. P-values are adjusted with Benjamini-Hochberg. The
    enrichment ratio per paired sample (IP/INPUT) is expressed on a log2
    scale.

- pseudocount:

  A small value added to counts before computing log2(IP/INPUT) in
  `"paired.ttest"`, to avoid log(0). Default is `1`. Set to `0` if
  `counts_mat` already contains normalised values (e.g., RPKM/TPM) that
  are guaranteed to be \> 0.

- return_long:

  Logical. Only used when `test_method = "paired.ttest"`. If `TRUE`,
  returns a named list with `$results` (the DE tibble) and `$long_data`
  (the per-gene, per-animal paired table used for the test). Default is
  `FALSE`.

- ngenes.out:

  Number of top genes (sorted by p-value) to include in the output when
  `kable.out = TRUE`. Default is `20`.

- kable.out:

  Logical. If `TRUE`, returns a `kableExtra` HTML table of the top
  `ngenes.out` genes instead of the full tibble. Default is `FALSE`.

## Value

When `kable.out = FALSE` (default), a tibble with one row per gene
sorted by p-value, containing `Gene`, `logFC`, test-statistic column(s),
`PValue`, `FDR`, treatment label, optionally a region label (when
`region_col` is set), and `diffexpressed` (`"UP"`, `"DOWN"`, `"NO"`).
When `kable.out = TRUE`, an HTML `kableExtra` table of the top
`ngenes.out` genes. For `"paired.ttest"` with `return_long = TRUE`, a
named list with `$results` and `$long_data`.

## Details

- `"LRT"` and `"QLF"` use edgeR's GLM framework, account for the paired
  animal structure via a block variable, and are recommended when you
  have **4 or more replicates** and raw count data.

- `"paired.ttest"` runs a per-gene **paired t-test** between IP and
  INPUT counts across experimental repeats, following the method
  described in Tan et al. (2016) *Cell* 167, 47–59
  [doi:10.1016/j.cell.2016.08.028](https://doi.org/10.1016/j.cell.2016.08.028)
  . This approach is recommended when you have **only 3 replicates** (as
  is the typical minimum in PhosphoTRAP experiments), where GLM-based
  dispersion estimates are unreliable. The enrichment ratio (IP / INPUT)
  is computed per paired sample and expressed on the log2 scale, as
  described in the same paper. **Note:** if you pass normalised values
  (e.g., RPKM/TPM) instead of raw counts, set `pseudocount = 0`.

## Automatic column-name parsing

When `sample_df = NULL`, the column names of `counts_mat` (excluding the
first, gene-ID column) are parsed to build sample metadata
automatically. Each column name must encode three pieces of information,
in any order and with any combination of separators (`_`, `-`, `.`,
space) or no separator at all:

- Treatment:

  One or more letters identifying the experimental group.

- Replicate number:

  A digit identifying the biological replicate.

- Fraction:

  `ip_level` or `input_level` (case-insensitive); `"in"` is also
  accepted as a short alias for the INPUT fraction.

Valid column name examples (default `ip_level = "IP"`,
`input_level = "INPUT"`):

|                   |                                  |
|-------------------|----------------------------------|
| **Column name**   | **Parsed as**                    |
| `b1input`, `b1ip` | treatment = `b`, block = `1`     |
| `nb2INPUT`        | treatment = `nb`, block = `2`    |
| `Nb_IP_1`         | treatment = `Nb`, block = `1`    |
| `B_3_INPUT`       | treatment = `B`, block = `3`     |
| `PB.2.ip`         | treatment = `PB`, block = `2`    |
| `Trim_10-INPUT`   | treatment = `Trim`, block = `10` |
| `SOL1INPUT`       | treatment = `SOL`, block = `1`   |

## References

Tan, C.L., Cooke, E.K., Leib, D.E., Lin, Y.C., Daly, G.E., Zimmerman,
C.A., and Knight, Z.A. (2016). Warm-Sensitive Neurons that Control Body
Temperature. *Cell* 167, 47–59.
[doi:10.1016/j.cell.2016.08.028](https://doi.org/10.1016/j.cell.2016.08.028)

## Examples

``` r
if (FALSE) { # \dontrun{
## ---- Option A: simplest call — auto-parse from column names ---------------
# counts_mat is a data frame where:
#   col 1      = gene IDs
#   col 2+     = samples named like "b1input", "b1ip", "nb2input", etc.
counts <- read.table("counts.txt", header = TRUE)

# single treatment in the matrix — treatment_name auto-detected
res <- ptrap_de(counts_mat = counts, test_method = "paired.ttest")

# multiple treatments — specify which one to analyze
res_b <- ptrap_de(counts_mat = counts, treatment_name = "b")

## ---- Option B: provide sample_df explicitly (original workflow) -----------
res_lrt <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  region_name    = "POA",
  treatment_name = "pb"
)

# quasi-likelihood F-test
res_qlf <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  region_name    = "POA",
  treatment_name = "pb",
  test_method    = "QLF"
)

# paired t-test — also return long-format paired table
res_pt <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  region_name    = "POA",
  treatment_name = "pb",
  test_method    = "paired.ttest",
  return_long    = TRUE
)
res_pt$results    # DE tibble
res_pt$long_data  # per-gene, per-animal pairs
} # }
```
