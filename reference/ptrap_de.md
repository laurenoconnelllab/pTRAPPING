# Differential expression analysis for PhosphoTRAP data

Compares IP vs INPUT fractions for a specified treatment condition using
one of six statistical approaches (see `test_method`). Both `sample_df`
and `gene_ids` are optional: the function can derive sample metadata
automatically from the column names of `counts_mat`, and gene
identifiers from its first column.

- `"LRT"` and `"QLF"` use edgeR's GLM framework with a multi-factor
  design (`~ fraction + block`), accounting for the paired animal
  structure via the blocking variable. Both methods automatically apply
  **edgeR's TMM (trimmed mean of M values) normalization** via
  `normLibSizes()`, which normalises effective library sizes (not the
  counts directly) for use in the GLM. `norm.method` is automatically
  forced to `"none"` for these methods — CPM and RPKM values, when
  needed, are computed from the TMM-adjusted effective library sizes.
  Recommended with **4 or more replicates** and raw count data.

- `"deseq"` uses the **DESeq2** pipeline with a multi-factor design
  (`~ block + fraction`; fraction is last so `results()` auto-tests it).
  DESeq2 applies its own **median of ratios normalization** internally
  via `DESeq()`, so `norm.method` is automatically forced to `"none"`.
  Recommended when you prefer the negative-binomial shrinkage estimators
  of DESeq2 over edgeR.

- `"voom"` uses the **limma-voom** pipeline with a paired design
  (`~ fraction + block`) and empirical Bayes moderation. Voom internally
  transforms counts to **log2-CPM** (using TMM-adjusted library sizes)
  and computes **precision weights** from the mean-variance trend.
  Because this transformation is integral to voom, `norm.method` is
  automatically forced to `"none"`.

- `"paired.ttest"` runs a per-gene **paired t-test between IP and
  INPUT** values across the experimental repeats, following Tan et
  al. (2016) *Cell* 167, 47-59
  [doi:10.1016/j.cell.2016.08.028](https://doi.org/10.1016/j.cell.2016.08.028)
  : *"p value for each gene was calculated as the paired t-test between
  input and immunoprecipitated RPKM values from the three experimental
  repeats."* The count scale for the test is set by `norm.method` (Tan
  et al. used `"RPKM"`; all four `norm.method` options are valid).
  `logFC` is reported as \\\log_2(\bar{\mathrm{IP}} /
  \bar{\mathrm{INPUT}})\\. Recommended for PhosphoTRAP experiments with
  **3 replicates**, where GLM-based dispersion estimates are unreliable.

- `"unpaired.ttest"` compares **fold enrichments (IP/INPUT) between two
  treatment groups** using a Welch (unpaired) t-test on per-animal
  log2(FE) values. Use this to test whether ribosomal association
  differs between conditions. Requires exactly two treatments; use
  `treatment_name` and `control_name` to specify which is which. All
  four `norm.method` options are valid.

## Usage

``` r
ptrap_de(
  counts_mat,
  sample_df = NULL,
  gene_ids = NULL,
  region_name = NULL,
  treatment_name = NULL,
  control_name = NULL,
  sample_col = "sample",
  fraction_col = "fraction",
  block_col = "tube",
  region_col = NULL,
  treatment_col = "Treatment",
  ip_level = "IP",
  input_level = "INPUT",
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  test_method = c("LRT", "QLF", "paired.ttest", "unpaired.ttest", "voom", "deseq"),
  norm.method = c("none", "CPM", "RPKM", "mratios"),
  gene.length = NULL,
  prior.count = 1,
  return_long = FALSE,
  ngenes.out = 20,
  kable.out = FALSE
)
```

## Arguments

- counts_mat:

  A counts matrix in one of two formats:

  - **Matrix** — numeric, genes x samples; column names are sample IDs;
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
  samples). For `"unpaired.ttest"`, this is the numerator condition in
  the differential fold enrichment (DFE = mean_FE_treatment /
  mean_FE_control).

  When `NULL` (default), the value is derived automatically from
  `sample_df` (or from the parsed column names when `sample_df = NULL`):
  if a single treatment is found the function proceeds silently; if
  multiple treatments are found an error asks the user to specify which
  one to analyse. For `"unpaired.ttest"` with exactly two treatments and
  both `treatment_name` and `control_name` as `NULL`, assignments are
  made alphabetically (first = control, second = treatment) with a
  message.

- control_name:

  The **reference / control** treatment used only by `"unpaired.ttest"`.
  Serves as the denominator when computing differential fold enrichment
  (DFE = mean_FE_treatment / mean_FE_control). When `NULL` (default),
  the function auto-detects the control as the non-`treatment_name`
  group when exactly two treatments are present.

- sample_col:

  Name of the column in `sample_df` whose values match the column names
  of `counts_mat`. Ignored when `sample_df = NULL` (auto-set
  internally). Default is `"sample"`.

- fraction_col:

  Name of the column in `sample_df` that distinguishes IP from INPUT
  fractions. Default is `"fraction"`.

- block_col:

  Name of the column in `sample_df` used as the blocking / pairing
  variable (e.g., individual animal or tube). For `"LRT"` / `"QLF"`,
  `"deseq"`, and `"voom"` it enters the design matrix; for
  `"paired.ttest"` and `"unpaired.ttest"` it aligns each animal's IP
  with its own INPUT. Default is `"tube"`.

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

  Maximum FDR (or adjusted p-value) allowed to classify a gene as
  differentially expressed. Default is `0.05`.

- test_method:

  Statistical method to use. One of:

  - `"LRT"` — edgeR likelihood ratio test via `glmFit` + `glmLRT`
    (default). TMM normalization applied automatically; `norm.method` is
    forced to `"none"`. Suitable for \>= 4 replicates and raw count
    data.

  - `"QLF"` — edgeR quasi-likelihood F-test via `glmQLFit` +
    `glmQLFTest`. TMM normalization applied automatically; `norm.method`
    is forced to `"none"`. More conservative than LRT.

  - `"deseq"` — DESeq2 pipeline with a multi-factor design
    (`~ block + fraction`). Median of ratios normalization is applied
    automatically by `DESeq()`; `norm.method` is forced to `"none"`.
    Results are obtained via
    `results(dds, name = "fraction_IP_vs_INPUT")`.

  - `"voom"` — limma-voom pipeline with a paired design
    (`~ fraction + block`) and empirical Bayes moderation via `lmFit` +
    `eBayes` + `topTable`. Log2-CPM + voom precision weights applied
    automatically; `norm.method` is forced to `"none"`.

  - `"paired.ttest"` — per-gene paired t-test between IP and INPUT
    values across replicates, following Tan et al. (2016). The count
    scale is set by `norm.method` (Tan et al. used `"RPKM"`; all four
    options are valid). Best suited for n = 3 replicates. P-values
    adjusted with BH.

  - `"unpaired.ttest"` — per-gene Welch unpaired t-test comparing
    log2(FE) values between two treatment groups. The count scale is set
    by `norm.method` (all four options are valid). Requires exactly two
    treatments; use `treatment_name` and `control_name` to specify
    groups.

- norm.method:

  Normalization method applied to counts before the t-test
  (`"paired.ttest"` / `"unpaired.ttest"` only). **Ignored and
  automatically set to `"none"` for `"LRT"`, `"QLF"`, `"deseq"`, and
  `"voom"`**, which each apply their own normalization internally (a
  message is shown if you supply a non-`"none"` value). One of:

  - `"none"` (default) — uses the TMM-adjusted raw counts stored in the
    DGEList (i.e., the count matrix is not rescaled, but effective
    library sizes reflect TMM factors). **Note:** TMM normalises
    *effective library sizes*, not the counts themselves; CPM and RPKM
    values are derived from these TMM-adjusted sizes.

  - `"CPM"` — counts per million, computed from TMM-adjusted effective
    library sizes via
    [`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html).
    Recommended for comparisons between replicates of the same sample
    group.

  - `"RPKM"` — reads per kilobase per million (CPM further divided by
    gene length in kilobases), via
    [`edgeR::rpkm()`](https://rdrr.io/pkg/edgeR/man/cpm.html). Requires
    `gene.length`. Recommended for within-sample gene-to-gene
    comparisons (as used in Tan et al. 2016). **Not** suitable for
    between-sample comparisons.

  - `"mratios"` — DESeq2 **median of ratios** normalization via
    `estimateSizeFactors()` followed by
    `counts(dds, normalized = TRUE)`. Produces normalized count values
    that are suitable for between-sample comparisons. Recommended when
    you want DESeq2-style normalization for the t-test without running
    the full DESeq2 pipeline.

- gene.length:

  Named numeric vector of gene lengths in base pairs (names must match
  gene IDs). Required when `norm.method = "RPKM"`; ignored otherwise.
  Default is `NULL`.

- prior.count:

  A small count added to IP and INPUT values before computing log2
  ratios (FE = (IP + prior.count) / (INPUT + prior.count)), to avoid
  log(0). Default is `1`. Set to `0` if values are already normalised
  and guaranteed positive. Used only by `"paired.ttest"` and
  `"unpaired.ttest"`.

- return_long:

  Logical. Only used when `test_method` is `"paired.ttest"` or
  `"unpaired.ttest"`. If `TRUE`, returns a named list with `$results`
  (the DE tibble) and `$long_data` (the per-gene, per-animal table used
  for the test). Default is `FALSE`.

- ngenes.out:

  Number of top genes (sorted by p-value) to include in the output when
  `kable.out = TRUE`. Default is `20`.

- kable.out:

  Logical. If `TRUE`, returns a `kableExtra` HTML table of the top
  `ngenes.out` genes instead of the full tibble. Default is `FALSE`.

## Value

When `kable.out = FALSE` (default), a tibble with one row per gene
sorted by p-value. Columns depend on `test_method`:

- **`"LRT"` / `"QLF"`**: `Gene`, `logFC`, `logCPM`, `LR` or `F`,
  `PValue`, `FDR`, treatment label, optional region label,
  `diffexpressed` (`"UP"`, `"DOWN"`, `"NO"`).

- **`"deseq"`**: `Gene`, `baseMean`, `logFC` (log2FoldChange from
  DESeq2), `lfcSE`, `stat`, `PValue`, `FDR` (Benjamini-Hochberg adjusted
  p-value from `results()`), treatment label, optional region label,
  `diffexpressed`. Genes with `NA` p-values (outliers or low counts
  flagged by DESeq2) are retained and sorted to the bottom.

- **`"paired.ttest"`**: `Gene`, `logFC` (\\\log_2(\bar{\mathrm{IP}} /
  \bar{\mathrm{INPUT}})\\), `t_statistic` (paired t-test on IP \\-\\
  INPUT differences), `PValue`, `FE_<block>` columns (one per animal
  giving the linear IP/INPUT fold enrichment for that animal), `FDR`,
  treatment label, optional region label, `diffexpressed`. With
  `return_long = TRUE`, returns a named list: `$results` (the tibble
  above) and `$long_data` (a per-gene, per-animal tibble with columns
  `Gene`, block column, `ip_count`, `input_count`, `FE`).

- **`"unpaired.ttest"`**: `Gene`, `logFC` (log2(mean_FE_treatment /
  mean_FE_control)), `diff_FE`, `mean_FE_<treatment>`,
  `mean_FE_<control>`, `t_statistic`, `df` (Welch degrees of freedom),
  `PValue`, `FDR`, treatment label, optional region label,
  `diffexpressed`. With `return_long = TRUE`, returns a named list:
  `$results` and `$long_data` (per-gene, per-animal, per-treatment
  tibble with columns `Gene`, block column, treatment column, `FE`,
  `log2_FE`).

- **`"voom"`**: `Gene`, `logFC`, `AveExpr`, `t`, `B`, `PValue`, `FDR`,
  treatment label, optional region label, `diffexpressed`.

When `kable.out = TRUE`, an HTML `kableExtra` table of the top
`ngenes.out` genes.

## Details

Perform differential expression analysis for PhosphoTRAP data using
edgeR, a paired t-test, an unpaired t-test, or voom/limma

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
Temperature. *Cell* 167, 47-59.
[doi:10.1016/j.cell.2016.08.028](https://doi.org/10.1016/j.cell.2016.08.028)

Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of
fold change and dispersion for RNA-seq data with DESeq2. *Genome
Biology* 15, 550.
[doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

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

# unpaired t-test — compare fold enrichment between two treatments
res_upt <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  treatment_name = "pb",
  control_name   = "nb",
  test_method    = "unpaired.ttest"
)

# voom/limma pipeline
res_voom <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  region_name    = "POA",
  treatment_name = "pb",
  test_method    = "voom"
)

# paired t-test with CPM normalization
res_cpm <- ptrap_de(
  counts_mat     = counts_mat,
  sample_df      = sample_df,
  gene_ids       = gene_ids,
  treatment_name = "pb",
  test_method    = "paired.ttest",
  norm.method    = "CPM"
)
} # }
```
