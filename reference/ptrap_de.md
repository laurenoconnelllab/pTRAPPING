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
  `normLibSizes()`, which adjusts effective library sizes (not the count
  matrix itself) for use in the GLM as offsets. `norm.method` is not
  applicable and is ignored. Recommended with **4 or more replicates**
  and raw count data.

- `"deseq"` uses the **DESeq2** pipeline with a multi-factor design
  (`~ block + fraction`; fraction goes last so `results()` auto-tests
  it). DESeq2 applies the **median of ratios normalization** internally
  via `DESeq()`. `norm.method` is not applicable and is ignored.
  Recommended when you prefer the negative-binomial shrinkage estimators
  of DESeq2 over edgeR. **Important:** by default, `logFC` is the
  **maximum likelihood estimate (MLE)** from `results()`, which can be
  noisy and large in magnitude for low-count genes. Set
  `shrink.lfc = TRUE` to apply empirical Bayes LFC shrinkage via
  `lfcShrink()` (requires the `apeglm` package), which is strongly
  recommended before ranking genes or producing volcano plots for
  publication.

- `"voom"` uses the **limma-voom** pipeline with a paired design
  (`~ fraction + block`) and empirical Bayes moderation. Voom internally
  transforms counts to **log2-CPM** (using TMM-adjusted library sizes)
  and derives **precision weights** from the mean-variance trend; these
  two steps are integral to the method. `norm.method` is not applicable
  and is ignored.

- `"paired.ttest"` runs a per-gene **paired t-test between IP and
  INPUT** values across the experimental repeats, following Tan et
  al. (2016) *Cell* 167, 47-59
  [doi:10.1016/j.cell.2016.08.028](https://doi.org/10.1016/j.cell.2016.08.028)
  : *"p value for each gene was calculated as the paired t-test between
  input and immunoprecipitated RPKM values from the three experimental
  repeats."* The count scale for the test is controlled by `norm.method`
  (Tan et al. used `"RPKM"`; all three options are valid; default is
  `"CPM"`). `logFC` is reported as \\\log_2(\bar{\mathrm{IP}} /
  \bar{\mathrm{INPUT}})\\. Recommended for PhosphoTRAP experiments with
  **3 replicates**, where GLM-based dispersion estimates are unreliable.

- `"unpaired.ttest"` compares **fold enrichments (IP/INPUT) between two
  treatment groups** using a Welch (unpaired) t-test on per-animal
  log2(FE) values. Use this to test whether ribosomal association
  differs between conditions. Requires exactly two treatments; use
  `treatment_name` and `control_name` to specify which is which. The
  count scale is controlled by `norm.method` (default is `"CPM"`; all
  three options are valid).

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
  norm.method = c("CPM", "RPKM", "mratios", "none"),
  gene.length = NULL,
  prior.count = 1,
  shrink.lfc = FALSE,
  return_long = FALSE,
  ngenes.out = 20,
  kable.out = FALSE
)
```

## Arguments

- counts_mat:

  A counts matrix in one of two formats:

  - **Matrix** – numeric, genes x samples; column names are sample IDs;
    gene IDs are in `rownames` or supplied via `gene_ids`.

  - **Data frame / tibble** – first column is a character vector of gene
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

  The brain region to subset and analyse (e.g., `"POA"`). When
  `region_col` is `NULL` (default) but `region_name` is supplied, the
  function automatically scans `sample_df` for a non-structural column
  whose values include `region_name`, and uses it as `region_col` if
  exactly one such column is found (a message is shown). Can be `NULL`
  when the data contain only a single region.

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

  Name of the column in `sample_df` containing brain region labels
  (e.g., `"BrainRegion"`). When `NULL` (default) and `region_name` is
  supplied, the column is auto-detected (see `region_name`). Set
  explicitly when the auto-detection is ambiguous or when you prefer to
  be explicit. Leave `NULL` when the data contain only one brain region.

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

  - `"LRT"` – edgeR likelihood ratio test via `glmFit` + `glmLRT`
    (default). TMM normalisation is applied internally; `norm.method` is
    not applicable. Suitable for \>= 4 replicates and raw count data.

  - `"QLF"` – edgeR quasi-likelihood F-test via `glmQLFit` +
    `glmQLFTest`. TMM normalisation is applied internally; `norm.method`
    is not applicable. More conservative than `"LRT"`.

  - `"deseq"` – DESeq2 pipeline with a multi-factor design
    (`~ block + fraction`). Median of ratios normalisation is applied
    internally by `DESeq()`; `norm.method` is not applicable. Results
    obtained via `results(dds, name = "fraction_IP_vs_INPUT")`. `logFC`
    is MLE by default; use `shrink.lfc = TRUE` for shrunken estimates.

  - `"voom"` – limma-voom pipeline with a paired design
    (`~ fraction + block`) and empirical Bayes moderation via `lmFit` +
    `eBayes` + `topTable`. Log2-CPM transformation and voom precision
    weights are applied internally; `norm.method` is not applicable.

  - `"paired.ttest"` – per-gene paired t-test between IP and INPUT
    values across replicates, following Tan et al. (2016). The count
    scale is controlled by `norm.method` (default `"CPM"`; Tan et al.
    used `"RPKM"`). Best suited for n = 3 replicates. P-values adjusted
    with BH.

  - `"unpaired.ttest"` – per-gene Welch unpaired t-test comparing
    log2(FE) values between two treatment groups. The count scale is
    controlled by `norm.method` (default `"CPM"`). Requires exactly two
    treatments; use `treatment_name` and `control_name` to specify which
    group is which.

- norm.method:

  Count normalisation method used by the t-test branches
  (`"paired.ttest"` and `"unpaired.ttest"`). Ignored for `"LRT"`,
  `"QLF"`, `"voom"`, and `"deseq"`, which each handle normalisation
  internally (a message is shown if you explicitly supply `norm.method`
  for those methods). Default is `"CPM"`. One of:

  - `"CPM"` (default) – counts per million, computed from edgeR's
    TMM-adjusted effective library sizes (`lib.size x norm.factors`) via
    [`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html). Because
    TMM adjusts *effective library sizes* rather than the count matrix
    directly, CPM values here reflect both sequencing depth and TMM
    normalisation. Suitable for comparisons between replicates of the
    same sample group.

  - `"RPKM"` – reads per kilobase per million: CPM further divided by
    gene length in kilobases, via
    [`edgeR::rpkm()`](https://rdrr.io/pkg/edgeR/man/cpm.html). Requires
    `gene.length`. Suitable for within-sample comparisons between genes
    (the scale used by Tan et al. 2016 for PhosphoTRAP). **Not**
    recommended for between-sample comparisons.

  - `"mratios"` – DESeq2 **median of ratios** normalisation, computed
    via `estimateSizeFactors()` and retrieved with
    `counts(dds, normalized = TRUE)`. Unlike CPM, this method directly
    rescales the count matrix by sample-specific size factors, making it
    suitable for between-sample comparisons. Use this when you want
    DESeq2-style normalisation for the t-test without running the full
    DESeq2 pipeline.

  - `"none"` – no normalisation is applied; the count matrix is used as
    supplied. Intended for pre-normalised matrices (e.g., TPM, FPKM, or
    any user-normalised values) where an additional normalisation step
    would be redundant or distorting.

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

- shrink.lfc:

  Logical. Only used when `test_method = "deseq"`. If `TRUE`, log2 fold
  changes are shrunk using empirical Bayes estimation via
  `DESeq2::lfcShrink(type = "apeglm")`, which requires the
  [apeglm](https://bioconductor.org/packages/apeglm/) package. Shrunken
  LFCs are more reliable for ranking genes and for visualisation because
  they pull noisy, high-variance estimates (typically from low-count
  genes) towards zero. P-values and FDR are not affected – they always
  come from `results()`. Default is `FALSE` (MLE fold changes are
  returned), but `TRUE` is strongly recommended for any result used in a
  figure or table.

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

- **`"deseq"`**: `Gene`, `baseMean`, `logFC` (MLE log2FoldChange by
  default; empirical Bayes shrunken estimate when `shrink.lfc = TRUE`),
  `lfcSE`, `stat` (Wald statistic; present with MLE, absent when
  `shrink.lfc = TRUE` as apeglm replaces it with a posterior estimate),
  `PValue`, `FDR` (Benjamini-Hochberg from `results()`; unaffected by
  shrinkage), treatment label, optional region label, `diffexpressed`.
  Genes with `NA` p-values (low-count outliers flagged by DESeq2) are
  retained and sorted to the bottom.

- **`"paired.ttest"`**: always returns a **named list** with two
  components:

  - `$results` — tibble with `Gene`, `logFC` (\\\log_2(\bar{\mathrm{IP}}
    / \bar{\mathrm{INPUT}})\\), `t_statistic`, `PValue`, `FDR`,
    treatment label, optional region label, `diffexpressed`.

  - `$fe` — wide tibble with `Gene` plus one `FE_<block>` column per
    animal, containing the per-animal linear IP/INPUT fold enrichment
    values used in the test.

  - `$long_data` — (only when `return_long = TRUE`) a per-gene,
    per-animal tibble with columns `Gene`, block column, `ip_count`,
    `input_count`, `FE`.

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
## ---- Option A: simplest call -- auto-parse from column names ---------------
# counts_mat is a data frame where:
#   col 1      = gene IDs
#   col 2+     = samples named like "b1input", "b1ip", "nb2input", etc.
counts <- read.table("counts.txt", header = TRUE)

# single treatment in the matrix -- treatment_name auto-detected
res <- ptrap_de(counts_mat = counts, test_method = "paired.ttest")

# multiple treatments -- specify which one to analyze
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

# paired t-test -- also return long-format paired table
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

# unpaired t-test -- compare fold enrichment between two treatments
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
