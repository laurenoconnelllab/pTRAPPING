# Differential expression analysis for PhosphoTRAP data

In PhosphoTRAP and TRAP-seq experiments, RNA from a target cell
population is physically isolated by immunoprecipitation — the **IP**
fraction. Each sample also provides an **INPUT**: total RNA before the
pulldown, serving as the background reference. `ptrap_de()` asks:
**which genes are enriched (or depleted) in the IP fraction relative to
input?** A positive log2 fold change (logFC \> 0) means a gene's RNA is
more abundant in the IP; a negative logFC means it is less abundant.

Six statistical methods are available via `test_method`. The right
choice depends on your replicate count and experimental question:

- `"LRT"` and `"QLF"` use **edgeR**'s generalized linear model (GLM)
  framework with a design that pairs each animal's IP with its own input
  (`~ fraction + block`). Both apply **TMM normalization** (trimmed mean
  of M values) internally, which corrects for differences in sample
  composition by rescaling effective library sizes rather than the count
  matrix itself. `"QLF"` is more conservative than `"LRT"`.
  `norm.method` is not applicable and is ignored. Recommended with **4
  or more replicates** and raw count data.

- `"deseq"` uses the **DESeq2** pipeline with the same paired design
  (`~ block + fraction`). Normalization (median of ratios) is handled
  internally; `norm.method` is ignored. **Important:** by default,
  `logFC` is the raw maximum-likelihood estimate (MLE), which can appear
  unrealistically large for genes with very few reads. Set
  `shrink.lfc = TRUE` to pull these noisy estimates towards zero —
  strongly recommended before ranking genes or producing publication
  figures.

- `"voom"` uses the **limma-voom** pipeline, which bridges count data
  and linear models: it converts counts to log2-CPM and assigns each
  observation a precision weight based on its mean-variance
  relationship, then fits a linear model with empirical Bayes
  moderation. `norm.method` is not applicable and is ignored.

- `"paired.ttest"` runs a per-gene **paired t-test** between IP and
  INPUT values, treating each animal like a matched before/after
  measurement. This is the approach from Tan et al. (2016) *Cell* 167,
  47–59
  [doi:10.1016/j.cell.2016.08.028](https://doi.org/10.1016/j.cell.2016.08.028)
  , who used `"RPKM"` normalization; the default here is `"CPM"`.
  Recommended for **3-replicate** experiments, where GLM-based methods
  struggle to estimate variability reliably. `logFC` is
  \\\log_2(\bar{\mathrm{IP}} / \bar{\mathrm{INPUT}})\\.

- `"unpaired.ttest"` answers a *different question* than the other
  methods: not "is IP enriched over input?" but "does the *degree* of
  enrichment differ between two treatment groups (e.g., pair-bonded vs.
  non-bonded animals)?" It runs a per-gene Welch t-test on per-animal
  log2(IP/INPUT) fold enrichment values. Requires exactly two
  treatments; use `treatment_name` and `control_name` to specify which
  is which. Count scale is controlled by `norm.method`.

Both `sample_df` and `gene_ids` are optional: sample metadata can be
parsed automatically from column names of `counts_mat`, and gene IDs
from its first column or row names.

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
  covariates = NULL,
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  test_method = c("LRT", "QLF", "paired.ttest", "unpaired.ttest", "voom", "deseq"),
  norm.method = c("CPM", "RPKM", "mratios", "none"),
  gene.length = NULL,
  prior.count = 1,
  shrink.lfc = FALSE,
  return_long = FALSE,
  ngenes.out = 20,
  genes.filter = NULL,
  kable.out = FALSE,
  filter = TRUE
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

  **Auto-parsing mode** (`sample_df = NULL`): `region_name` is
  additionally used during column-name parsing to extract the region
  token directly from each column name (e.g., `"POA"` in `b1inputPOA`).
  See the *Automatic column-name parsing* section for details.

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

  Name of the column in `sample_df` that identifies which IP and INPUT
  samples come from the **same animal** (the blocking / pairing
  variable, e.g., tube ID or animal number). Linking each IP to its own
  input control — like a matched before/after comparison within a
  subject — removes animal-to-animal variability from the test. For
  GLM-based methods (`"LRT"`, `"QLF"`, `"deseq"`, `"voom"`) it enters
  the design formula; for t-test methods it aligns paired samples by
  position. Default is `"tube"`.

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

- covariates:

  Optional character vector of covariate names to include in the design
  formula. Default is `NULL`.

- lfc_threshold:

  Minimum absolute log2 fold change required to classify a gene as
  differentially expressed. A threshold of `1` means a gene must be at
  least 2× more (or less) abundant in IP vs. input. Default is `1`.

- fdr_threshold:

  Maximum FDR (false discovery rate — the expected proportion of
  significant findings that are false alarms, adjusted for multiple
  testing) allowed to classify a gene as differentially expressed.
  Default is `0.05`.

- test_method:

  Statistical method to use. See `@description` for a full explanation
  of each. One of:

  - `"LRT"` (default) — edgeR likelihood ratio test (`glmFit` +
    `glmLRT`). Recommended for \>= 4 replicates and raw count data. TMM
    normalisation applied internally; `norm.method` ignored.

  - `"QLF"` — edgeR quasi-likelihood F-test (`glmQLFit` + `glmQLFTest`).
    More conservative than `"LRT"`. TMM normalisation applied
    internally; `norm.method` ignored.

  - `"deseq"` — DESeq2 pipeline (`~ block + fraction`). Median-of-ratios
    normalisation applied internally; `norm.method` ignored. `logFC` is
    MLE by default; use `shrink.lfc = TRUE` for publication-ready
    estimates.

  - `"voom"` — limma-voom pipeline (`~ fraction + block`) with empirical
    Bayes moderation (`lmFit` + `eBayes` + `topTable`). Log2-CPM
    transformation and precision weights applied internally;
    `norm.method` ignored.

  - `"paired.ttest"` — per-gene paired t-test (IP vs. INPUT across
    replicates). Best for n = 3 replicates. Count scale set by
    `norm.method` (default `"CPM"`; Tan et al. 2016 used `"RPKM"`).
    P-values adjusted with Benjamini-Hochberg (BH).

  - `"unpaired.ttest"` — per-gene Welch t-test comparing log2(IP/INPUT)
    between two treatment groups. Count scale set by `norm.method`.
    Requires exactly two treatments (set via `treatment_name` and
    `control_name`).

- norm.method:

  Count normalisation method used only by the t-test branches
  (`"paired.ttest"` and `"unpaired.ttest"`). Ignored for `"LRT"`,
  `"QLF"`, `"voom"`, and `"deseq"`, which each handle normalisation
  internally (a message is shown if you set this for those methods). The
  goal of normalisation is to make counts comparable across samples by
  accounting for differences in sequencing depth and library
  composition. Default is `"CPM"`. One of:

  - `"CPM"` (default) — counts per million, scaled by TMM-adjusted
    library sizes via
    [`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html). Suitable
    for comparing the same gene across replicates of the same group.

  - `"RPKM"` — CPM divided by gene length in kilobases, via
    [`edgeR::rpkm()`](https://rdrr.io/pkg/edgeR/man/cpm.html). Requires
    `gene.length`. This was the scale used by Tan et al. 2016 for
    PhosphoTRAP. Useful for within-sample comparisons between genes of
    different lengths, but **not** recommended for between-sample
    comparisons.

  - `"mratios"` — DESeq2 median-of-ratios normalisation: each sample's
    counts are divided by a size factor estimated from the geometric
    mean of all genes. This directly rescales the count matrix and is
    suitable for between-sample comparisons. Use when you want
    DESeq2-style normalisation without running the full DESeq2 pipeline.

  - `"none"` — no normalisation; the count matrix is used as supplied.
    For pre-normalised inputs (e.g., TPM, FPKM) where an additional
    normalisation step would be redundant.

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

  Logical; only applies to `test_method = "deseq"`. Raw (MLE)
  fold-change estimates from DESeq2 can be unrealistically large for
  genes with very few reads, making them noisy to rank or plot. Setting
  `shrink.lfc = TRUE` corrects this by applying empirical Bayes
  shrinkage via `DESeq2::lfcShrink(type = "apeglm")` (requires the
  `apeglm` package): estimates from low-count genes are pulled towards
  zero, giving more trustworthy rankings and cleaner volcano plots.
  P-values and FDR are not affected — they always come from `results()`.
  Default is `FALSE`, but `TRUE` is strongly recommended for any result
  used in a figure or table.

- return_long:

  Logical. Only used when `test_method` is `"paired.ttest"` or
  `"unpaired.ttest"`. If `TRUE`, returns a named list with `$results`
  (the DE tibble) and `$long_data` (the per-gene, per-animal table used
  for the test). Default is `FALSE`.

- ngenes.out:

  Number of top genes (sorted by p-value) to include in the output when
  `kable.out = TRUE`. Ignored when `genes.filter` is supplied. Default
  is `20`.

- genes.filter:

  Optional character vector of gene names (matching the `Gene` column of
  the output) to retain. When supplied, only the listed genes are
  returned and `ngenes.out` is ignored. Default is `NULL`.

- kable.out:

  Logical. If `TRUE`, returns a `kableExtra` HTML table of the top
  `ngenes.out` genes instead of the full tibble. Default is `FALSE`.

- filter:

  Logical. If `TRUE` (default), low-expression genes are removed using
  [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
  before fitting any model. This is appropriate for raw count matrices.
  Set `filter = FALSE` when passing a pre-normalised matrix (e.g., RPKM,
  TPM) via `norm.method = "none"`, because `filterByExpr` expects
  integer counts and will incorrectly discard genes whose normalised
  values fall below the count threshold.

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

**Region-aware parsing**: When `region_name` is also supplied and
`sample_df = NULL`, the parser uses `region_name` to identify and strip
the brain-region token from each column name. Matching is
case-insensitive and the extracted region is stored in a `region` column
of the auto-built sample metadata (which is then used for the usual
region-filtering step). Column names that do *not* contain the supplied
`region_name` will retain the unrecognised token as part of the
treatment string; those samples are automatically excluded when the
function filters for `treatment_name`.

Region-aware examples (`region_name = "POA"`, `ip_level = "IP"`,
`input_level = "INPUT"`):

|                 |               |            |
|-----------------|---------------|------------|
| **Column name** | **treatment** | **region** |
| `b1inputPOA`    | `b`           | `POA`      |
| `a2ipPOA`       | `a`           | `POA`      |
| `b_1_input_POA` | `b`           | `POA`      |

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
