# --- internal helpers ---------------------------------------------------------
# (not exported; used only by ptrap_de)

# Parse a single sample column name into treatment / block / fraction.
# Strategy: strip separators (_, -, .), split at every alpha<->digit boundary,
# then classify each token as fraction keyword | numeric block | treatment.
#
# Examples (all produce the same output fields):
#   "b1input"   -> treatment="b",  block="1", fraction=input_level
#   "nb1ip"     -> treatment="nb", block="1", fraction=ip_level
#   "Nb_IP_1"   -> treatment="Nb", block="1", fraction=ip_level
#   "B.3.INPUT" -> treatment="B",  block="3", fraction=input_level
.parse_one_sample <- function(nm, ip_level, input_level) {
  # Step 1: Remove separators.
  cleaned <- gsub("[-_. ]+", "", nm)

  # Step 2: Inject a canonical separator ("_") just before each fraction
  # keyword, so they are always split off as their own token regardless of
  # whether they are fused with the treatment name.
  # INPUT must be handled before IP so the longer keyword wins.
  # Use case-insensitive replacement via perl = TRUE + (?i).
  cleaned <- gsub(
    paste0("(?i)(?=", input_level, ")"),
    "_",
    cleaned,
    perl = TRUE
  )
  cleaned <- gsub(
    paste0("(?i)(?=", ip_level, ")"),
    "_",
    cleaned,
    perl = TRUE
  )

  # Step 3: Re-strip any double separators produced above.
  cleaned <- gsub("_+", "_", cleaned)
  cleaned <- gsub("^_|_$", "", cleaned)

  # Step 4: Split on the canonical separator.
  parts <- unlist(strsplit(cleaned, "_"))

  # Step 5: Further split each part at digit<->letter transitions.
  parts <- unlist(lapply(parts, function(p) {
    unlist(strsplit(p, "(?<=\\d)(?=\\D)|(?<=\\D)(?=\\d)", perl = TRUE))
  }))
  parts <- parts[nzchar(parts)]

  parts_l <- tolower(parts)

  # fraction: matches ip_level, input_level, or "in" (short alias for INPUT)
  frac_kws <- unique(tolower(c(ip_level, input_level, "in")))
  frac_idx <- which(parts_l %in% frac_kws)
  # block: purely numeric token
  block_idx <- which(grepl("^\\d+$", parts))
  # treatment: everything that is neither fraction nor block
  treat_idx <- setdiff(seq_along(parts), c(frac_idx, block_idx))

  if (length(frac_idx) == 0L) {
    stop(
      "Cannot identify IP/INPUT fraction in column name '",
      nm,
      "'. ",
      "Column names must contain '",
      ip_level,
      "' or '",
      input_level,
      "' (case-insensitive). See ?ptrap_de for naming conventions."
    )
  }
  if (length(block_idx) == 0L) {
    stop(
      "Cannot identify a replicate number in column name '",
      nm,
      "'. ",
      "Column names must contain a digit identifying the biological replicate. ",
      "See ?ptrap_de for naming conventions."
    )
  }

  fraction <- if (parts_l[frac_idx[1L]] == tolower(ip_level)) {
    ip_level
  } else {
    input_level
  }
  block <- parts[block_idx[1L]]
  treatment <- paste(parts[treat_idx], collapse = "")

  list(sample = nm, treatment = treatment, block = block, fraction = fraction)
}

# Build a tibble of sample metadata by parsing every column name.
.build_sample_df_from_cols <- function(col_names, ip_level, input_level) {
  parsed <- lapply(
    col_names,
    .parse_one_sample,
    ip_level = ip_level,
    input_level = input_level
  )
  tibble::tibble(
    sample = vapply(parsed, `[[`, character(1L), "sample"),
    treatment = vapply(parsed, `[[`, character(1L), "treatment"),
    block = vapply(parsed, `[[`, character(1L), "block"),
    fraction = vapply(parsed, `[[`, character(1L), "fraction")
  )
}

# Build a design formula, inserting covariate terms additively.
# For DESeq2 (method = "deseq"): fraction_col must be last so that results()
#   auto-detects the contrast. Order: block + covariates + fraction.
# For edgeR / voom: order does not matter; covariates appended at the end.
#   Order: fraction + block + covariates.
.build_formula <- function(method, fraction_col, block_col, covariates) {
  cov_terms <- if (!is.null(covariates)) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }
  rhs_parts <- if (method == "deseq") {
    c(block_col, cov_terms, fraction_col) # fraction last
  } else {
    c(fraction_col, block_col, cov_terms) # fraction first
  }
  rhs_parts <- rhs_parts[!vapply(rhs_parts, is.null, logical(1L))]
  as.formula(paste("~", paste(rhs_parts, collapse = " + ")))
}


# --- main function ------------------------------------------------------------

#' Perform differential expression analysis for PhosphoTRAP data using edgeR,
#' a paired t-test, an unpaired t-test, or voom/limma
#'
#' @title Differential expression analysis for PhosphoTRAP data
#'
#' @description
#' Compares IP vs INPUT fractions for a specified treatment condition using
#' one of six statistical approaches (see `test_method`). Both `sample_df`
#' and `gene_ids` are optional: the function can derive sample metadata
#' automatically from the column names of `counts_mat`, and gene identifiers
#' from its first column.
#'
#' * `"LRT"` and `"QLF"` use edgeR's GLM framework with a multi-factor
#'   design (`~ fraction + block`), accounting for the paired animal structure
#'   via the blocking variable. Both methods automatically apply **edgeR's TMM
#'   (trimmed mean of M values) normalization** via `normLibSizes()`, which
#'   adjusts effective library sizes (not the count matrix itself) for use in
#'   the GLM as offsets. `norm.method` is not applicable and is ignored.
#'   Recommended with **4 or more replicates** and raw count data.
#' * `"deseq"` uses the **DESeq2** pipeline with a multi-factor design
#'   (`~ block + fraction`; fraction goes last so `results()` auto-tests it).
#'   DESeq2 applies the **median of ratios normalization** internally via
#'   `DESeq()`. `norm.method` is not applicable and is ignored. Recommended
#'   when you prefer the negative-binomial shrinkage estimators of DESeq2 over
#'   edgeR. **Important:** by default, `logFC` is the **maximum likelihood
#'   estimate (MLE)** from `results()`, which can be noisy and large in
#'   magnitude for low-count genes. Set `shrink.lfc = TRUE` to apply
#'   empirical Bayes LFC shrinkage via `lfcShrink()` (requires the
#'   `apeglm` package), which is strongly recommended before ranking genes
#'   or producing volcano plots for publication.
#' * `"voom"` uses the **limma-voom** pipeline with a paired design
#'   (`~ fraction + block`) and empirical Bayes moderation. Voom internally
#'   transforms counts to **log2-CPM** (using TMM-adjusted library sizes) and
#'   derives **precision weights** from the mean-variance trend; these two
#'   steps are integral to the method. `norm.method` is not applicable and is
#'   ignored.
#' * `"paired.ttest"` runs a per-gene **paired t-test between IP and INPUT**
#'   values across the experimental repeats, following Tan et al. (2016)
#'   *Cell* 167, 47-59 \doi{10.1016/j.cell.2016.08.028}: *"p value for each
#'   gene was calculated as the paired t-test between input and
#'   immunoprecipitated RPKM values from the three experimental repeats."*
#'   The count scale for the test is controlled by `norm.method` (Tan et al.
#'   used `"RPKM"`; all three options are valid; default is `"CPM"`). `logFC`
#'   is reported as \eqn{\log_2(\bar{\mathrm{IP}} / \bar{\mathrm{INPUT}})}.
#'   Recommended for PhosphoTRAP experiments with **3 replicates**, where
#'   GLM-based dispersion estimates are unreliable.
#' * `"unpaired.ttest"` compares **fold enrichments (IP/INPUT) between two
#'   treatment groups** using a Welch (unpaired) t-test on per-animal log2(FE)
#'   values. Use this to test whether ribosomal association differs between
#'   conditions. Requires exactly two treatments; use `treatment_name` and
#'   `control_name` to specify which is which. The count scale is controlled
#'   by `norm.method` (default is `"CPM"`; all three options are valid).
#'
#' @section Automatic column-name parsing:
#' When `sample_df = NULL`, the column names of `counts_mat` (excluding the
#' first, gene-ID column) are parsed to build sample metadata automatically.
#' Each column name must encode three pieces of information, in any order and
#' with any combination of separators (`_`, `-`, `.`, space) or no separator
#' at all:
#' \describe{
#'   \item{Treatment}{One or more letters identifying the experimental group.}
#'   \item{Replicate number}{A digit identifying the biological replicate.}
#'   \item{Fraction}{`ip_level` or `input_level` (case-insensitive); `"in"`
#'     is also accepted as a short alias for the INPUT fraction.}
#' }
#' Valid column name examples (default `ip_level = "IP"`,
#' `input_level = "INPUT"`):
#' \tabular{ll}{
#'   **Column name** \tab **Parsed as** \cr
#'   `b1input`, `b1ip`    \tab treatment = `b`, block = `1` \cr
#'   `nb2INPUT`           \tab treatment = `nb`, block = `2` \cr
#'   `Nb_IP_1`            \tab treatment = `Nb`, block = `1` \cr
#'   `B_3_INPUT`          \tab treatment = `B`, block = `3` \cr
#'   `PB.2.ip`            \tab treatment = `PB`, block = `2` \cr
#'   `Trim_10-INPUT`      \tab treatment = `Trim`, block = `10` \cr
#'   `SOL1INPUT`          \tab treatment = `SOL`, block = `1` \cr
#' }
#'
#' @param counts_mat A counts matrix in one of two formats:
#'   * **Matrix** -- numeric, genes x samples; column names are sample IDs;
#'     gene IDs are in `rownames` or supplied via `gene_ids`.
#'   * **Data frame / tibble** -- first column is a character vector of gene
#'     IDs; remaining columns are numeric counts with sample names as column
#'     names. When `sample_df = NULL`, column names must follow the naming
#'     convention described in the *Automatic column-name parsing* section.
#' @param sample_df Optional data frame of sample metadata. When `NULL`
#'   (default), metadata is parsed automatically from the column names of
#'   `counts_mat`. When provided, the arguments `sample_col`, `fraction_col`,
#'   `block_col`, `region_col`, and `treatment_col` specify which columns to
#'   use (falling back to their defaults if the names match).
#' @param gene_ids Optional character vector of gene identifiers corresponding
#'   to the rows of `counts_mat`. When `NULL` (default), gene IDs are
#'   extracted from the first column of `counts_mat` (if it is a data frame)
#'   or from `rownames(counts_mat)` (if it is a matrix).
#' @param region_name The brain region to subset and analyse (e.g., `"POA"`).
#'   When `region_col` is `NULL` (default) but `region_name` is supplied, the
#'   function automatically scans `sample_df` for a non-structural column whose
#'   values include `region_name`, and uses it as `region_col` if exactly one
#'   such column is found (a message is shown). Can be `NULL` when the data
#'   contain only a single region.
#' @param treatment_name The treatment condition whose samples will be
#'   **subsetted** for the IP vs INPUT comparison (e.g., `"pb"` to analyse
#'   only the pair-bonded samples). For `"unpaired.ttest"`, this is the
#'   numerator condition in the differential fold enrichment
#'   (DFE = mean_FE_treatment / mean_FE_control).
#'
#'   When `NULL` (default), the value is derived automatically from
#'   `sample_df` (or from the parsed column names when `sample_df = NULL`):
#'   if a single treatment is found the function proceeds silently; if
#'   multiple treatments are found an error asks the user to specify which
#'   one to analyse. For `"unpaired.ttest"` with exactly two treatments and
#'   both `treatment_name` and `control_name` as `NULL`, assignments are made
#'   alphabetically (first = control, second = treatment) with a message.
#' @param control_name The **reference / control** treatment used only by
#'   `"unpaired.ttest"`. Serves as the denominator when computing differential
#'   fold enrichment (DFE = mean_FE_treatment / mean_FE_control). When `NULL`
#'   (default), the function auto-detects the control as the non-`treatment_name`
#'   group when exactly two treatments are present.
#' @param sample_col Name of the column in `sample_df` whose values match the
#'   column names of `counts_mat`. Ignored when `sample_df = NULL` (auto-set
#'   internally). Default is `"sample"`.
#' @param fraction_col Name of the column in `sample_df` that distinguishes IP
#'   from INPUT fractions. Default is `"fraction"`.
#' @param block_col Name of the column in `sample_df` used as the blocking /
#'   pairing variable (e.g., individual animal or tube). For `"LRT"` / `"QLF"`,
#'   `"deseq"`, and `"voom"` it enters the design matrix; for `"paired.ttest"`
#'   and `"unpaired.ttest"` it aligns each animal's IP with its own INPUT.
#'   Default is `"tube"`.
#' @param region_col Name of the column in `sample_df` containing brain region
#'   labels (e.g., `"BrainRegion"`). When `NULL` (default) and `region_name`
#'   is supplied, the column is auto-detected (see `region_name`). Set
#'   explicitly when the auto-detection is ambiguous or when you prefer to be
#'   explicit. Leave `NULL` when the data contain only one brain region.
#' @param treatment_col Name of the column in `sample_df` containing treatment
#'   labels. Default is `"Treatment"`.
#' @param ip_level The value in `fraction_col` that identifies the IP fraction.
#'   Default is `"IP"`.
#' @param input_level The value in `fraction_col` that identifies the INPUT
#'   fraction (reference level). Default is `"INPUT"`.
#' @param covariates Optional character vector of covariate names to include in
#'   the design formula. Default is `NULL`.
#' @param lfc_threshold Minimum absolute log2 fold change required to classify
#'   a gene as differentially expressed. Default is `1`.
#' @param fdr_threshold Maximum FDR (or adjusted p-value) allowed to classify
#'   a gene as differentially expressed. Default is `0.05`.
#' @param test_method Statistical method to use. One of:
#'   * `"LRT"` -- edgeR likelihood ratio test via `glmFit` + `glmLRT`
#'     (default). TMM normalisation is applied internally; `norm.method` is
#'     not applicable. Suitable for >= 4 replicates and raw count data.
#'   * `"QLF"` -- edgeR quasi-likelihood F-test via `glmQLFit` + `glmQLFTest`.
#'     TMM normalisation is applied internally; `norm.method` is not
#'     applicable. More conservative than `"LRT"`.
#'   * `"deseq"` -- DESeq2 pipeline with a multi-factor design
#'     (`~ block + fraction`). Median of ratios normalisation is applied
#'     internally by `DESeq()`; `norm.method` is not applicable. Results
#'     obtained via `results(dds, name = "fraction_IP_vs_INPUT")`. `logFC`
#'     is MLE by default; use `shrink.lfc = TRUE` for shrunken estimates.
#'   * `"voom"` -- limma-voom pipeline with a paired design
#'     (`~ fraction + block`) and empirical Bayes moderation via
#'     `lmFit` + `eBayes` + `topTable`. Log2-CPM transformation and voom
#'     precision weights are applied internally; `norm.method` is not
#'     applicable.
#'   * `"paired.ttest"` -- per-gene paired t-test between IP and INPUT values
#'     across replicates, following Tan et al. (2016). The count scale is
#'     controlled by `norm.method` (default `"CPM"`; Tan et al. used
#'     `"RPKM"`). Best suited for n = 3 replicates. P-values adjusted with BH.
#'   * `"unpaired.ttest"` -- per-gene Welch unpaired t-test comparing log2(FE)
#'     values between two treatment groups. The count scale is controlled by
#'     `norm.method` (default `"CPM"`). Requires exactly two treatments; use
#'     `treatment_name` and `control_name` to specify which group is which.
#' @param norm.method Count normalisation method used by the t-test branches
#'   (`"paired.ttest"` and `"unpaired.ttest"`). Ignored for `"LRT"`, `"QLF"`,
#'   `"voom"`, and `"deseq"`, which each handle normalisation internally (a
#'   message is shown if you explicitly supply `norm.method` for those
#'   methods). Default is `"CPM"`. One of:
#'   * `"CPM"` (default) -- counts per million, computed from edgeR's
#'     TMM-adjusted effective library sizes (`lib.size x norm.factors`) via
#'     `edgeR::cpm()`. Because TMM adjusts *effective library sizes* rather
#'     than the count matrix directly, CPM values here reflect both
#'     sequencing depth and TMM normalisation. Suitable for comparisons
#'     between replicates of the same sample group.
#'   * `"RPKM"` -- reads per kilobase per million: CPM further divided by gene
#'     length in kilobases, via `edgeR::rpkm()`. Requires `gene.length`.
#'     Suitable for within-sample comparisons between genes (the scale used
#'     by Tan et al. 2016 for PhosphoTRAP). **Not** recommended for
#'     between-sample comparisons.
#'   * `"mratios"` -- DESeq2 **median of ratios** normalisation, computed via
#'     `estimateSizeFactors()` and retrieved with `counts(dds, normalized = TRUE)`.
#'     Unlike CPM, this method directly rescales the count matrix by
#'     sample-specific size factors, making it suitable for between-sample
#'     comparisons. Use this when you want DESeq2-style normalisation for the
#'     t-test without running the full DESeq2 pipeline.
#'   * `"none"` -- no normalisation is applied; the count matrix is used
#'     as supplied. Intended for pre-normalised matrices (e.g., TPM, FPKM,
#'     or any user-normalised values) where an additional normalisation step
#'     would be redundant or distorting.
#' @param gene.length Named numeric vector of gene lengths in base pairs
#'   (names must match gene IDs). Required when `norm.method = "RPKM"`;
#'   ignored otherwise. Default is `NULL`.
#' @param prior.count A small count added to IP and INPUT values before
#'   computing log2 ratios (FE = (IP + prior.count) / (INPUT + prior.count)),
#'   to avoid log(0). Default is `1`. Set to `0` if values are already
#'   normalised and guaranteed positive. Used only by `"paired.ttest"` and
#'   `"unpaired.ttest"`.
#' @param shrink.lfc Logical. Only used when `test_method = "deseq"`. If
#'   `TRUE`, log2 fold changes are shrunk using empirical Bayes estimation via
#'   `DESeq2::lfcShrink(type = "apeglm")`, which requires the
#'   [apeglm](https://bioconductor.org/packages/apeglm/) package. Shrunken
#'   LFCs are more reliable for ranking genes and for visualisation because
#'   they pull noisy, high-variance estimates (typically from low-count genes)
#'   towards zero. P-values and FDR are not affected -- they always come from
#'   `results()`. Default is `FALSE` (MLE fold changes are returned), but
#'   `TRUE` is strongly recommended for any result used in a figure or table.
#' @param return_long Logical. Only used when `test_method` is
#'   `"paired.ttest"` or `"unpaired.ttest"`. If `TRUE`, returns a named list
#'   with `$results` (the DE tibble) and `$long_data` (the per-gene,
#'   per-animal table used for the test). Default is `FALSE`.
#' @param ngenes.out Number of top genes (sorted by p-value) to include in the
#'   output when `kable.out = TRUE`. Default is `20`.
#' @param kable.out Logical. If `TRUE`, returns a `kableExtra` HTML table of
#'   the top `ngenes.out` genes instead of the full tibble. Default is
#'   `FALSE`.
#'
#' @return When `kable.out = FALSE` (default), a tibble with one row per gene
#'   sorted by p-value. Columns depend on `test_method`:
#'
#'   * **`"LRT"` / `"QLF"`**: `Gene`, `logFC`, `logCPM`, `LR` or `F`,
#'     `PValue`, `FDR`, treatment label, optional region label,
#'     `diffexpressed` (`"UP"`, `"DOWN"`, `"NO"`).
#'   * **`"deseq"`**: `Gene`, `baseMean`, `logFC` (MLE log2FoldChange by
#'     default; empirical Bayes shrunken estimate when `shrink.lfc = TRUE`),
#'     `lfcSE`, `stat` (Wald statistic; present with MLE, absent when
#'     `shrink.lfc = TRUE` as apeglm replaces it with a posterior estimate),
#'     `PValue`, `FDR` (Benjamini-Hochberg from `results()`; unaffected by
#'     shrinkage), treatment label, optional region label, `diffexpressed`.
#'     Genes with `NA` p-values (low-count outliers flagged by DESeq2) are
#'     retained and sorted to the bottom.
#'   * **`"paired.ttest"`**: always returns a **named list** with two
#'     components:
#'     - `$results` — tibble with `Gene`,
#'       `logFC` (\eqn{\log_2(\bar{\mathrm{IP}} / \bar{\mathrm{INPUT}})}),
#'       `t_statistic`, `PValue`, `FDR`, treatment label, optional region
#'       label, `diffexpressed`.
#'     - `$fe` — wide tibble with `Gene` plus one `FE_<block>` column per
#'       animal, containing the per-animal linear IP/INPUT fold enrichment
#'       values used in the test.
#'     - `$long_data` — (only when `return_long = TRUE`) a per-gene,
#'       per-animal tibble with columns `Gene`, block column, `ip_count`,
#'       `input_count`, `FE`.
#'   * **`"unpaired.ttest"`**: `Gene`, `logFC` (log2(mean_FE_treatment /
#'     mean_FE_control)), `diff_FE`, `mean_FE_<treatment>`,
#'     `mean_FE_<control>`, `t_statistic`, `df` (Welch degrees of freedom),
#'     `PValue`, `FDR`, treatment label, optional region label,
#'     `diffexpressed`. With `return_long = TRUE`, returns a named list:
#'     `$results` and `$long_data` (per-gene, per-animal, per-treatment
#'     tibble with columns `Gene`, block column, treatment column, `FE`,
#'     `log2_FE`).
#'   * **`"voom"`**: `Gene`, `logFC`, `AveExpr`, `t`, `B`, `PValue`, `FDR`,
#'     treatment label, optional region label, `diffexpressed`.
#'
#'   When `kable.out = TRUE`, an HTML `kableExtra` table of the top
#'   `ngenes.out` genes.
#'
#' @references
#' Tan, C.L., Cooke, E.K., Leib, D.E., Lin, Y.C., Daly, G.E., Zimmerman,
#' C.A., and Knight, Z.A. (2016). Warm-Sensitive Neurons that Control Body
#' Temperature. *Cell* 167, 47-59.
#' \doi{10.1016/j.cell.2016.08.028}
#'
#' Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold
#' change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15,
#' 550. \doi{10.1186/s13059-014-0550-8}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## ---- Option A: simplest call -- auto-parse from column names ---------------
#' # counts_mat is a data frame where:
#' #   col 1      = gene IDs
#' #   col 2+     = samples named like "b1input", "b1ip", "nb2input", etc.
#' counts <- read.table("counts.txt", header = TRUE)
#'
#' # single treatment in the matrix -- treatment_name auto-detected
#' res <- ptrap_de(counts_mat = counts, test_method = "paired.ttest")
#'
#' # multiple treatments -- specify which one to analyze
#' res_b <- ptrap_de(counts_mat = counts, treatment_name = "b")
#'
#' ## ---- Option B: provide sample_df explicitly (original workflow) -----------
#' res_lrt <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb"
#' )
#'
#' # quasi-likelihood F-test
#' res_qlf <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb",
#'   test_method    = "QLF"
#' )
#'
#' # paired t-test -- also return long-format paired table
#' res_pt <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb",
#'   test_method    = "paired.ttest",
#'   return_long    = TRUE
#' )
#' res_pt$results    # DE tibble
#' res_pt$long_data  # per-gene, per-animal pairs
#'
#' # unpaired t-test -- compare fold enrichment between two treatments
#' res_upt <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   treatment_name = "pb",
#'   control_name   = "nb",
#'   test_method    = "unpaired.ttest"
#' )
#'
#' # voom/limma pipeline
#' res_voom <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   region_name    = "POA",
#'   treatment_name = "pb",
#'   test_method    = "voom"
#' )
#'
#' # paired t-test with CPM normalization
#' res_cpm <- ptrap_de(
#'   counts_mat     = counts_mat,
#'   sample_df      = sample_df,
#'   gene_ids       = gene_ids,
#'   treatment_name = "pb",
#'   test_method    = "paired.ttest",
#'   norm.method    = "CPM"
#' )
#' }
#'
#' @importFrom edgeR DGEList filterByExpr normLibSizes estimateDisp glmFit glmLRT glmQLFit glmQLFTest topTags cpm rpkm
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink estimateSizeFactors counts
#' @importFrom dplyr filter mutate case_when relocate slice_head arrange bind_rows bind_cols rename
#' @importFrom rlang .data :=
#' @importFrom tibble as_tibble tibble rownames_to_column
#' @importFrom stats as.formula model.matrix p.adjust pt sd var
#' @importFrom knitr kable
#' @importFrom kableExtra kable_classic row_spec column_spec
#' @importFrom cli cli_abort

ptrap_de <- function(
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
  test_method = c(
    "LRT",
    "QLF",
    "paired.ttest",
    "unpaired.ttest",
    "voom",
    "deseq"
  ),
  norm.method = c("CPM", "RPKM", "mratios", "none"),
  gene.length = NULL,
  prior.count = 1,
  shrink.lfc = FALSE,
  return_long = FALSE,
  ngenes.out = 20,
  kable.out = FALSE,
  filter = TRUE
) {
  # ---- Step 1: match arguments -----------------------------------------------
  # Capture whether norm.method was explicitly supplied before match.arg()
  # consumes the default, so we can warn only when the user actually set it.
  norm_method_set <- !missing(norm.method)
  test_method <- match.arg(test_method)
  norm.method <- match.arg(norm.method)

  # ---- Step 1b: inform users when norm.method is irrelevant ------------------
  # LRT, QLF, voom, and deseq each handle normalisation internally; norm.method
  # is silently ignored for these methods. A message is shown only when the
  # user explicitly passed a norm.method argument.
  if (norm_method_set && test_method %in% c("LRT", "QLF")) {
    message(
      "norm.method = '",
      norm.method,
      "' is ignored for test_method = '",
      test_method,
      "'. edgeR's ",
      test_method,
      " applies TMM (trimmed mean ",
      "of M values) normalisation internally via normLibSizes(): TMM adjusts ",
      "effective library sizes (not the count matrix itself), and those ",
      "adjusted sizes enter the GLM as offsets."
    )
  }
  if (norm_method_set && test_method == "deseq") {
    message(
      "norm.method = '",
      norm.method,
      "' is ignored for test_method = ",
      "'deseq'. DESeq2 applies the median of ratios normalisation method ",
      "automatically inside DESeq()."
    )
  }
  if (norm_method_set && test_method == "voom") {
    message(
      "norm.method = '",
      norm.method,
      "' is ignored for test_method = ",
      "'voom'. limma-voom transforms counts to log2-CPM (using TMM-adjusted ",
      "library sizes) and derives precision weights from the mean-variance ",
      "trend; these steps are integral to the voom pipeline."
    )
  }

  # ---- Step 2: resolve gene IDs and ensure counts_mat is a numeric matrix ----
  if (is.data.frame(counts_mat)) {
    first_col <- counts_mat[[1L]]
    if (is.character(first_col) || is.factor(first_col)) {
      # first column contains gene identifiers
      if (is.null(gene_ids)) {
        gene_ids <- as.character(first_col)
      }
      counts_mat <- as.matrix(counts_mat[, -1L, drop = FALSE])
      mode(counts_mat) <- "numeric"
    } else {
      counts_mat <- as.matrix(counts_mat)
    }
  }

  if (is.null(gene_ids)) {
    rn <- rownames(counts_mat)
    if (!is.null(rn)) {
      gene_ids <- rn
    } else {
      stop(
        "Cannot determine gene IDs. Please either:\n",
        "  (a) include gene IDs as the first column of `counts_mat`, or\n",
        "  (b) pass them explicitly via `gene_ids`."
      )
    }
  }

  # ---- Step 3: build sample metadata if not supplied -------------------------
  if (is.null(sample_df)) {
    sample_df <- .build_sample_df_from_cols(
      colnames(counts_mat),
      ip_level = ip_level,
      input_level = input_level
    )
    # override column name arguments to match the auto-built data frame
    sample_col <- "sample"
    fraction_col <- "fraction"
    block_col <- "block"
    treatment_col <- "treatment"

    message(
      "Auto-parsed sample metadata from column names.\n",
      "  Treatments : ",
      paste(sort(unique(sample_df$treatment)), collapse = ", "),
      "\n",
      "  Blocks     : ",
      paste(sort(unique(sample_df$block)), collapse = ", "),
      "\n",
      "  Fractions  : ",
      paste(sort(unique(sample_df$fraction)), collapse = ", ")
    )
  }

  # ---- Step 4: resolve treatment_name (and control_name for unpaired.ttest) --
  if (test_method == "unpaired.ttest") {
    unique_tx <- unique(sample_df[[treatment_col]])

    if (is.null(treatment_name) && is.null(control_name)) {
      if (length(unique_tx) == 2L) {
        sorted_tx <- sort(unique_tx)
        control_name <- sorted_tx[1L]
        treatment_name <- sorted_tx[2L]
        message(
          "Auto-assigned treatments alphabetically: ",
          "control = '",
          control_name,
          "', treatment = '",
          treatment_name,
          "'."
        )
      } else {
        stop(
          "For 'unpaired.ttest', exactly two treatments are required. ",
          "Found: ",
          paste(unique_tx, collapse = ", "),
          ". ",
          "Please specify both `treatment_name` and `control_name`."
        )
      }
    } else if (is.null(treatment_name) && !is.null(control_name)) {
      remaining <- setdiff(unique_tx, control_name)
      if (length(remaining) == 1L) {
        treatment_name <- remaining
        message("Auto-assigned treatment_name = '", treatment_name, "'.")
      } else {
        stop(
          "Cannot auto-resolve `treatment_name`. ",
          "After excluding control_name = '",
          control_name,
          "', ",
          "found treatments: ",
          paste(remaining, collapse = ", "),
          ". ",
          "Please specify `treatment_name` explicitly."
        )
      }
    } else if (!is.null(treatment_name) && is.null(control_name)) {
      remaining <- setdiff(unique_tx, treatment_name)
      if (length(remaining) == 1L) {
        control_name <- remaining
        message("Auto-assigned control_name = '", control_name, "'.")
      } else {
        stop(
          "Cannot auto-resolve `control_name`. ",
          "After excluding treatment_name = '",
          treatment_name,
          "', ",
          "found treatments: ",
          paste(remaining, collapse = ", "),
          ". ",
          "Please specify `control_name` explicitly."
        )
      }
    }
    # Both are now set; validate they differ
    if (identical(treatment_name, control_name)) {
      stop(
        "`treatment_name` and `control_name` must be different; both are '",
        treatment_name,
        "'."
      )
    }
  } else {
    # All methods other than unpaired.ttest: resolve single treatment_name
    if (is.null(treatment_name)) {
      unique_tx <- unique(sample_df[[treatment_col]])
      if (length(unique_tx) == 1L) {
        treatment_name <- unique_tx
        message("Using the only treatment found: '", treatment_name, "'")
      } else {
        stop(
          "Multiple treatments found (",
          paste(unique_tx, collapse = ", "),
          "). Please specify `treatment_name`."
        )
      }
    }
  }

  # ---- Step 5: resolve region_col / region_name ------------------------------

  # 5a. region_name supplied but region_col not: try to auto-detect the column
  #     by scanning sample_df for a column that contains region_name as a value.
  if (!is.null(region_name) && is.null(region_col)) {
    structural <- c(sample_col, fraction_col, block_col, treatment_col)
    candidate_cols <- Filter(
      function(col) {
        col %in%
          names(sample_df) &&
          !col %in% structural &&
          any(as.character(sample_df[[col]]) == region_name)
      },
      names(sample_df)
    )
    if (length(candidate_cols) == 1L) {
      region_col <- candidate_cols
      message(
        "Auto-detected region column '",
        region_col,
        "' ",
        "from region_name = '",
        region_name,
        "'."
      )
    } else if (length(candidate_cols) == 0L) {
      warning(
        "`region_name = '",
        region_name,
        "'` was supplied but no column in ",
        "`sample_df` contains '",
        region_name,
        "' as a value. ",
        "Region filtering will NOT be applied. If your data span multiple ",
        "brain regions, set `region_col` to the appropriate column name."
      )
    } else {
      stop(
        "`region_name = '",
        region_name,
        "'` was supplied but multiple columns ",
        "contain '",
        region_name,
        "' as a value: ",
        paste(candidate_cols, collapse = ", "),
        ". ",
        "Please specify `region_col` explicitly."
      )
    }
  }

  # 5b. region_col supplied but region_name not: auto-resolve if only one region
  if (!is.null(region_col) && is.null(region_name)) {
    unique_reg <- unique(sample_df[[region_col]])
    if (length(unique_reg) == 1L) {
      region_name <- unique_reg
      message("Using the only region found: '", region_name, "'")
    } else {
      stop(
        "Multiple regions found (",
        paste(unique_reg, collapse = ", "),
        "). Please specify `region_name`."
      )
    }
  }

  # ---- Step 6: subset region_samples -----------------------------------------
  if (test_method == "unpaired.ttest") {
    if (!is.null(region_col)) {
      region_samples <- sample_df |>
        filter(
          .data[[region_col]] == region_name,
          .data[[treatment_col]] %in% c(treatment_name, control_name)
        )
    } else {
      region_samples <- sample_df |>
        filter(.data[[treatment_col]] %in% c(treatment_name, control_name))
    }
  } else {
    if (!is.null(region_col)) {
      region_samples <- sample_df |>
        filter(
          .data[[region_col]] == region_name,
          .data[[treatment_col]] == treatment_name
        )
    } else {
      region_samples <- sample_df |>
        filter(.data[[treatment_col]] == treatment_name)
    }
  }

  region_samples <- region_samples |>
    mutate(
      !!fraction_col := factor(
        .data[[fraction_col]],
        levels = c(input_level, ip_level)
      )
    )

  if (nrow(region_samples) == 0L) {
    stop(
      "No samples found for treatment '",
      treatment_name,
      "'",
      if (!is.null(region_col)) {
        paste0(" in region '", region_name, "'")
      } else {
        ""
      },
      ". Check that these values exist in `sample_df`."
    )
  }

  # ---- Step 7: validate sample-column matching --------------------------------
  missing_samples <- setdiff(region_samples[[sample_col]], colnames(counts_mat))
  if (length(missing_samples) > 0L) {
    stop(
      "The following samples are in `sample_df` but not in `counts_mat` columns: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  # ---- Step 8: subset counts and set row names --------------------------------
  counts_region <- counts_mat[, region_samples[[sample_col]], drop = FALSE]
  rownames(counts_region) <- gene_ids

  # ---- Step 8b: validate covariates ------------------------------------------
  if (!is.null(covariates)) {
    if (!is.character(covariates)) {
      cli::cli_abort(
        "{.arg covariates} must be a character vector of column names in \\
        {.arg sample_df}, e.g. {.code covariates = c(\"batch\", \"sex\")}."
      )
    }
    missing_covs <- setdiff(covariates, names(region_samples))
    if (length(missing_covs) > 0L) {
      cli::cli_abort(
        c(
          "The following {.arg covariates} column{?s} {?was/were} not found \\
          in the sample metadata:",
          "x" = "{.val {missing_covs}}",
          "i" = "Available columns: {.val {names(region_samples)}}"
        )
      )
    }
    structural_cols <- c(fraction_col, block_col, treatment_col)
    clashing <- intersect(covariates, structural_cols)
    if (length(clashing) > 0L) {
      cli::cli_abort(
        c(
          "{.arg covariates} must not include columns already used as \\
          structural terms in the model.",
          "x" = "{.val {clashing}} {?is/are} already used as \\
          {.arg fraction_col}, {.arg block_col}, or {.arg treatment_col}."
        )
      )
    }
    if (test_method %in% c("paired.ttest", "unpaired.ttest")) {
      message(
        "covariates = c(",
        paste0("'", covariates, "'", collapse = ", "),
        ") ",
        "cannot be applied to test_method = '",
        test_method,
        "', which does ",
        "not use a design matrix. Covariates will be ignored. ",
        "To control for covariates, use test_method = 'LRT', 'QLF', ",
        "'voom', or 'deseq'."
      )
      covariates <- NULL
    }
  }

  # ---- Step 9: common DGEList preprocessing (ALL methods) --------------------
  dge <- DGEList(counts = counts_region)
  dge$genes <- data.frame(Gene = gene_ids)

  if (filter) {
    keep <- filterByExpr(dge, group = region_samples[[fraction_col]])
    dge <- dge[keep, , keep.lib.sizes = FALSE]
  }

  # TMM normalization -- skipped for "deseq" (DESeq2 handles its own normalization)
  if (test_method != "deseq") {
    dge <- normLibSizes(dge)
  }

  # ---- Step 10: compute working counts for t-test branches -------------------
  if (test_method %in% c("paired.ttest", "unpaired.ttest")) {
    if (norm.method == "CPM") {
      # CPM computed from TMM-adjusted effective library sizes
      wc <- cpm(dge, log = FALSE)
    } else if (norm.method == "RPKM") {
      # RPKM = CPM / (gene length in kb); uses TMM-adjusted library sizes
      if (is.null(gene.length)) {
        stop(
          "'gene.length' must be a named numeric vector when ",
          "norm.method = \"RPKM\"."
        )
      }
      gl <- gene.length[rownames(dge)]
      if (anyNA(gl)) {
        stop(
          "Some gene IDs in the filtered data are missing from ",
          "names(gene.length)."
        )
      }
      wc <- rpkm(dge, gene.length = gl, log = FALSE)
    } else if (norm.method == "mratios") {
      # "mratios": DESeq2 median of ratios -- size factors only, no full pipeline
      coldata_mr <- as.data.frame(region_samples)
      rownames(coldata_mr) <- region_samples[[sample_col]]
      dds_mr <- DESeqDataSetFromMatrix(
        countData = round(dge$counts),
        colData = coldata_mr,
        design = ~1 # design not needed for normalisation only
      )
      dds_mr <- estimateSizeFactors(dds_mr)
      wc <- counts(dds_mr, normalized = TRUE)
    } else {
      # "none": use the count matrix as-is, without any normalisation
      wc <- dge$counts
    }
  }

  # ---- paired t-test branch (Tan et al. 2016) --------------------------------
  if (test_method == "paired.ttest") {
    # separate and sort IP / INPUT by block to ensure correct pairing
    ip_samples <- region_samples |>
      filter(.data[[fraction_col]] == ip_level) |>
      arrange(.data[[block_col]])

    input_samples <- region_samples |>
      filter(.data[[fraction_col]] == input_level) |>
      arrange(.data[[block_col]])

    if (
      !identical(
        sort(as.character(ip_samples[[block_col]])),
        sort(as.character(input_samples[[block_col]]))
      )
    ) {
      stop(
        "IP and INPUT samples do not share the same '",
        block_col,
        "' values. ",
        "Each animal/tube must have exactly one IP and one INPUT sample."
      )
    }

    ip_mat <- wc[, ip_samples[[sample_col]], drop = FALSE]
    input_mat <- wc[, input_samples[[sample_col]], drop = FALSE]
    n_reps <- ncol(ip_mat)
    block_ids <- as.character(ip_samples[[block_col]])

    # Paired t-test between IP and INPUT values across n experimental repeats,
    # following Tan et al. 2016 (Cell 167, 47-59):
    # "p value for each gene was calculated as the paired t test between input
    #  and immunoprecipitated RPKM values from the three experimental repeats."
    # norm.method controls the count scale (Tan et al. used "RPKM").
    paired_d <- ip_mat - input_mat # internal: pairing by column (animal)
    t_stat <- rowMeans(paired_d) / (apply(paired_d, 1L, sd) / sqrt(n_reps))
    p_val <- 2 * pt(-abs(t_stat), df = n_reps - 1L)

    # logFC for reporting: mean of per-animal log2 fold enrichments.
    # Each animal contributes equally, consistent with the paired design.
    logFC_vec <- rowMeans(log2(
      (ip_mat + prior.count) / (input_mat + prior.count)
    ))

    # per-animal fold enrichment: FE_<block_id> — kept separate from results
    FE_mat <- (ip_mat + prior.count) / (input_mat + prior.count)
    FE_cols <- as.data.frame(FE_mat)
    colnames(FE_cols) <- paste0("FE_", block_ids)

    # fe tibble: Gene + one FE_<block> column per animal
    fe <- bind_cols(tibble(Gene = rownames(ip_mat)), as_tibble(FE_cols)) |>
      mutate(
        log2_mean_FE = log2(rowMeans(as.matrix(FE_cols)))
      )

    # assemble clean results tibble (no FE columns)
    results <- tibble(
      Gene = rownames(ip_mat),
      logFC = as.numeric(logFC_vec),
      t_statistic = as.numeric(t_stat),
      PValue = as.numeric(p_val)
    ) |>
      mutate(
        FDR = p.adjust(.data$PValue, method = "BH"),
        !!treatment_col := treatment_name,
        diffexpressed = case_when(
          .data$logFC > lfc_threshold & .data$FDR < fdr_threshold ~ "UP",
          .data$logFC < -lfc_threshold & .data$FDR < fdr_threshold ~ "DOWN",
          TRUE ~ "NO"
        )
      ) |>
      arrange(.data$PValue)

    # long-format paired table (returned when return_long = TRUE)
    long_data <- bind_rows(
      lapply(seq_len(n_reps), function(j) {
        tibble(
          Gene = rownames(ip_mat),
          !!block_col := block_ids[j],
          ip_count = as.numeric(ip_mat[, j]),
          input_count = as.numeric(input_mat[, j]),
          FE = as.numeric(FE_mat[, j])
        )
      })
    ) |>
      arrange(.data$Gene)

    if (!is.null(region_col)) {
      results <- results |> mutate(!!region_col := region_name)
    }

    if (kable.out) {
      return(
        results |>
          slice_head(n = ngenes.out) |>
          kable(
            digits = 3L,
            table.attr = 'data-quarto-disable-processing="true"',
            "html"
          ) |>
          kable_classic(full_width = FALSE, html_font = "Cambria") |>
          row_spec(0L, italic = TRUE, bold = TRUE) |>
          column_spec(1L, italic = FALSE, bold = TRUE)
      )
    }

    # Always return a named list: $results (DE tibble) + $fe (per-animal FE).
    # $long_data is added when return_long = TRUE.
    out <- list(results = results, fe = fe)
    if (return_long) {
      out$long_data <- long_data
    }
    return(out)
  }

  # ---- unpaired t-test branch ------------------------------------------------
  if (test_method == "unpaired.ttest") {
    # Separate samples by treatment group
    trt_samples <- region_samples |>
      filter(.data[[treatment_col]] == treatment_name)
    ctrl_samples <- region_samples |>
      filter(.data[[treatment_col]] == control_name)

    # Within each group, split into IP and INPUT, sort by block_col
    trt_ip <- trt_samples |>
      filter(.data[[fraction_col]] == ip_level) |>
      arrange(.data[[block_col]])

    trt_input <- trt_samples |>
      filter(.data[[fraction_col]] == input_level) |>
      arrange(.data[[block_col]])

    ctrl_ip <- ctrl_samples |>
      filter(.data[[fraction_col]] == ip_level) |>
      arrange(.data[[block_col]])

    ctrl_input <- ctrl_samples |>
      filter(.data[[fraction_col]] == input_level) |>
      arrange(.data[[block_col]])

    # Validate block IDs match within each group
    if (
      !identical(
        sort(as.character(trt_ip[[block_col]])),
        sort(as.character(trt_input[[block_col]]))
      )
    ) {
      stop(
        "For treatment '",
        treatment_name,
        "', IP and INPUT samples do not ",
        "share the same '",
        block_col,
        "' values. ",
        "Each animal/tube must have exactly one IP and one INPUT sample."
      )
    }
    if (
      !identical(
        sort(as.character(ctrl_ip[[block_col]])),
        sort(as.character(ctrl_input[[block_col]]))
      )
    ) {
      stop(
        "For control '",
        control_name,
        "', IP and INPUT samples do not ",
        "share the same '",
        block_col,
        "' values. ",
        "Each animal/tube must have exactly one IP and one INPUT sample."
      )
    }

    # Extract matrices from working counts
    trt_ip_mat <- wc[, trt_ip[[sample_col]], drop = FALSE]
    trt_input_mat <- wc[, trt_input[[sample_col]], drop = FALSE]
    ctrl_ip_mat <- wc[, ctrl_ip[[sample_col]], drop = FALSE]
    ctrl_input_mat <- wc[, ctrl_input[[sample_col]], drop = FALSE]

    # Compute FE (linear) and log2_FE per animal per group
    FE_trt <- (trt_ip_mat + prior.count) / (trt_input_mat + prior.count)
    FE_ctrl <- (ctrl_ip_mat + prior.count) / (ctrl_input_mat + prior.count)
    log2_FE_trt <- log2(FE_trt)
    log2_FE_ctrl <- log2(FE_ctrl)

    # Vectorised Welch unpaired t-test on log2_FE per gene
    n_t <- ncol(log2_FE_trt)
    n_c <- ncol(log2_FE_ctrl)

    # Power warning: Welch t-test with n <= 3 per group leaves only 1-2
    # degrees of freedom per group. Over thousands of genes, BH correction
    # becomes extremely conservative and will often yield no significant hits.
    if (n_t < 4L || n_c < 4L) {
      message(
        "Low-power warning ('unpaired.ttest'): '",
        treatment_name,
        "' has n = ",
        n_t,
        " and '",
        control_name,
        "' has n = ",
        n_c,
        " replicates.\n",
        "  A Welch t-test with ",
        min(n_t, n_c),
        " replicates per group ",
        "yields only ",
        min(n_t, n_c) - 1L,
        " degree(s) of freedom per group,\n",
        "  making BH correction over ",
        nrow(log2_FE_trt),
        " genes highly conservative -- few or no genes may reach\n",
        "  significance at the default thresholds. Consider:\n",
        "    * relaxing fdr_threshold (e.g. 0.10 or 0.20) or lfc_threshold;\n",
        "    * increasing biological replicates (n >= 4 recommended);\n",
        "    * using test_method = 'LRT', 'QLF', 'voom', or 'deseq', which\n",
        "      borrow information across genes to improve power at small n."
      )
    }

    mean_t <- rowMeans(log2_FE_trt)
    mean_c <- rowMeans(log2_FE_ctrl)
    var_t <- apply(log2_FE_trt, 1L, var)
    var_c <- apply(log2_FE_ctrl, 1L, var)
    se <- sqrt(var_t / n_t + var_c / n_c)
    t_stat <- (mean_t - mean_c) / se
    df_w <- (var_t / n_t + var_c / n_c)^2 /
      ((var_t / n_t)^2 / (n_t - 1L) + (var_c / n_c)^2 / (n_c - 1L))
    p_val <- 2 * pt(-abs(t_stat), df = df_w)

    # Mean linear FE per condition and differential FE (descriptive columns)
    mean_FE_trt <- rowMeans(FE_trt)
    mean_FE_ctrl <- rowMeans(FE_ctrl)
    diff_FE <- mean_FE_trt / mean_FE_ctrl
    # logFC: mean(log2_FE_treatment) - mean(log2_FE_control), consistent with
    # the Welch t-test which is run on per-animal log2 FE values.
    logFC <- mean_t - mean_c

    # Build results tibble
    results <- tibble(
      Gene = rownames(trt_ip_mat),
      logFC = as.numeric(logFC),
      diff_FE = as.numeric(diff_FE),
      !!paste0("mean_FE_", treatment_name) := as.numeric(mean_FE_trt),
      !!paste0("mean_FE_", control_name) := as.numeric(mean_FE_ctrl),
      t_statistic = as.numeric(t_stat),
      df = as.numeric(df_w),
      PValue = as.numeric(p_val)
    ) |>
      mutate(
        FDR = p.adjust(.data$PValue, method = "BH"),
        !!treatment_col := treatment_name,
        diffexpressed = case_when(
          .data$logFC > lfc_threshold & .data$FDR < fdr_threshold ~ "UP",
          .data$logFC < -lfc_threshold & .data$FDR < fdr_threshold ~ "DOWN",
          TRUE ~ "NO"
        )
      ) |>
      arrange(.data$PValue)

    # Build long_data (per-animal log2_FE values for both groups)
    long_data <- bind_rows(
      lapply(seq_len(n_t), function(j) {
        tibble(
          Gene = rownames(FE_trt),
          !!block_col := as.character(trt_ip[[block_col]][j]),
          !!treatment_col := treatment_name,
          FE = as.numeric(FE_trt[, j]),
          log2_FE = as.numeric(log2_FE_trt[, j])
        )
      }),
      lapply(seq_len(n_c), function(j) {
        tibble(
          Gene = rownames(FE_ctrl),
          !!block_col := as.character(ctrl_ip[[block_col]][j]),
          !!treatment_col := control_name,
          FE = as.numeric(FE_ctrl[, j]),
          log2_FE = as.numeric(log2_FE_ctrl[, j])
        )
      })
    ) |>
      arrange(.data$Gene)

    if (!is.null(region_col)) {
      results <- results |> mutate(!!region_col := region_name)
    }

    if (kable.out) {
      return(
        results |>
          slice_head(n = ngenes.out) |>
          kable(
            digits = 3L,
            table.attr = 'data-quarto-disable-processing="true"',
            "html"
          ) |>
          kable_classic(full_width = FALSE, html_font = "Cambria") |>
          row_spec(0L, italic = TRUE, bold = TRUE) |>
          column_spec(1L, italic = FALSE, bold = TRUE)
      )
    }

    if (return_long) {
      return(list(results = results, long_data = long_data))
    }
    return(results)
  }

  # ---- DESeq2 branch ---------------------------------------------------------
  if (test_method == "deseq") {
    # Build colData from region_samples; rownames must match countData columns
    coldata_ds <- as.data.frame(region_samples)
    rownames(coldata_ds) <- region_samples[[sample_col]]

    # Ensure block is a factor so DESeq2 treats it as categorical
    coldata_ds[[block_col]] <- factor(coldata_ds[[block_col]])

    # fraction already factored in Step 6: levels = c(input_level, ip_level)
    # Block first, optional covariates next, fraction last ->
    # results() auto-tests the last variable
    design_ds <- .build_formula("deseq", fraction_col, block_col, covariates)

    dds <- DESeqDataSetFromMatrix(
      countData = round(dge$counts), # DESeq2 requires integer counts
      colData = coldata_ds,
      design = design_ds
    )
    dds <- DESeq(dds, quiet = TRUE)

    # Coefficient name follows DESeq2 convention: <var>_<level>_vs_<ref>
    res_name <- paste0(fraction_col, "_", ip_level, "_vs_", input_level)
    res_deseq <- results(dds, name = res_name)

    # Optional: empirical Bayes LFC shrinkage via apeglm (recommended for
    # publication figures and gene ranking). P-values / FDR are unaffected.
    if (shrink.lfc) {
      if (!requireNamespace("apeglm", quietly = TRUE)) {
        stop(
          "shrink.lfc = TRUE requires the 'apeglm' package.\n",
          "Install it with: BiocManager::install('apeglm')"
        )
      }
      res_deseq <- lfcShrink(
        dds,
        coef = res_name,
        type = "apeglm",
        quiet = TRUE
      )
    }

    results <- as.data.frame(res_deseq) |>
      rownames_to_column("Gene") |>
      as_tibble() |>
      rename(logFC = "log2FoldChange", PValue = "pvalue", FDR = "padj") |>
      mutate(
        !!treatment_col := treatment_name,
        diffexpressed = case_when(
          !is.na(.data$logFC) &
            .data$logFC > lfc_threshold &
            !is.na(.data$FDR) &
            .data$FDR < fdr_threshold ~ "UP",
          !is.na(.data$logFC) &
            .data$logFC < -lfc_threshold &
            !is.na(.data$FDR) &
            .data$FDR < fdr_threshold ~ "DOWN",
          TRUE ~ "NO"
        )
      ) |>
      arrange(.data$PValue) |>
      relocate("Gene")

    if (!is.null(region_col)) {
      results <- results |> mutate(!!region_col := region_name)
    }

    if (kable.out) {
      return(
        results |>
          slice_head(n = ngenes.out) |>
          kable(
            digits = 2L,
            table.attr = 'data-quarto-disable-processing="true"',
            "html"
          ) |>
          kable_classic(full_width = FALSE, html_font = "Cambria") |>
          row_spec(0L, italic = TRUE, bold = TRUE) |>
          column_spec(1L, italic = FALSE, bold = TRUE)
      )
    }

    return(results)
  }

  # ---- voom branch -----------------------------------------------------------
  if (test_method == "voom") {
    design_formula <- .build_formula(
      "voom",
      fraction_col,
      block_col,
      covariates
    )
    design <- model.matrix(design_formula, data = region_samples)

    # For RPKM: apply rpkm() and pass matrix; for CPM/none: pass dge directly
    if (norm.method == "RPKM") {
      if (is.null(gene.length)) {
        stop(
          "'gene.length' must be a named numeric vector when ",
          "norm.method = \"RPKM\"."
        )
      }
      gl <- gene.length[rownames(dge)]
      if (anyNA(gl)) {
        stop("Some gene IDs missing from names(gene.length).")
      }
      warning(
        "norm.method = 'RPKM' with voom is non-standard; ",
        "voom is designed for count-scale data."
      )
      v <- voom(rpkm(dge, gene.length = gl, log = FALSE), design, plot = FALSE)
    } else {
      v <- voom(dge, design, plot = FALSE)
    }

    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    lrt_coef <- paste0(fraction_col, ip_level)

    results <- topTable(
      fit,
      coef = lrt_coef,
      number = Inf,
      adjust.method = "BH"
    ) |>
      as_tibble() |>
      rename(PValue = "P.Value", FDR = "adj.P.Val") |>
      mutate(
        !!treatment_col := treatment_name,
        diffexpressed = case_when(
          .data$logFC > lfc_threshold & .data$FDR < fdr_threshold ~ "UP",
          .data$logFC < -lfc_threshold & .data$FDR < fdr_threshold ~ "DOWN",
          TRUE ~ "NO"
        )
      )

    if (!is.null(region_col)) {
      results <- results |> mutate(!!region_col := region_name)
    }

    results <- results |> relocate("Gene")

    if (kable.out) {
      return(
        results |>
          slice_head(n = ngenes.out) |>
          kable(
            digits = 2L,
            table.attr = 'data-quarto-disable-processing="true"',
            "html"
          ) |>
          kable_classic(full_width = FALSE, html_font = "Cambria") |>
          row_spec(0L, italic = TRUE, bold = TRUE) |>
          column_spec(1L, italic = FALSE, bold = TRUE)
      )
    }

    return(results)
  }

  # ---- edgeR branch (LRT / QLF) ----------------------------------------------
  design_formula <- .build_formula("edger", fraction_col, block_col, covariates)
  design <- model.matrix(design_formula, data = region_samples)

  dge <- estimateDisp(dge, design)
  lrt_coef <- paste0(fraction_col, ip_level)

  if (test_method == "LRT") {
    fit <- glmFit(dge, design)
    test <- glmLRT(fit, coef = lrt_coef)
  } else {
    fit <- glmQLFit(dge, design)
    test <- glmQLFTest(fit, coef = lrt_coef)
  }

  # topTags already carries the Gene column from dge$genes -- do not re-add it
  results <- topTags(test, n = Inf)$table |>
    as_tibble() |>
    mutate(
      !!treatment_col := treatment_name,
      diffexpressed = case_when(
        .data$logFC > lfc_threshold & .data$FDR < fdr_threshold ~ "UP",
        .data$logFC < -lfc_threshold & .data$FDR < fdr_threshold ~ "DOWN",
        TRUE ~ "NO"
      )
    )

  if (!is.null(region_col)) {
    results <- results |> mutate(!!region_col := region_name)
  }

  results <- results |> relocate("Gene")

  if (kable.out) {
    return(
      results |>
        slice_head(n = ngenes.out) |>
        kable(
          digits = 2L,
          table.attr = 'data-quarto-disable-processing="true"',
          "html"
        ) |>
        kable_classic(full_width = FALSE, html_font = "Cambria") |>
        row_spec(0L, italic = TRUE, bold = TRUE) |>
        column_spec(1L, italic = FALSE, bold = TRUE)
    )
  }

  return(results)
}
