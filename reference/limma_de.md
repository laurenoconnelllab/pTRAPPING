# Differential expression analysis using limma

Runs a complete limma workflow on a feature-by-sample expression or
abundance matrix: builds a design matrix from a group vector, fits a
linear model, applies contrasts, and tests for differential expression
using either
[`limma::treat()`](https://rdrr.io/pkg/limma/man/ebayes.html) or
[`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html). Returns
a list containing per-contrast DE tables, a combined table of all
significant features, a summary of features with inconsistent directions
across contrasts, and direction counts.

## Usage

``` r
limma_de(
  expr_mat,
  group,
  contrast_mat,
  coefs = NULL,
  model_type = c("pairwise", "reference"),
  test_method = c("treat", "eBayes"),
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  feature_col = "feature",
  ngenes_out = 20,
  kable_out = FALSE,
  kable_coef = 1
)
```

## Arguments

- expr_mat:

  A numeric matrix with features (genes / proteins) as rows and samples
  as columns. Row names must be feature identifiers and column names
  must be sample IDs.

- group:

  A character or factor vector of group labels, one element per column
  of `expr_mat`. The order must match the columns of `expr_mat`.

- contrast_mat:

  A contrast matrix produced by
  [`limma::makeContrasts()`](https://rdrr.io/pkg/limma/man/makeContrasts.html),
  **or** a character vector of contrast strings (e.g.
  `c("groupA-groupB", "groupB-groupC")`) which will be passed to
  [`limma::makeContrasts()`](https://rdrr.io/pkg/limma/man/makeContrasts.html)
  internally. Contrast names must match the column names of the design
  matrix (e.g. `"groupA"`, `"groupB"`, ...).

- coefs:

  An integer vector specifying which columns of `contrast_mat` to
  extract results for. If `NULL` (default), all contrasts are used.

- model_type:

  One of `"pairwise"` (no intercept, `~ 0 + group`, recommended for
  multi-group designs with explicit contrasts) or `"reference"` (first
  level as reference, `~ group`). See the *Choosing model_type* section.
  Default is `"pairwise"`.

- test_method:

  One of `"treat"` or `"eBayes"`. See the *Choosing test_method* section
  for guidance. Default is `"treat"`.

- lfc_threshold:

  Minimum absolute log2 fold change. Passed to
  `treat(lfc = lfc_threshold)` when `test_method = "treat"`, and used as
  a filter threshold for both methods. Default is `1`.

- fdr_threshold:

  Maximum adjusted p-value for a feature to be classified as
  differentially expressed. Default is `0.05`.

- feature_col:

  Name of the column that will hold feature identifiers (i.e. the row
  names of `expr_mat`) in the output tables. Default is `"feature"`.

- ngenes_out:

  Number of top features to display when `kable_out = TRUE`. Default is
  `20`.

- kable_out:

  Logical. If `TRUE`, returns an HTML `kableExtra` table for a single
  contrast (indexed by `kable_coef`) instead of the full list. Default
  is `FALSE`.

- kable_coef:

  Integer giving the position within `coefs` of the contrast to display
  when `kable_out = TRUE`. Default is `1` (first selected contrast).

## Value

When `kable_out = FALSE` (default), a named list with:

- `de_list`:

  A named list of tibbles, one per selected contrast, containing only DE
  features (filtered by `lfc_threshold` and `fdr_threshold`). Each
  tibble includes `feature_col`, `logFC`, `adj.P.Val`, `contrast`, and
  `diffexpressed` (`"UP"` or `"DOWN"`).

- `all_de`:

  All per-contrast DE tibbles from `de_list` combined into one tibble
  with
  [`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html).
  A feature appears once per contrast in which it is significant.

- `shared`:

  Tibble of features that are DE in **opposite directions** across
  different contrasts (min logFC \< 0 and max logFC \> 0). Useful for
  flagging features whose direction of change is not consistent.

- `unique`:

  Tibble counting consistently DE features by direction (`"UP"` or
  `"DOWN"`) — features absent from `shared`.

- `fit`:

  The fitted `MArrayLM` object returned by `treat()` or `eBayes()`, for
  downstream use (e.g. calling `topTreat()` or `topTable()` with
  additional contrast columns).

When `kable_out = TRUE`, an HTML `kableExtra` table of the top
`ngenes_out` features for the selected contrast.

## Choosing `model_type`

- `"pairwise"`:

  Fits `model.matrix(~ 0 + group)`. No intercept; every group gets its
  own coefficient. Use when you want explicit pairwise contrasts (e.g.
  `groupA - groupB`). Requires a contrast matrix.

- `"reference"`:

  Fits `model.matrix(~ group)`. The first factor level is the intercept
  / reference group; other coefficients represent differences from it.
  Simpler when you have one natural reference condition and only a few
  comparisons.

## Choosing `test_method`

- `"treat"`:

  Calls [`limma::treat()`](https://rdrr.io/pkg/limma/man/ebayes.html),
  which tests against a **minimum fold-change threshold**
  (`lfc_threshold`) rather than zero. Recommended for proteomics or
  metabolomics data where you want to ensure biological as well as
  statistical significance. Results are extracted with
  [`limma::topTreat()`](https://rdrr.io/pkg/limma/man/toptable.html).

- `"eBayes"`:

  Calls [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html),
  the standard empirical Bayes moderation. Tests against a null of logFC
  = 0. Appropriate for RNA-seq (after `voom`) or any analysis where no
  minimum fold-change is required. Results are extracted with
  [`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html).

## Examples

``` r
if (FALSE) { # \dontrun{
# --- pairwise design with treat() (proteomics / metabolomics) ---
group_vec <- rep(c("groupA", "groupB", "groupC"), each = 3)

# build the contrast matrix for selected pairwise comparisons
df_design <- data.frame(group = factor(group_vec))
design_mat <- model.matrix(~ 0 + group, data = df_design)
cont_mat <- limma::makeContrasts(
  groupA - groupB,
  groupA - groupC,
  levels = design_mat
)

res <- limma_de(
  expr_mat      = my_matrix,
  group         = group_vec,
  contrast_mat  = cont_mat,
  coefs         = 1:2,
  test_method   = "treat",
  lfc_threshold = 1,
  fdr_threshold = 0.05
)

res$de_list    # named list, one tibble per contrast
res$all_de     # all contrasts combined
res$shared  # features with inconsistent direction
res$unique  # count of consistently UP / DOWN features

# --- single comparison, reference group design with eBayes() ---
res2 <- limma_de(
  expr_mat     = my_matrix,
  group        = group_vec,
  contrast_mat = "groupgroupB",   # coefficient name in ~ group design
  model_type   = "reference",
  test_method  = "eBayes"
)
} # }
```
