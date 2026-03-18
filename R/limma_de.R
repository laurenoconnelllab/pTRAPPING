#' Differential expression analysis using limma
#'
#' Runs a complete limma workflow on a feature-by-sample expression or
#' abundance matrix: builds a design matrix from a group vector, fits a linear
#' model, applies contrasts, and tests for differential expression using either
#' [limma::treat()] or [limma::eBayes()]. Returns a list containing
#' per-contrast DE tables, a combined table of all significant features, a
#' summary of features with inconsistent directions across contrasts, and
#' direction counts.
#'
#' @section Choosing `model_type`:
#' \describe{
#'   \item{`"pairwise"`}{Fits `model.matrix(~ 0 + group)`. No intercept; every
#'     group gets its own coefficient. Use when you want explicit pairwise
#'     contrasts (e.g. `groupA - groupB`). Requires a contrast matrix.}
#'   \item{`"reference"`}{Fits `model.matrix(~ group)`. The first factor level
#'     is the intercept / reference group; other coefficients represent
#'     differences from it. Simpler when you have one natural reference
#'     condition and only a few comparisons.}
#' }
#'
#' @section Choosing `test_method`:
#' \describe{
#'   \item{`"treat"`}{Calls [limma::treat()], which tests against a **minimum
#'     fold-change threshold** (`lfc_threshold`) rather than zero. Recommended
#'     for proteomics or metabolomics data where you want to ensure biological
#'     as well as statistical significance. Results are extracted with
#'     [limma::topTreat()].}
#'   \item{`"eBayes"`}{Calls [limma::eBayes()], the standard empirical Bayes
#'     moderation. Tests against a null of logFC = 0. Appropriate for RNA-seq
#'     (after `voom`) or any analysis where no minimum fold-change is required.
#'     Results are extracted with [limma::topTable()].}
#' }
#'
#' @param expr_mat A numeric matrix with features (genes / proteins) as rows
#'   and samples as columns. Row names must be feature identifiers and column
#'   names must be sample IDs.
#' @param group A character or factor vector of group labels, one element per
#'   column of `expr_mat`. The order must match the columns of `expr_mat`.
#' @param contrast_mat A contrast matrix produced by
#'   [limma::makeContrasts()], **or** a character vector of contrast strings
#'   (e.g. `c("groupA-groupB", "groupB-groupC")`) which will be passed to
#'   [limma::makeContrasts()] internally. Contrast names must match the column
#'   names of the design matrix (e.g. `"groupA"`, `"groupB"`, ...).
#' @param coefs An integer vector specifying which columns of `contrast_mat` to
#'   extract results for. If `NULL` (default), all contrasts are used.
#' @param model_type One of `"pairwise"` (no intercept, `~ 0 + group`,
#'   recommended for multi-group designs with explicit contrasts) or
#'   `"reference"` (first level as reference, `~ group`). See the
#'   *Choosing model_type* section. Default is `"pairwise"`.
#' @param test_method One of `"treat"` or `"eBayes"`. See the *Choosing
#'   test_method* section for guidance. Default is `"treat"`.
#' @param lfc_threshold Minimum absolute log2 fold change. Passed to
#'   `treat(lfc = lfc_threshold)` when `test_method = "treat"`, and used as a
#'   filter threshold for both methods. Default is `1`.
#' @param fdr_threshold Maximum adjusted p-value for a feature to be
#'   classified as differentially expressed. Default is `0.05`.
#' @param feature_col Name of the column that will hold feature identifiers
#'   (i.e. the row names of `expr_mat`) in the output tables. Default is
#'   `"feature"`.
#' @param ngenes_out Number of top features to display when
#'   `kable_out = TRUE`. Default is `20`.
#' @param kable_out Logical. If `TRUE`, returns an HTML `kableExtra` table for
#'   a single contrast (indexed by `kable_coef`) instead of the full list.
#'   Default is `FALSE`.
#' @param kable_coef Integer giving the position within `coefs` of the
#'   contrast to display when `kable_out = TRUE`. Default is `1` (first
#'   selected contrast).
#'
#' @return When `kable_out = FALSE` (default), a named list with:
#' \describe{
#'   \item{`de_list`}{A named list of tibbles, one per selected contrast,
#'     containing only DE features (filtered by `lfc_threshold` and
#'     `fdr_threshold`). Each tibble includes `feature_col`, `logFC`,
#'     `adj.P.Val`, `contrast`, and `diffexpressed` (`"UP"` or `"DOWN"`).}
#'   \item{`all_de`}{All per-contrast DE tibbles from `de_list` combined into
#'     one tibble with [dplyr::bind_rows()]. A feature appears once per
#'     contrast in which it is significant.}
#'   \item{`shared`}{Tibble of features that are DE in **opposite directions**
#'     across different contrasts (min logFC < 0 and max logFC > 0). Useful
#'     for flagging features whose direction of change is not consistent.}
#'   \item{`unique`}{Tibble counting consistently DE features by direction
#'     (`"UP"` or `"DOWN"`) — features absent from `shared`.}
#'   \item{`fit`}{The fitted `MArrayLM` object returned by `treat()` or
#'     `eBayes()`, for downstream use (e.g. calling `topTreat()` or
#'     `topTable()` with additional contrast columns).}
#' }
#' When `kable_out = TRUE`, an HTML `kableExtra` table of the top
#' `ngenes_out` features for the selected contrast.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # --- pairwise design with treat() (proteomics / metabolomics) ---
#' group_vec <- rep(c("groupA", "groupB", "groupC"), each = 3)
#'
#' # build the contrast matrix for selected pairwise comparisons
#' df_design <- data.frame(group = factor(group_vec))
#' design_mat <- model.matrix(~ 0 + group, data = df_design)
#' cont_mat <- limma::makeContrasts(
#'   groupA - groupB,
#'   groupA - groupC,
#'   levels = design_mat
#' )
#'
#' res <- limma_de(
#'   expr_mat      = my_matrix,
#'   group         = group_vec,
#'   contrast_mat  = cont_mat,
#'   coefs         = 1:2,
#'   test_method   = "treat",
#'   lfc_threshold = 1,
#'   fdr_threshold = 0.05
#' )
#'
#' res$de_list    # named list, one tibble per contrast
#' res$all_de     # all contrasts combined
#' res$shared  # features with inconsistent direction
#' res$unique  # count of consistently UP / DOWN features
#'
#' # --- single comparison, reference group design with eBayes() ---
#' res2 <- limma_de(
#'   expr_mat     = my_matrix,
#'   group        = group_vec,
#'   contrast_mat = "groupgroupB",   # coefficient name in ~ group design
#'   model_type   = "reference",
#'   test_method  = "eBayes"
#' )
#' }
#'
#' @importFrom limma lmFit contrasts.fit treat eBayes makeContrasts
#'   topTreat topTable
#' @importFrom dplyr filter mutate case_when bind_rows group_by summarise
#'   distinct count if_else across slice_head all_of
#' @importFrom tibble as_tibble tibble
#' @importFrom rlang .data
#' @importFrom stats model.matrix
#' @importFrom knitr kable
#' @importFrom kableExtra kable_classic row_spec column_spec

limma_de <- function(
  expr_mat,
  group,
  contrast_mat,
  coefs         = NULL,
  model_type    = c("pairwise", "reference"),
  test_method   = c("treat", "eBayes"),
  lfc_threshold = 1,
  fdr_threshold = 0.05,
  feature_col   = "feature",
  ngenes_out    = 20,
  kable_out     = FALSE,
  kable_coef    = 1
) {
  model_type  <- match.arg(model_type)
  test_method <- match.arg(test_method)

  # --- build design matrix ---------------------------------------------------
  # "pairwise": ~ 0 + group  → one coefficient per group, no intercept.
  #             Needed when you want explicit A-vs-B contrasts.
  # "reference": ~ group     → first level is the intercept/reference group.
  #             Simpler when you have one natural baseline condition.
  df_design <- data.frame(group = factor(group))

  design <- if (model_type == "pairwise") {
    model.matrix(~ 0 + group, data = df_design)
  } else {
    model.matrix(~ group, data = df_design)
  }

  # --- build contrast matrix if character strings were supplied --------------
  # Users can pass makeContrasts() output directly, or a character vector of
  # contrast expressions (e.g. c("groupA-groupB", "groupA-groupC")).
  if (is.character(contrast_mat)) {
    contrast_mat <- makeContrasts(
      contrasts = contrast_mat,
      levels    = design
    )
  }

  # --- resolve which contrasts to extract ------------------------------------
  if (is.null(coefs)) {
    coefs <- seq_len(ncol(contrast_mat))
  }

  # --- fit linear model and apply contrasts ----------------------------------
  fit      <- lmFit(expr_mat, design)
  fit_cont <- contrasts.fit(fit, contrast_mat)

  # --- apply test method -----------------------------------------------------
  # treat():  tests H0: |logFC| <= lfc_threshold. Use for proteomics /
  #           metabolomics where a minimum effect size is biologically required.
  #           Extract results with topTreat().
  # eBayes(): tests H0: logFC = 0. Standard choice for RNA-seq (post-voom)
  #           or when no minimum fold-change is required.
  #           Extract results with topTable().
  if (test_method == "treat") {
    fit_test <- treat(fit_cont, lfc = lfc_threshold)
    top_fn   <- function(k) {
      topTreat(fit_test, coef = k, adjust.method = "BH", number = Inf)
    }
  } else {
    fit_test <- eBayes(fit_cont)
    top_fn   <- function(k) {
      topTable(fit_test, coef = k, adjust.method = "BH", number = Inf)
    }
  }

  # --- extract and filter DE results for each selected contrast --------------
  de_list <- lapply(coefs, function(k) {
    contrast_name <- colnames(contrast_mat)[k]
    top_fn(k) |>
      as_tibble(rownames = feature_col) |>
      mutate(
        contrast      = contrast_name,
        diffexpressed = case_when(
          .data$logFC >  lfc_threshold & .data$adj.P.Val < fdr_threshold ~ "UP",
          .data$logFC < -lfc_threshold & .data$adj.P.Val < fdr_threshold ~ "DOWN",
          TRUE ~ "NO"
        )
      ) |>
      filter(.data$diffexpressed != "NO")
  })
  names(de_list) <- colnames(contrast_mat)[coefs]

  # --- combine all DE results into one table ---------------------------------
  all_de <- bind_rows(de_list)

  # --- shared & consistent summaries ----------------------------------------
  # shared:     features DE in opposite directions in different contrasts
  #             (min logFC < 0 AND max logFC > 0). Flag these — their direction
  #             of change is not consistent across the experiment.
  # consistent: features excluded from shared, counted by direction.
  if (nrow(all_de) > 0) {

    direction_summary <- all_de |>
      distinct(across(all_of(c(feature_col, "logFC"))), .keep_all = TRUE) |>
      group_by(across(all_of(feature_col))) |>
      summarise(
        min_fc      = min(.data$logFC),
        max_fc      = max(.data$logFC),
        n_contrasts = dplyr::n(),
        .groups     = "drop"
      )

    shared <- direction_summary |>
      filter(.data$min_fc < 0 & .data$max_fc > 0)

    unique <- direction_summary |>
      filter(!(.data$min_fc < 0 & .data$max_fc > 0)) |>
      mutate(direction = if_else(.data$max_fc > 0, "UP", "DOWN")) |>
      count(.data$direction, name = "n_features")

  } else {
    shared <- tibble()
    unique <- tibble()
  }

  # --- kable output for a single contrast ------------------------------------
  if (kable_out) {
    k <- coefs[[kable_coef]]
    return(
      top_fn(k) |>
        as_tibble(rownames = feature_col) |>
        slice_head(n = ngenes_out) |>
        kable(
          digits     = 2,
          table.attr = 'data-quarto-disable-processing="true"',
          "html"
        ) |>
        kable_classic(full_width = FALSE, html_font = "Cambria") |>
        row_spec(0, italic = TRUE, bold = TRUE) |>
        column_spec(1, italic = FALSE, bold = TRUE)
    )
  }

  # --- return list -----------------------------------------------------------
  list(
    de_list = de_list,
    all_de  = all_de,
    shared  = shared,
    unique  = unique,
    fit     = fit_test
  )
}
