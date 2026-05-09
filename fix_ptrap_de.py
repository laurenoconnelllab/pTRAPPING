"""Apply all changes to R/ptrap_de.R:
1. Replace .parse_one_sample() with the surround-substitute fix.
2. Add genes.filter param to roxygen, signature, and all 5 kable branches.
"""

import re

PATH = "R/ptrap_de.R"

with open(PATH, "r") as f:
    src = f.read()

# ---------------------------------------------------------------------------
# CHANGE 1: Replace .parse_one_sample() body (steps 2-3 only)
# Old: two-step lookahead injection
# New: single surround-substitution + "in" alias handler
# ---------------------------------------------------------------------------

OLD_PARSER = '''\
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

  # Step 3: Re-strip any double separators produced above.'''

NEW_PARSER = '''\
  # Step 2: Surround the fraction keyword with "_" on both sides, isolating it
  # as its own token regardless of its position in the string.
  # input_level is listed first in the alternation so PCRE's left-to-right
  # greedy matching picks the longer keyword when ip_level is a substring
  # of input_level (e.g. "IP" inside "INPUT").
  cleaned <- gsub(
    paste0("(?i)(", input_level, "|", ip_level, ")"),
    "_\\\\1_",
    cleaned,
    perl = TRUE
  )
  # Handle "in" alias: only when NOT already part of input_level.
  # After the substitution above, input_level tokens are flanked by "_", so
  # any bare "in" remaining is a genuine alias.
  cleaned <- gsub(
    "(?i)(?<![A-Za-z])(in)(?![A-Za-z])",
    "_\\\\1_",
    cleaned,
    perl = TRUE
  )

  # Step 3: Re-strip any double separators produced above.'''

assert OLD_PARSER in src, "Could not find OLD_PARSER block"
src = src.replace(OLD_PARSER, NEW_PARSER, 1)

# Also update the header comment block to document all orderings
OLD_HEADER = '''\
# Examples (all produce the same output fields):
#   "b1input"   -> treatment="b",  block="1", fraction=input_level
#   "nb1ip"     -> treatment="nb", block="1", fraction=ip_level
#   "Nb_IP_1"   -> treatment="Nb", block="1", fraction=ip_level
#   "B.3.INPUT" -> treatment="B",  block="3", fraction=input_level'''

NEW_HEADER = '''\
# Supports ALL orderings of treatment (T) / replicate (R) / fraction (F):
#   T-R-F: "PACAP1IP",      "PACAP_1_IP",    "PACAP.1-IP"
#   T-F-R: "PACAP_Input_1", "PACAP_IP_1",    "PACAPInput1"
#   R-T-F: "1_PACAP_IP",    "1PACAPINPUT"
#   R-F-T: "1_IP_PACAP",    "1_INPUT_PACAP"
#   F-T-R: "IP_PACAP_1",    "INPUT_PACAP_1"
#   F-R-T: "IP_1_PACAP",    "INPUT_1_PACAP"
#
# "in" (case-insensitive) is accepted as a short alias for input_level.'''

assert OLD_HEADER in src, "Could not find OLD_HEADER block"
src = src.replace(OLD_HEADER, NEW_HEADER, 1)

# ---------------------------------------------------------------------------
# CHANGE 2: Add genes.filter roxygen param
# ---------------------------------------------------------------------------

OLD_ROXYGEN = '''\
#\' @param ngenes.out Number of top genes (sorted by p-value) to include in the
#\'   output when `kable.out = TRUE`. Default is `20`.
#\' @param kable.out Logical. If `TRUE`, returns a `kableExtra` HTML table of
#\'   the top `ngenes.out` genes instead of the full tibble. Default is
#\'   `FALSE`.'''

NEW_ROXYGEN = '''\
#\' @param ngenes.out Number of top genes (sorted by p-value) to include in the
#\'   output when `kable.out = TRUE`. Ignored when `genes.filter` is supplied.
#\'   Default is `20`.
#\' @param genes.filter Optional character vector of gene names (matching the
#\'   `Gene` column of the output) to retain in the results. When supplied, only
#\'   the listed genes are returned and `ngenes.out` is ignored. Default is
#\'   `NULL` (no filtering).
#\' @param kable.out Logical. If `TRUE`, returns a `kableExtra` HTML table of
#\'   the top `ngenes.out` genes instead of the full tibble. Default is
#\'   `FALSE`.'''

assert OLD_ROXYGEN in src, "Could not find OLD_ROXYGEN block"
src = src.replace(OLD_ROXYGEN, NEW_ROXYGEN, 1)

# ---------------------------------------------------------------------------
# CHANGE 3: Add genes.filter to function signature
# ---------------------------------------------------------------------------

OLD_SIG = '  ngenes.out = 20,\n  kable.out = FALSE,'
NEW_SIG = '  ngenes.out = 20,\n  genes.filter = NULL,\n  kable.out = FALSE,'

assert OLD_SIG in src, "Could not find OLD_SIG"
src = src.replace(OLD_SIG, NEW_SIG, 1)

# ---------------------------------------------------------------------------
# CHANGE 4: Update all 5 kable branches
# Each branch has a unique digits= value and surrounding context.
# Pattern to find: the `if (kable.out)` block with slice_head
# We replace the kable.out block in each branch.
# ---------------------------------------------------------------------------

def replace_kable_branch(src, region_col_block, digits):
    """
    region_col_block: the text that appears just before `if (kable.out)` 
                      (used to uniquely identify which branch we're in)
    digits: the digits= value used in kable() for this branch
    """
    old_block = (
        region_col_block +
        "\n"
        "\n"
        "    if (kable.out) {\n"
        "      return(\n"
        "        results |\n"
        "          slice_head(n = ngenes.out) |\n"
        "          kable(\n"
        f"            digits = {digits}L,\n"
        '            table.attr = \'data-quarto-disable-processing="true"\',\n'
        '            "html"\n'
        "          ) |\n"
        "          kable_classic(full_width = FALSE, html_font = \"Cambria\") |\n"
        "          row_spec(0L, italic = TRUE, bold = TRUE) |\n"
        "          column_spec(1L, italic = FALSE, bold = TRUE)\n"
        "      )\n"
        "    }"
    )
    new_block = (
        region_col_block +
        "\n"
        "\n"
        "    if (!is.null(genes.filter)) {\n"
        "      results <- dplyr::filter(results, .data$Gene %in% genes.filter)\n"
        "    }\n"
        "\n"
        "    if (kable.out) {\n"
        "      tbl <- if (is.null(genes.filter)) slice_head(results, n = ngenes.out) else results\n"
        "      return(\n"
        "        tbl |\n"
        "          kable(\n"
        f"            digits = {digits}L,\n"
        '            table.attr = \'data-quarto-disable-processing="true"\',\n'
        '            "html"\n'
        "          ) |\n"
        "          kable_classic(full_width = FALSE, html_font = \"Cambria\") |\n"
        "          row_spec(0L, italic = TRUE, bold = TRUE) |\n"
        "          column_spec(1L, italic = FALSE, bold = TRUE)\n"
        "      )\n"
        "    }"
    )
    assert old_block in src, f"Could not find kable branch with digits={digits} and context:\n{region_col_block}"
    return src.replace(old_block, new_block, 1)


# paired.ttest branch (digits=3, ends with # Always return named list comment)
PAIRED_CTX = "    if (!is.null(region_col)) {\n      results <- results |> mutate(!!region_col := region_name)\n    }"
src = replace_kable_branch(src, PAIRED_CTX, 3)

# unpaired.ttest branch (digits=3, same region_col pattern but different surroundings)
# After first replacement the paired one is gone; second occurrence is unpaired
src = replace_kable_branch(src, PAIRED_CTX, 3)

# deseq branch (digits=2, same region_col pattern - third occurrence)
src = replace_kable_branch(src, PAIRED_CTX, 2)

# voom branch: uses `results <- results |> relocate("Gene")` before kable
VOOM_CTX = '    results <- results |> relocate("Gene")'
old_voom = (
    VOOM_CTX +
    "\n"
    "\n"
    "    if (kable.out) {\n"
    "      return(\n"
    "        results |\n"
    "          slice_head(n = ngenes.out) |\n"
    "          kable(\n"
    "            digits = 2L,\n"
    '            table.attr = \'data-quarto-disable-processing="true"\',\n'
    '            "html"\n'
    "          ) |\n"
    "          kable_classic(full_width = FALSE, html_font = \"Cambria\") |\n"
    "          row_spec(0L, italic = TRUE, bold = TRUE) |\n"
    "          column_spec(1L, italic = FALSE, bold = TRUE)\n"
    "      )\n"
    "    }"
)
new_voom = (
    VOOM_CTX +
    "\n"
    "\n"
    "    if (!is.null(genes.filter)) {\n"
    "      results <- dplyr::filter(results, .data$Gene %in% genes.filter)\n"
    "    }\n"
    "\n"
    "    if (kable.out) {\n"
    "      tbl <- if (is.null(genes.filter)) slice_head(results, n = ngenes.out) else results\n"
    "      return(\n"
    "        tbl |\n"
    "          kable(\n"
    "            digits = 2L,\n"
    '            table.attr = \'data-quarto-disable-processing="true"\',\n'
    '            "html"\n'
    "          ) |\n"
    "          kable_classic(full_width = FALSE, html_font = \"Cambria\") |\n"
    "          row_spec(0L, italic = TRUE, bold = TRUE) |\n"
    "          column_spec(1L, italic = FALSE, bold = TRUE)\n"
    "      )\n"
    "    }"
)
assert old_voom in src, "Could not find voom kable branch"
src = src.replace(old_voom, new_voom, 1)

# edgeR branch (LRT/QLF) - top-level (no indentation), uses 2-space indent
EDGER_CTX = '  results <- results |> relocate("Gene")'
old_edger = (
    EDGER_CTX +
    "\n"
    "\n"
    "  if (kable.out) {\n"
    "    return(\n"
    "      results |\n"
    "        slice_head(n = ngenes.out) |\n"
    "        kable(\n"
    "          digits = 2L,\n"
    '          table.attr = \'data-quarto-disable-processing="true"\',\n'
    '          "html"\n'
    "        ) |\n"
    "        kable_classic(full_width = FALSE, html_font = \"Cambria\") |\n"
    "        row_spec(0L, italic = TRUE, bold = TRUE) |\n"
    "        column_spec(1L, italic = FALSE, bold = TRUE)\n"
    "    )\n"
    "  }\n"
    "\n"
    "  return(results)\n"
    "}"
)
new_edger = (
    EDGER_CTX +
    "\n"
    "\n"
    "  if (!is.null(genes.filter)) {\n"
    "    results <- dplyr::filter(results, .data$Gene %in% genes.filter)\n"
    "  }\n"
    "\n"
    "  if (kable.out) {\n"
    "    tbl <- if (is.null(genes.filter)) slice_head(results, n = ngenes.out) else results\n"
    "    return(\n"
    "      tbl |\n"
    "        kable(\n"
    "          digits = 2L,\n"
    '          table.attr = \'data-quarto-disable-processing="true"\',\n'
    '          "html"\n'
    "        ) |\n"
    "        kable_classic(full_width = FALSE, html_font = \"Cambria\") |\n"
    "        row_spec(0L, italic = TRUE, bold = TRUE) |\n"
    "        column_spec(1L, italic = FALSE, bold = TRUE)\n"
    "    )\n"
    "  }\n"
    "\n"
    "  return(results)\n"
    "}"
)
assert old_edger in src, "Could not find edgeR kable branch"
src = src.replace(old_edger, new_edger, 1)

with open(PATH, "w") as f:
    f.write(src)

print("All changes applied successfully.")
print(f"Total characters: {len(src)}")
print(f"Total lines: {src.count(chr(10))}")
