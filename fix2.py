import re
PATH = "R/ptrap_de.R"
with open(PATH, "r") as f:
    src = f.read()

# change 1: parser
c1_old = ("  # Step 2: Inject a canonical separator (\"_\") just before each fraction\n"
  "  # keyword, so they are always split off as their own token regardless of\n"
  "  # whether they are fused with the treatment name.\n"
  "  # INPUT must be handled before IP so the longer keyword wins.\n"
  "  # Use case-insensitive replacement via perl = TRUE + (?i).\n"
  "  cleaned <- gsub(\n"
  "    paste0(\"(?i)(?=\", input_level, \")\"),\n"
  "    \"_\",\n"
  "    cleaned,\n"
  "    perl = TRUE\n"
  "  )\n"
  "  cleaned <- gsub(\n"
  "    paste0(\"(?i)(?=\", ip_level, \")\"),\n"
  "    \"_\",\n"
  "    cleaned,\n"
  "    perl = TRUE\n"
  "  )\n"
  "\n"
  "  # Step 3: Re-strip any double separators produced above.")
c1_new = ("  # Step 2: Surround fraction keyword with \"_\" on both sides.\n"
  "  # input_level listed first so PCRE picks the longer keyword (prevents\n"
  "  # \"IP\" matching inside \"INPUT\").\n"
  "  cleaned <- gsub(\n"
  "    paste0(\"(?i)(\", input_level, \"|\", ip_level, \")\"),\n"
  "    \"_\\\\1_\",\n"
  "    cleaned,\n"
  "    perl = TRUE\n"
  "  )\n"
  "  # Handle \"in\" alias (bare, not inside another word).\n"
  "  cleaned <- gsub(\n"
  "    \"(?i)(?<![A-Za-z])(in)(?![A-Za-z])\",\n"
  "    \"_\\\\1_\",\n"
  "    cleaned,\n"
  "    perl = TRUE\n"
  "  )\n"
  "\n"
  "  # Step 3: Re-strip any double separators produced above.")
assert c1_old in src, "c1_old not found"
src = src.replace(c1_old, c1_new, 1)
print("1 parser OK")

# change 2: header
c2_old = ("# Examples (all produce the same output fields):\n"
  "#   \"b1input\"   -> treatment=\"b\",  block=\"1\", fraction=input_level\n"
  "#   \"nb1ip\"     -> treatment=\"nb\", block=\"1\", fraction=ip_level\n"
  "#   \"Nb_IP_1\"   -> treatment=\"Nb\", block=\"1\", fraction=ip_level\n"
  "#   \"B.3.INPUT\" -> treatment=\"B\",  block=\"3\", fraction=input_level")
c2_new = ("# All orderings of T / R / F are supported:\n"
  "#   T-R-F: \"PACAP1IP\"    T-F-R: \"PACAP_Input_1\"  R-T-F: \"1_PACAP_IP\"\n"
  "#   R-F-T: \"1_IP_PACAP\"  F-T-R: \"IP_PACAP_1\"     F-R-T: \"IP_1_PACAP\"\n"
  "# \"in\" (case-insensitive) is accepted as a short alias for input_level.")
assert c2_old in src, "c2_old not found"
src = src.replace(c2_old, c2_new, 1)
print("2 header OK")

# change 3: roxygen
c3_old = ("#' @param ngenes.out Number of top genes (sorted by p-value) to include in the\n"
  "#'   output when `kable.out = TRUE`. Default is `20`.\n"
  "#' @param kable.out")
c3_new = ("#' @param ngenes.out Number of top genes (sorted by p-value) to include in the\n"
  "#'   output when `kable.out = TRUE`. Ignored when `genes.filter` is supplied.\n"
  "#'   Default is `20`.\n"
  "#' @param genes.filter Optional character vector of gene names (matching the\n"
  "#'   `Gene` column of the output) to retain. When supplied, only the listed\n"
  "#'   genes are returned and `ngenes.out` is ignored. Default is `NULL`.\n"
  "#' @param kable.out")
assert c3_old in src, "c3_old not found"
src = src.replace(c3_old, c3_new, 1)
print("3 roxygen OK")

# change 4: signature
c4_old = "  ngenes.out = 20,\n  kable.out = FALSE,"
c4_new = "  ngenes.out = 20,\n  genes.filter = NULL,\n  kable.out = FALSE,"
assert c4_old in src, "c4_old not found"
src = src.replace(c4_old, c4_new, 1)
print("4 signature OK")

# change 5: patch kable branches
KABLE_PAT = re.compile(
    r'(?P<indent>[ ]*)if \(kable\.out\) \{\n'
    r'(?P=indent)  return\(\n'
    r'(?P=indent)    results \|>\n'
    r'(?P=indent)      slice_head\(n = ngenes\.out\) \|>\n'
    r'(?P=indent)      kable\(\n'
    r'(?P=indent)        digits = (?P<digits>\d+)L,\n'
    r'(?P=indent)        table\.attr = \'data-quarto-disable-processing="true"\',\n'
    r'(?P=indent)        "html"\n'
    r'(?P=indent)      \) \|>\n'
    r'(?P=indent)      kable_classic\(full_width = FALSE, html_font = "Cambria"\) \|>\n'
    r'(?P=indent)      row_spec\(0L, italic = TRUE, bold = TRUE\) \|>\n'
    r'(?P=indent)      column_spec\(1L, italic = FALSE, bold = TRUE\)\n'
    r'(?P=indent)  \)\n'
    r'(?P=indent)\}'
)

def replacement(m):
    ind = m.group("indent")
    d = m.group("digits")
    return (
        f"{ind}if (!is.null(genes.filter)) {{\n"
        f"{ind}  results <- dplyr::filter(results, .data$Gene %in% genes.filter)\n"
        f"{ind}}}\n\n"
        f"{ind}if (kable.out) {{\n"
        f"{ind}  tbl <- if (is.null(genes.filter)) slice_head(results, n = ngenes.out) else results\n"
        f"{ind}  return(\n"
        f"{ind}    tbl |>\n"
        f"{ind}      kable(\n"
        f"{ind}        digits = {d}L,\n"
        f"{ind}        table.attr = 'data-quarto-disable-processing=\"true\"',\n"
        f'{ind}        "html"\n'
        f"{ind}      ) |>\n"
        f'{ind}      kable_classic(full_width = FALSE, html_font = "Cambria") |>\n'
        f"{ind}      row_spec(0L, italic = TRUE, bold = TRUE) |>\n"
        f"{ind}      column_spec(1L, italic = FALSE, bold = TRUE)\n"
        f"{ind}  )\n"
        f"{ind}}}"
    )

new_src, count = KABLE_PAT.subn(replacement, src)
assert count == 5, f"Expected 5, got {count}"
src = new_src
print(f"5 kable OK (replaced {count})")

with open(PATH, "w") as f:
    f.write(src)
print(f"Done. Lines: {src.count(chr(10))}")
