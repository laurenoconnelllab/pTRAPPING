# Lab meeting 2024-06-20

?ptrap_de()
?ptrap_volcano()

head(warm_counts)


warm_counts |>
  ptrap_de(
    norm.method = "none"
  ) |>
  ptrap_volcano(genes.annot = c("Adcyap1", "Bdnf"))


# Pre-computed RPKM values — already normalised; use with norm.method = "none"
counts_rpkm <- read.delim(
  system.file("extdata", "TAN_etal_2016_RPKM.txt", package = "pTRAPPING")
)

dim(counts_rpkm) # 8863 genes × 9 columns (1 gene-ID + 8 samples)
names(counts_rpkm) # column names encode treatment, fraction, and replicate


pacap_de <- counts_rpkm |>
  filter(dplyr::if_any(dplyr::where(is.numeric), ~ .x > 1)) |> #keep counts > 1
  mutate(gene = make.unique(Gene)) |>
  tibble::column_to_rownames("gene") |>
  #round() |>
  ptrap_de(
    test_method = "paired.ttest",
    norm.method = "none",
    treatment_name = "PACAP",
    filter = FALSE,
    lfc_threshold = 0.3,
    prior.count = 0
  )


bdnf_de <- counts_rpkm |>
  filter(dplyr::if_any(dplyr::where(is.numeric), ~ .x > 1)) |> #keep counts > 1
  mutate(gene = make.unique(Gene)) |>
  tibble::column_to_rownames("gene") |>
  #round() |>
  ptrap_de(
    test_method = "paired.ttest",
    norm.method = "none",
    treatment_name = "BDNF",
    filter = FALSE,
    lfc_threshold = 0.3,
    prior.count = 0
  )

ptrap_volcano2(
  de_result_1 = pacap_de$results,
  de_result_2 = bdnf_de$results,
  fdr = FALSE,
  title = "PACAP vs BDNF — Preoptic Area",
  interactive = TRUE
)
