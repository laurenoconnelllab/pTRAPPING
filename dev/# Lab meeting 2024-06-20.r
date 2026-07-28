# Lab meeting 2024-06-20

library(poisonfrogs)

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
    prior.count = 0,
    genes.filter = c(
      "Adcyap1",
      "Bdnf",
      "Ucn3",
      "Gng8",
      "Fosl2",
      "Junb",
      "Trappc12",
      "Gfap"
    )
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
    prior.count = 0,
    genes.filter = c(
      "Adcyap1",
      "Bdnf",
      "Ucn3",
      "Gng8",
      "Fosl2",
      "Junb",
      "Trappc12",
      "Gfap"
    )
  )


bind_rows(pacap_de$results, bdnf_de$results) |>
  ggplot(aes(
    x = treatment,
    y = Gene,
    size = -log10(FDR),
    fill = logFC
  )) +
  geom_point(
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(2, 12)
  ) +
  scale_fill_gradientn(
    colours = c("#25599b", "#ABD0DD", "#F2F9FE", "#F88705", "#B8351F")
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Treatment",
    y = "Gene"
  )

ptrap_volcano()


bind_rows(pacap_de$results, bdnf_de$results) |>
  ptrap_bubble()


filterByExpr()
