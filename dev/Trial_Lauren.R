library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(ggpubr)


#Data sets

bolt_counts <- read_tsv(
  "/Users/camilorl/Library/CloudStorage/Dropbox/LOBSU_postdoc/Figures4Lauren/bolt_counts_matrix.txt"
)

food_counts <- read_tsv(
  "/Users/camilorl/Library/CloudStorage/Dropbox/LOBSU_postdoc/Figures4Lauren/food_counts_matrix.txt"
)

female_counts <- read_tsv(
  "/Users/camilorl/Library/CloudStorage/Dropbox/LOBSU_postdoc/Figures4Lauren/female_counts_matrix.txt"
)

tad_counts <- read_tsv(
  "/Users/camilorl/Library/CloudStorage/Dropbox/LOBSU_postdoc/Figures4Lauren/tad_counts_matrix.txt"
)


bolt_de <- bolt_counts |>
  rename_with(
    ~ sub(
      "^(bolt)_(input|ip)_([0-9]+)$",
      "\\1_\\3_\\2",
      .x,
      ignore.case = TRUE
    ),
    .cols = -gene
  ) |>
  column_to_rownames("gene") |>
  round() |>
  ptrap_de(
    test_method = "deseq"
  )

bolt_de_tt <- bolt_counts |>
  rename_with(
    ~ sub(
      "^(bolt)_(input|ip)_([0-9]+)$",
      "\\1_\\3_\\2",
      .x,
      ignore.case = TRUE
    ),
    .cols = -gene
  ) |>
  column_to_rownames("gene") |>
  round() |>
  ptrap_de(
    test_method = "paired.ttest",
    norm.method = "none",
    filter = FALSE,
    lfc_threshold = 0.3,
    prior.count = 0.1
  )


food_de <- food_counts |>
  rename_with(
    ~ sub(
      "^(food)_(input|ip)_([0-9]+)$",
      "\\1_\\3_\\2",
      .x,
      ignore.case = TRUE
    ),
    .cols = -gene
  ) |>
  column_to_rownames("gene") |>
  round() |>
  ptrap_de(
    test_method = "deseq"
  )

food_de_tt <- food_counts |>
  rename_with(
    ~ sub(
      "^(food)_(input|ip)_([0-9]+)$",
      "\\1_\\3_\\2",
      .x,
      ignore.case = TRUE
    ),
    .cols = -gene
  ) |>
  column_to_rownames("gene") |>
  round() |>
  ptrap_de(
    test_method = "paired.ttest",
    norm.method = "none",
    filter = TRUE,
    lfc_threshold = 0.3,
    prior.count = 0
  )


female_de <- female_counts |>
  rename_with(
    ~ sub(
      "^(female)_(input|ip)_([0-9]+)$",
      "\\1_\\3_\\2",
      .x,
      ignore.case = TRUE
    ),
    .cols = -gene
  ) |>
  column_to_rownames("gene") |>
  round() |>
  ptrap_de(
    test_method = "deseq"
  )

female_de_tt <- female_counts |>
  rename_with(
    ~ sub(
      "^(female)_(input|ip)_([0-9]+)$",
      "\\1_\\3_\\2",
      .x,
      ignore.case = TRUE
    ),
    .cols = -gene
  ) |>
  column_to_rownames("gene") |>
  round() |>
  ptrap_de(
    test_method = "paired.ttest",
    norm.method = "none",
    filter = FALSE,
    lfc_threshold = 0.3,
    prior.count = 0
  )

tad_de <- tad_counts |>
  rename_with(
    ~ sub(
      "^(tad)_(input|ip)_([0-9]+)$",
      "\\1_\\3_\\2",
      .x,
      ignore.case = TRUE
    ),
    .cols = -gene
  ) |>
  column_to_rownames("gene") |>
  round() |>
  ptrap_de(
    test_method = "deseq"
  )


bolt_de_tt$results |>
  filter(
    Gene %in%
      c(
        "npy",
        "pomc",
        "CARTPT",
        "AGRP",
        "HCRT",
        "gal",
        "fmr1-b",
        "fmr1-a",
        "fmr1"
      )
  )

food_de |>
  filter(
    Gene %in%
      c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal", "fmr1-b", "fmr1-a")
  )

female_de_tt$results |>
  filter(
    Gene %in%
      c(
        "npy",
        "pomc",
        "CARTPT",
        "AGRP",
        "HCRT",
        "gal",
        "fmr1-b",
        "fmr1-a",
        "fmr1"
      )
  )

tad_de |>
  filter(
    Gene %in%
      c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal", "fmr1-b", "fmr1-a")
  )

bolt_food_plot <- ptrap_volcano2(
  bolt_de,
  food_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c(
    "npy",
    "pomc",
    "CARTPT",
    "AGRP",
    "HCRT",
    "gal",
    "fmr1-b",
    "fmr1-a"
  )
)


female_food_plot <- ptrap_volcano2(
  female_de_tt$results,
  food_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("npy", "pomc", "CARTPT", "HCRT", "gal", "fmr1-b", "fmr1-a")
)

female_bolt_plot <- ptrap_volcano2(
  female_de,
  bolt_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal")
)

tad_food_plot <- ptrap_volcano2(
  tad_de,
  food_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal")
)

tad_bolt_plot <- ptrap_volcano2(
  tad_de,
  bolt_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal")
)

tad_female_plot <- ptrap_volcano2(
  tad_de,
  female_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal")
)


ggarrange(
  bolt_food_plot,
  female_food_plot,
  female_bolt_plot,
  tad_food_plot,
  tad_bolt_plot,
  tad_female_plot,
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "top"
)

### pTRAPPING demonstration

## Data from Tan et al., RPKM transformed counts
tan_rpkm <- read_tsv(
  "dev/TAN_etal_2016_RPKM.txt"
)

## Data from Tan et al, raw (no-transformed) counts
tan_raw <- read_tsv("dev/TAN_etal_2016_raw.txt")

#The simpliest usage of the functions

###_________DE genes____________###

tan_raw |>
  ptrap_de(treatment_name = "PACAP")

###_________DE genes in kable output____________###
tan_raw |>
  ptrap_de(
    treatment_name = "PACAP",
    kable.out = T,
    genes.filter = c(
      "Ucn3",
      "Fosl2",
      "Bdnf",
      "Adcyap1",
      "Junb",
      "Gadd45b"
    )
  )

###_________volcano plot for one condition____________###
tan_raw |>
  ptrap_de(treatment_name = "PACAP") |>
  ptrap_volcano()

###_________volcano plot for one condition with annotated genes____________###
tan_raw |>
  ptrap_de(treatment_name = "PACAP") |>
  ptrap_volcano(
    genes.annot = c(
      "Ucn3",
      "Fosl2",
      "Bdnf",
      "Adcyap1",
      "Junb",
      "Gadd45b"
    )
  )

###______Dual volcano plot comparing 2 conditions_____###

tan_dbnf <- tan_raw |>
  ptrap_de(treatment_name = "BDNF")

tan_pacap <- tan_raw |>
  ptrap_de(treatment_name = "PACAP")

ptrap_volcano2(tan_dbnf, tan_pacap)


###______Dual volcano plot with annotated genes_____###
ptrap_volcano2(
  tan_dbnf,
  tan_pacap,
  genes.annot = c(
    "Ucn3",
    "Fosl2",
    "Bdnf",
    "Adcyap1",
    "Junb",
    "Gadd45b"
  )
)

# ptrap_de():

## Show parsing sample metadata from column names

### --------------parsing 1------------------ ###

parsing1 <- read_tsv("dev/parser_examples/parsing1.txt")

parsing1 |>
  colnames()

parsing1 |>
  ptrap_de(treatment_name = "ABC")

### --------------parsing 2------------------- ###

parsing2 <- read_tsv("dev/parser_examples/parsing2.txt")

parsing2 |>
  colnames()

parsing2 |>
  ptrap_de(treatment_name = "cont")

### --------------parsing 3------------------- ###

parsing3 <- read_tsv("dev/parser_examples/parsing3.txt")

parsing3 |>
  ptrap_de(treatment_name = "a")


### --------------parsing 4: including samples from more than one brain region------------------- ###

parsing4 <- read_tsv("dev/parser_examples/parsing4.txt")

parsing4 |>
  ptrap_de(treatment_name = "a", region_name = "STR")


## test method

###__________edgeR pipeline: LRT__________###

tan_raw |>
  ptrap_de(treatment_name = "PACAP", test_method = "LRT") |>
  ptrap_volcano(genes.annot = c("Pwwp2a"))

###__________edgeR pipeline: QLF__________###

tan_raw |>
  ptrap_de(treatment_name = "PACAP", test_method = "QLF") |>
  ptrap_volcano(genes.annot = c("Pwwp2a"))

###__________ Deseq2 pipeline_________##
tan_raw |>
  ptrap_de(
    treatment_name = "PACAP",
    test_method = "deseq",
    shrink.lfc = TRUE
  ) |>
  ptrap_volcano(genes.annot = c("Pwwp2a"))


###__________limma-voom pipeline_________##
tan_raw |>
  ptrap_de(treatment_name = "PACAP", test_method = "voom") |>
  ptrap_volcano(genes.annot = c("Pwwp2a"))

###__________paired t-test_____________###

##Paired t-test pipeline with no transformation, no filtering
pacap_de <- tan_rpkm |>
  filter(if_any(where(is.numeric), ~ .x > 1)) |> #keep counts > 1
  mutate(gene = make.unique(Gene)) |>
  column_to_rownames("gene") |>
  #round() |>
  ptrap_de(
    test_method = "paired.ttest",
    norm.method = "none",
    treatment_name = "PACAP",
    filter = FALSE,
    lfc_threshold = 0.3,
    prior.count = 0
  )

pacap_de$fe |>
  mutate(log2FE_1 = log2(FE_1), log2FE_2 = log2(FE_2))

pacap_de$results |>
  ptrap_volcano(
    fdr = FALSE,
    gene_col = "Gene",
    treatment_col = "treatment",
    genes.annot = c(
      "Adcyap1",
      "Ghrh",
      "Ucn3",
      "Nxph4",
      "Bdnf"
    )
  )


pacap_de$results |>
  filter(
    Gene %in%
      c(
        "Adcyap1",
        "Ghrh",
        "Ucn3",
        "Nxph4",
        "Bdnf"
      )
  )


bdnf_de <- tan_rpkm |>
  filter(if_any(where(is.numeric), ~ .x > 1)) |>
  mutate(gene = make.unique(Gene)) |>
  column_to_rownames("gene") |>
  #round() |>
  ptrap_de(
    test_method = "paired.ttest",
    norm.method = "none",
    treatment_name = "BDNF",
    filter = FALSE,
    lfc_threshold = 0.3,
    prior.count = 0
  )

bdnf_de$fe |>
  mutate(log2FE_1 = log2(FE_1), log2FE_2 = log2(FE_2))

bdnf_de$results |>
  filter(
    Gene %in%
      c(
        "Adcyap1",
        "Ghrh",
        "Ucn3",
        "Nxph4",
        "Bdnf"
      )
  )


ptrap_volcano2(
  bdnf_de$results,
  pacap_de$results,
  fdr = FALSE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("Adcyap1", "Bdnf", "Ucn3", "Nxph4", "Ghrh")
)

warm_de <- read.delim(
  "/Users/camilorodriguezlopez/Downloads/GSE80121_PhosphoTRAP_Processed_Data.txt",
  skip = 2, # skip the two descriptor rows
  header = TRUE,
  check.names = FALSE
) |>
  rename(Gene = `Gene Name`) |>
  setNames(c(
    "Gene",
    "WC1_INPUT",
    "WC1_IP",
    "WC2_INPUT",
    "WC2_IP",
    "WC3_INPUT",
    "WC3_IP"
  )) |>
  ptrap_de(
    sample_df = tibble::tibble(
      sample = c(
        "WC1_INPUT",
        "WC1_IP",
        "WC2_INPUT",
        "WC2_IP",
        "WC3_INPUT",
        "WC3_IP"
      ),
      fraction = c("INPUT", "IP", "INPUT", "IP", "INPUT", "IP"),
      block = c("1", "1", "2", "2", "3", "3"),
      treatment = "WarmChallenge"
    ),
    sample_col = "sample",
    fraction_col = "fraction",
    block_col = "block",
    treatment_col = "treatment",
    treatment_name = "WarmChallenge",
    test_method = "paired.ttest",
    norm.method = "none",
    filter = FALSE,
    prior.count = 0,
    lfc_threshold = 0.3
  )


warm_de$results |>
  ptrap_volcano(
    fdr = FALSE,
    gene_col = "Gene",
    treatment_col = "treatment",
    genes.annot = c(
      "Fosl2",
      "Egr2",
      "Bdnf",
      "Adcyap1",
      "Junb",
      "Fosb",
      "Rrad",
      "Gadd45b"
    ),
    log_base = 2
  )

left_join(
  pacap_de$results |> select(Gene, logFC) |> rename(logFC_PACAP = logFC),
  bdnf_de$results |> select(Gene, logFC) |> rename(logFC_BDNF = logFC),
  by = "Gene"
) |>
  mutate(
    highlight = Gene %in% c("Adcyap1", "Ucn3", "Ghrh", "Nxph4", "Bdnf"),
    label = if_else(highlight, Gene, NA_character_)
  ) |>
  ggplot(aes(x = logFC_PACAP, y = logFC_BDNF)) +
  geom_point(aes(colour = highlight), alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  ggrepel::geom_label_repel(
    aes(label = label),
    colour = "#2166ac",
    size = 3,
    na.rm = TRUE,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_colour_manual(
    values = c("FALSE" = "grey60", "TRUE" = "#2166ac"),
    guide = "none"
  ) +
  labs(
    x = expression(log[2] ~ "PACAP IP/IN"),
    y = expression(log[2] ~ "BDNF IP/IN")
  ) +
  theme_classic()


ptrap_volcano2(
  de_result_1 = bdnf_de$results,
  de_result_2 = pacap_de$results,
  fdr = FALSE, # no FDR filtering, just show enrichment
  lfc_threshold = 0,
  gene_col = "Gene",
  treatment_col = "treatment",
  genes.annot = c("Adcyap1", "Ucn3", "Ghrh", "Nxph4", "Bdnf")
)


pacap_fe <- pacap_de$fe # columns: Gene, FE_1, FE_2
bdnf_fe <- bdnf_de$fe

pacap_fe |>
  left_join(bdnf_fe, by = "Gene", suffix = c("_pacap", "_bdnf")) |>
  filter(
    log2(FE_1_pacap) > 1, # PACAP rep 1
    log2(FE_2_pacap) > 1, # PACAP rep 2
    log2(FE_1_bdnf) > 1, # BDNF rep 1
    log2(FE_2_bdnf) > 1 # BDNF rep 2
  ) |>
  mutate(
    log2_mean_FE_pacap = log2((FE_1_pacap + FE_2_pacap) / 2),
    log2_mean_FE_bdnf = log2((FE_1_bdnf + FE_2_bdnf) / 2)
  ) |>
  select(Gene, log2_mean_FE_pacap, log2_mean_FE_bdnf) |>
  arrange(desc(log2_mean_FE_pacap))


##

write.table(
  warm_de,
  file = "dev/warm_counts.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
