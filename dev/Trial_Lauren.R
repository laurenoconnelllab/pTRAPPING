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


bolt_de |>
  filter(
    Gene %in% c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal")
  )

food_de |>
  filter(Gene %in% c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal"))

female_de |>
  filter(Gene %in% c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal"))

tad_de |>
  filter(Gene %in% c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal"))

bolt_food_plot <- ptrap_volcano2(
  bolt_de,
  food_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal")
)


female_food_plot <- ptrap_volcano2(
  female_de,
  food_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("npy", "pomc", "CARTPT", "AGRP", "HCRT", "gal")
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
