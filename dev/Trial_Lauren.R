?ptrap_de()


#Data sets

#all data per region, all treatments

#R. imitator
counts_rimi_0
#R. variabilis
counts_rv_0


#one region - POA

#R. imitator
rimi_poa
#R. variabilis
rvar_poa


#metadata
samples_rimi
samples_rv


rimi_de <- ptrap_de(rimi_poa, treatment_name = "pb", test_method = "voom")
rvar_de <- ptrap_de(rvar_poa, treatment_name = "sol", test_method = "voom")

ptrap_volcano2(rimi_de, rvar_de, fdr = FALSE, genes.annot = c("TRH", "NPY"))


bolt_counts <- read_tsv(
  "/Users/camilorodriguezlopez/Library/CloudStorage/Dropbox/LOBSU_postdoc/Figures4Lauren/bolt_counts_matrix.txt"
)

food_counts <- read_tsv(
  "/Users/camilorodriguezlopez/Library/CloudStorage/Dropbox/LOBSU_postdoc/Figures4Lauren/food_counts_matrix.txt"
)

female_counts <- read_tsv(
  "/Users/camilorodriguezlopez/Library/CloudStorage/Dropbox/LOBSU_postdoc/Figures4Lauren/female_counts_matrix.txt"
)

tad_counts <- read_tsv(
  "/Users/camilorodriguezlopez/Library/CloudStorage/Dropbox/LOBSU_postdoc/Figures4Lauren/tad_counts_matrix.txt"
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
    Gene %in% c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP")
  ) |>

  food_de |>
  filter(Gene %in% c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP"))

female_de |>
  filter(Gene %in% c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP"))

tad_de |>
  filter(Gene %in% c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP"))

bolt_food_plot <- ptrap_volcano2(
  bolt_de,
  food_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP")
)


female_food_plot <- ptrap_volcano2(
  female_de,
  food_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP")
)

female_bolt_plot <- ptrap_volcano2(
  female_de,
  bolt_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP")
)

tad_food_plot <- ptrap_volcano2(
  tad_de,
  food_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP")
)

tad_bolt_plot <- ptrap_volcano2(
  tad_de,
  bolt_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP")
)

tad_female_plot <- ptrap_volcano2(
  tad_de,
  female_de,
  fdr = TRUE,
  point_alpha = 0.1,
  treatment_col = "treatment",
  genes.annot = c("Ucn", "Ucn3", "crhbp", "crhr1", "crhr2", "HCRT", "AGRP")
)
