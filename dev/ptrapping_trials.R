library(ggrepel)
library(tidyverse)
library(dplyr)
library(readr)
library(edgeR)
library(stringr)
# load counts and metadata
counts_rimi_0 <- read_tsv(
  "/Users/camilorodriguezlopez/Quarto_projects/pTRAPPING/dev/counts_Rimitator_gene.txt",
  comment = "#"
)

frogs_metadata <- read_csv(
  "/Users/camilorodriguezlopez/Quarto_projects/pTRAPPING/dev/frogs_metadata.csv"
)
# remove extra annotation columns to get counts matrix
counts_rimi_0 <- counts_rimi_0 %>%
  select(-Chr, -Start, -End, -Strand, -Length)

# store gene IDs separately for later use and remove from counts matrix
gene_ids <- counts_rimi_0$Geneid
counts_rimi <- counts_rimi_0 %>% select(-Geneid)

# create sample metadata dataframe by parsing sample names and joining with frog metadata
samples <- tibble(sample = colnames(counts)) |>
  mutate(
    tube = str_extract(sample, "\\d+") |> as.integer(),
    fraction = if_else(str_detect(sample, "IP"), "IP", "INPUT")
  ) |>
  left_join(frogs_metadata, by = "tube") |>
  # optional: standardize names
  rename(BrainRegion = `Brain region`)

# remove samples with known issues (e.g., low library size, outliers, etc.)
samples_rimi <- samples %>%
  filter(!tube %in% c(2, 17, 27))


# DE analysis using DESeq2 pipeline, with LFC shrinkage
rimim_pb_deseq <- ptrap_de(
  counts_mat = counts_rimi,
  sample_df = samples_rimi,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "pb",
  block_col = "tube",
  test_method = "deseq",
  shrink.lfc = TRUE
)

rimim_pb_deseq |>
  ptrap_volcano(fdr = TRUE, genes.annot = c("TRH", "NPY", "AVP"))

samples |>
  print(n = Inf)

rimim_sol_deseq <- ptrap_de(
  counts_mat = counts_rimi,
  sample_df = samples_rimi,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "sol",
  block_col = "tube",
  test_method = "deseq",
  shrink.lfc = TRUE
)

ptrap_volcano2(rimim_pb_deseq, rimim_sol_deseq, fdr = FALSE)

rimim_pb_deseq |>
  filter(PValue < 0.05) |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )


rimim_pb_deseq |>
  filter(PValue < 0.05) |>
  print(n = Inf) |>
  filter(Gene == "AANAT")


rimim_sol |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )

# DE analysis using paired t-test pipeline
rimim_pb_pttest <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "pb",
  block_col = "tube",
  test_method = "paired.ttest",
  norm.method = "mratios"
)

rimim_pb_pttest$results |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )


rimim_sol_pttest <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "sol",
  block_col = "tube",
  test_method = "paired.ttest",
  norm.method = "CPM"
)

rimim_sol_pttest$results |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )


ptrap_volcano2(
  rimim_pb_pttest$results,
  rimim_sol_pttest$results,
  fdr = FALSE,
  genes.annot = c("TRH", "NPY", "AVP")
)


# DE analysis using the voom + limma pipeline
rimim_pb_voom <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "pb",
  block_col = "tube",
  test_method = "voom"
)

ptrap_volcano(rimim_pb_voom, fdr = FALSE, genes.annot = c("TRH", "NPY", "AVP"))

rimim_pb_voom |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )

rimim_sol_voom <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "sol",
  block_col = "tube",
  test_method = "voom"
)

rimim_sol_voom |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )

ptrap_volcano2(
  rimim_pb_voom,
  rimim_sol_voom,
  fdr = FALSE,
  genes.annot = c("TRH", "NPY", "AVP")
)


# DE analysis using the edgeR pipeline - QLF test
rimim_pb_QLF <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "pb",
  block_col = "tube",
  test_method = "QLF"
)

rimim_pb_QLF |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )

rimim_sol_QLF <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "sol",
  block_col = "tube",
  test_method = "QLF"
)

rimim_sol_QLF |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )

ptrap_volcano2(
  rimim_pb_QLF,
  rimim_sol_QLF,
  fdr = TRUE,
  genes.annot = c("TRH", "NPY", "AVP")
)

# DE analysis using the unpaired t-test pipeline
rimim_unp_ttest <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "pb",
  block_col = "tube",
  test_method = "unpaired.ttest",
  norm.method = "CPM"
)

ptrap_volcano(rimim_unp_ttest, fdr = TRUE, genes.annot = c("TRH", "NPY", "AVP"))


#R. variabilis data

counts_rv_0 <- read_tsv(
  "/Users/camilorodriguezlopez/Quarto_projects/pTRAPPING/dev/counts_Rvariabilis_gene.txt",
  comment = "#"
)

frogs_metadata <- read_csv(
  "/Users/camilorodriguezlopez/Quarto_projects/pTRAPPING/dev/frogs_metadata.csv"
)
counts_rv_0 <- counts_rv_0 %>%
  select(-Chr, -Start, -End, -Strand, -Length)

gene_ids_rv <- counts_rv_0$Geneid
counts_rv <- counts_rv_0 %>% select(-Geneid)

samples_rv <- tibble(sample = colnames(counts_rv)) %>%
  mutate(
    tube = str_extract(sample, "\\d+") |> as.integer(),
    fraction = if_else(str_detect(sample, "IP"), "IP", "INPUT")
  ) |>
  left_join(frogs_metadata, by = "tube") |>
  rename(BrainRegion = `Brain region`)


rvar_sol_pttest <- ptrap_de(
  counts_mat = counts_rv,
  sample_df = samples_rv,
  gene_ids = gene_ids_rv,
  region_name = "POA",
  treatment_name = "sol",
  block_col = "tube",
  test_method = "paired.ttest",
  norm.method = "mratios"
)

ptrap_volcano2(
  rimim_pb_pttest$results,
  rvar_sol_pttest$results,
  fdr = FALSE,
  genes.annot = c("TRH", "NPY")
)

rimim_pb_pttest$results |>
  filter(PValue < 0.05) |>
  #filter(diffexpressed == "UP") |>
  print(n = Inf) |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )

rvar_sol_deseq |>
  filter(PValue < 0.05) |>
  print(n = Inf) |>
  filter(
    Gene %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )


gene_ids_rv |>
  data.frame() |>
  filter(
    gene_ids_rv %in%
      c(
        "VIP",
        "TRH",
        "PRRP",
        "PNOC",
        "PENKA",
        "PDYN",
        "OXT",
        "NPY",
        "GAL",
        "FAAH1",
        "DRD5",
        "CRFR1",
        "CART",
        "AVP",
        "HTR1E"
      )
  )


rvar_poa <- counts_rv_0 |>
  column_to_rownames("Geneid") |>
  select(
    samples_rv |>
      filter(BrainRegion == "POA") |>
      select(sample) |>
      pull()
  )

colnames(rvar_poa) <- samples_rv |>
  filter(BrainRegion == "POA") |>
  mutate(sample2 = paste(Treatment, tube, fraction, sep = "_")) |>
  select(sample2) |>
  pull()


# Counts matrix: column names can be anything
counts <- data.frame(
  Gene = c("Actb", "Gapdh", "Tuba1a"),
  s1_a = c(120, 340, 210),
  s1_b = c(980, 210, 870),
  s2_a = c(130, 360, 190),
  s2_b = c(910, 230, 800),
  s3_a = c(115, 320, 220),
  s3_b = c(950, 200, 850)
)

# sample_df: explicitly maps each column to its metadata
sample_df <- data.frame(
  sample = c("s1_a", "s1_b", "s2_a", "s2_b", "s3_a", "s3_b"),
  Treatment = rep("WT", 6),
  block = c("1", "1", "2", "2", "3", "3"),
  fraction = c("INPUT", "IP", "INPUT", "IP", "INPUT", "IP")
)

ptrap_de(
  counts_mat = counts,
  sample_df = sample_df,
  sample_col = "sample",
  treatment_col = "Treatment",
  block_col = "block",
  fraction_col = "fraction"
)
