library(ggrepel)
library(tidyverse)
library(readr)
library(edgeR)
library(stringr)
# load counts and metadata
counts_rimi <- read_tsv(
  "/Users/camilorl/Quarto_projects/PairBond_PhTRAP/counts_Rimitator_gene.txt",
  comment = "#"
)

frogs_metadata <- read_csv(
  "/Users/camilorl/Quarto_projects/PairBond_PhTRAP/frogs_metadata.csv"
)
# remove extra annotation columns to get counts matrix
counts <- counts_rimi %>%
  select(-Chr, -Start, -End, -Strand, -Length)

# store gene IDs separately for later use and remove from counts matrix
gene_ids <- counts$Geneid
counts <- counts %>% select(-Geneid)

# create sample metadata dataframe by parsing sample names and joining with frog metadata
samples <- tibble(sample = colnames(counts)) %>%
  mutate(
    tube = str_extract(sample, "\\d+") |> as.integer(),
    fraction = if_else(str_detect(sample, "IP"), "IP", "INPUT")
  ) |>
  left_join(frogs_metadata, by = "tube") |>
  # optional: standardize names
  rename(BrainRegion = `Brain region`)

# remove samples with known issues (e.g., low library size, outliers, etc.)
samples <- samples %>%
  filter(!tube %in% c(2, 17, 27))


# De analysis using pTRAPPING's wrapper function for DESeq2, with LFC shrinkage
rimim_pb <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "pb",
  block_col = "tube",
  test_method = "deseq",
  shrink.lfc = TRUE
)

rimim_pb |>
  ptrap_volcano(fdr = FALSE, genes.annot = c("TRH", "NPY", "AVP"))

samples |>
  print(n = Inf)

rimim_sol <- ptrap_de(
  counts_mat = counts,
  sample_df = samples,
  gene_ids = gene_ids,
  region_name = "POA",
  treatment_name = "sol",
  block_col = "tube",
  test_method = "deseq",
  shrink.lfc = TRUE
)

ptrap_volcano2(rimim_pb, rimim_sol)

rimim_pb |>
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

# DE analysis using edgeR's quasi-likelihood pipeline, with paired design and LFC shrinkage via glmTreat
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


str(rimim_pb_pttest)

rimim_pb_pttest$results
rimim_pb_pttest$fe


fff <- c(
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


class(fff)
