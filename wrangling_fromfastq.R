library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)

col_spec <- str_c("??", (str_dup('_?', 95)))
col_names <- c('gene_id', paste0('D25-21000', 1:9), paste0('D25-2100', 10:96))
seancounts <- read_delim("data/raw/allcountz_2.txt", col_types = col_spec, 
                         show_col_types = F, col_names = col_names) |>
  separate_wider_delim('gene_id', delim = '.', names = c('ens_id', 'version'), too_few = 'align_start') |>
  dplyr::select(!version) 


# remove weird ensembl genes <3

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes")
ensembl_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Get Ensembl gene IDs for all protein-coding genes
# We use the 'biotype' filter with the value 'protein_coding'
protein_coding_genes <- getBM(
  attributes = c("ensembl_gene_id", 'hgnc_symbol'),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl_dataset
) |>
  filter(hgnc_symbol != "")

symbols <- mapIds(org.Hs.eg.db, keys = seancounts$ens_id,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]

filtered <- protein_coding_genes |>
  left_join(y = seancounts, by = join_by('ensembl_gene_id' == 'ens_id'))


mapper <- read_csv("data/mapper.csv") |>
  dplyr::select(!uhh_idk) # idk what that column is saying tbh

data <- filtered |>
  dplyr::select(hgnc_symbol:`D25-210070`) |> # IGNORING TECH REPS: 71-96
  # fun nuance here is that some genes have multiple ensembl IDs. gotta love bioinformatics
  group_by(hgnc_symbol) |>
  summarise(across(everything(), sum)) |>
  pivot_longer(cols = starts_with("D25-"), names_to = 'core_id', values_to = 'counts') |>
  left_join(y = mapper) |>
  # IF WE WANT TO MERGE TECH REPS:
  #group_by(timepoint, donor, drug, symbols) |>
  # create new column: "tech rep average counts: but actually im summing <3 lol"
  # now core_id isn't meaningful </3
  #summarize(counts = sum(counts)) |>
  #ungroup() |>
  mutate(tdd_id = str_c(timepoint, donor, drug, sep = "_"),
         drug_broad = case_when(
          drug == 'dmso' ~ 'dmso',
          drug == 'ctrl' ~ 'ctrl',
          TRUE ~ 'jnki'
         ))


