library(tidyverse)

only96 <- read_delim("data/D25-210096-7213E_geneexp.txt", 
                     col_names = c("ensembl_id", 'counts', 'gene')) |>
  filter(!is.na(gene))

mapper <- read_csv("data/mapper.csv") |>
  dplyr::select(!uhh_idk) # idk what that column is saying tbh

count_matrix <- read_delim("data/geneexp.merge.txt") |>
  mutate(across(everything(), as.double)) |>
  bind_cols(only96 |> dplyr::select(!counts)) |>
  column_to_rownames('ensembl_id') |>
  dplyr::select(!gene)

colSums(count_matrix) |> hist()


data <- count_matrix |>
  bind_cols(only96 |> select(!counts)) |>
  pivot_longer(cols = starts_with("D25-"), names_to = 'core_id', values_to = 'counts') |>
  left_join(y = mapper) |>
# want to merge technical replicates, those aren't meaningfully different. i hope.
  group_by(timepoint, donor, drug, ensembl_id) |>
  # create new column: "tech rep average counts"
  # now core_id isn't meaningful </3
  # also round to nearest for DESEQ purposes
  summarize(tra_counts = round(mean(counts))) |>
  ungroup() |>
  mutate(tdd_id = str_c(timepoint, donor, drug, sep = "_"),
         drug_broad = case_when(
          drug == 'dmso' ~ 'dmso',
          drug == 'ctrl' ~ 'ctrl',
          TRUE ~ 'jnki'
         ))



write_csv(data, file = 'data/initial_counts_long.csv')

count_matrix_new <- data |>
  select(ensembl_id, tra_counts, tdd_id) |>
  pivot_wider(values_from = tra_counts, names_from = tdd_id) |>
  column_to_rownames('ensembl_id')



## get DESEQ variance stabilized data

library(DESeq2)


# Filter out lowly expressed genes
# Retain those with at least 1 CPM in 10 samples
cpm_counts=edgeR::cpm(count_matrix_new)
cpm_more1<-apply(cpm_counts, 1, function(x) sum(x>1))
idx_keep=which(cpm_more1>=10)
count_matrix_new_new <- count_matrix_new[idx_keep,]


# prep coldata
timepointt <- '6h'

coldata <- data |>
  select(tdd_id, drug_broad, timepoint, donor, drug) |>
  unique() |>
  mutate(across(everything(), as.factor)) |>
  filter(timepoint == timepointt) |>
  select(!timepoint)

rownames(coldata) <- coldata$tdd_id


# Normalize the filtered data using DESeq2's VST
dds <- DESeqDataSetFromMatrix(countData = count_matrix_new_new |> select(starts_with('6h')),
                                  colData = coldata,
                                  design = as.formula('~ drug_broad'))

dds <- DESeq(dds)
vsd.fixed <- varianceStabilizingTransformation(dds, blind=TRUE)
vst_expr <- assay(vsd.fixed)

vst_data <- vst_expr |>
  as.data.frame() |>
  rownames_to_column('ensembl_id') |>
  left_join(only96 |> select(!counts)) |>
  mutate(gene = case_when(
    # these are duplicates
    gene %in% c('CRYBG3', 'IDS', 'POLR2J4', 'ZNF33B') ~ ensembl_id,
    TRUE ~ gene
  )) |>
  column_to_rownames('gene') |>
  select(!ensembl_id) |>
  t() |>
  as.data.frame() |>
  rownames_to_column('tdd_id') |>
  separate_wider_delim(cols = 'tdd_id', delim = '_', names = c('timepoint', 'donor', 'drug')) |>
  mutate(drug_broad = case_when(
          drug == 'dmso' ~ 'dmso',
          drug == 'ctrl' ~ 'ctrl',
          TRUE ~ 'jnki'
         ), .after = drug) 


write_csv(vst_data, file = 'data/processed/deseq_varstabilized_data.csv')
