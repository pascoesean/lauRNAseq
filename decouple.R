# test decouplR to see if transcription factor inference results make sense

# should probably start with 7 day data
source('functions.R')
source('wrangling_fromfastq.R')
library(decoupleR)

net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)


#setup variables
contrast_vector <- c('drug', '3013', 'dmso') # control induces x log fold change over dmso
timepointt <- '7d'
deseq_formula <- '~ drug + donor'
filename <- str_c(contrast_vector[2], 'vs', contrast_vector[3], timepointt, sep='_')

# setup subsets of data:
data7 <- data |> filter(timepoint == '7d')

count_matrix_new_new <- prep_data_(data7)
dds <- form_dds_obj_(count_matrix_new_new, data7, deseq_formula)

res <- results(dds, contrast=contrast_vector, independentFiltering = TRUE)

stats <- res@listData$stat |> as.matrix()
rownames(stats) <- res@rownames
# Run ulm
contrast_acts <- decoupleR::run_ulm(mat = stats, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor='mor', 
                                    minsize = 5)


# screaming crying throwing up. i am so happy rn
# note: does not work at 6 hour timepoint

# ETS family induced by jnki https://www.nature.com/articles/1204036

contrast_acts$p_adjust <- stats::p.adjust(contrast_acts$p_value, method = 'BH')

contrast_acts |>
  filter(p_adjust < 0.001) |>
  ggplot(aes(x = fct_reorder(factor(source), score), y = score, fill = score)) +
  geom_col() +
  theme_minimal()


# ok i think i should loop through drugs and make a heatmap. werk

contrasts <- list(
  c('drug', 'ctrl', 'dmso'),
  c('drug', '3013', 'dmso'),
  c('drug', '2728', 'dmso'),
  c('drug', '939', 'dmso'),
  c('drug', '985', 'dmso'),
  c('drug', 'sp', 'dmso')
)

acts <- lapply(contrasts, get_tf_activities, data = data7, net = net, subtitle = '7 days')

tf_activities <- purrr::reduce(acts, bind_cols) |>
  dplyr::select(starts_with("source_ctrl") | starts_with("score") | starts_with("p_adjust")) |>
  rename('tf' = 'source_ctrl_vs_dmso') |>
  rename_with( ~ str_replace(.x, "p_adjust", "padj"), cols = starts_with('p_adjust')) |>
  pivot_longer(cols = -tf, names_to = c("metric", 'first', 'vs', 'second'), names_sep = '_') |>
  mutate(cond = str_c(first, vs, second, sep = "_")) |>
  dplyr::select(!c(first, vs, second)) |>
  pivot_wider(names_from = 'metric', values_from = 'value')
  
tfs_i_like <- filter(contrast_acts, p_adjust < 0.001)$source
tfs_i_like <- unique(filter(tf_activities, padj < 0.001)$tf)

tf_activities |>
  filter(tf %in% tfs_i_like) |>
  ggplot(aes(x = score, size = -log10(padj), y = tf, color = cond)) +
  geom_point() +
  theme_minimal()



# ok that plot is maybe kinda dumb. shouldve been making heatmaps lol

tfa_forheatmap <- tf_activities |>
  select(tf, )


# LOOK AT GENES: 

tf <- 'JUN'

tf_targets <- net %>%
      dplyr::filter(source == tf) %>%
      dplyr::arrange(target) %>%
      dplyr::mutate(ID = target, color = "3") %>%
      tibble::column_to_rownames('target')


res |> as.data.frame() |>
  rownames_to_column(var = 'gene') |>
  filter(gene %in% tf_targets$ID) |>
  mutate(`-log10(p_adj)` = -log10(padj)) |>
  ggplot(aes(x = log2FoldChange, y = `-log10(p_adj)`)) +
  geom_point(alpha = 0.1) +
  labs(title = str_c("Genes regulated by", tf, sep = " ")) +
  ggrepel::geom_label_repel(aes(label = gene))
