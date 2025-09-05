source('functions.R')
source('wrangling_fromfastq.R')

# RUN WRANGLING THRU LINE 48

#setup variables
contrast_vector <- c('drug', '3013', 'dmso') # control induces x log fold change over dmso
timepointt <- '7d'
deseq_formula <- '~ drug + donor'
filename <- str_c(contrast_vector[2], 'vs', contrast_vector[3], timepointt, sep='_')

################################################
# SANDBOX------
################################################

# JUST PLOT COUNTS!!
# bc this isn't tpm, its really only useful for looking for zeros (and deciding unfiltered deseq is bad)
geneee = "CDC20"
data |>
      filter(symbols == geneee) |>
      ggplot(aes(x = factor(donor), y = tra_counts, color = drug)) +
      geom_point(position=position_jitter(w=0.1,h=0)) + 
      scale_y_log10() +
      labs(title = str_c('raw counts of ', geneee)) +
      facet_wrap(~timepoint) +
      ggpubr::theme_pubr()



count_matrix_new_new <- prep_data_(data, cpm_threshold = 0)
dds <- form_dds_obj_(count_matrix_new_new, data, timepointt, deseq_formula)

#res <- results(dds, contrast=contrast_vector)



mod_mat <- model.matrix(design(dds), colData(dds))

dmso <- colMeans(mod_mat[dds$drug == "dmso", ])
jnki <- colMeans(mod_mat[dds$drug %in% c("3013", "2728", "sp", "939", "985"),])

res_broad <- results(dds, contrast = jnki - dmso, independentFiltering = FALSE)

# CHECK INDEPENDENT FILTERING: doesn't make sense if we have super high cpm thresholds, as we do
as_tibble(metadata(res_broad)$filterNumRej) %>%
  ggplot(aes(x = theta, y = numRej)) +
  geom_point() +
  geom_vline(xintercept = 0.484,
             color = 'red')


dmsoDEGs <- res_broad |>
  as.data.frame() |>
  rownames_to_column(var = 'gene') 

EnhancedVolcano::EnhancedVolcano(dmsoDEGs,
      lab = dmsoDEGs$gene,
      x = 'log2FoldChange',
      title = str_c(contrast_vector[2], ' vs ', contrast_vector[3]),
      subtitle = str_c('timepoint ', timepointt),
      y = 'padj', pCutoff = 0.1, FCcutoff = 0.75)

enrichment <- get_enrichment(degs, pvalcuttoff = 0.5,
                             contrast_vector=contrast_vector, timepointt=timepointt)

genee = "TRAPPC4"
plotCounts(dds, gene=genee, intgroup=c("drug", 'donor'), returnData=TRUE) |>
      ggplot(aes(x = drug, y = count, color = donor)) +
      geom_point(position=position_jitter(w=0.1,h=0)) + 
      scale_y_log10() +
      labs(title = str_c('vst counts of ', genee)) +
      ggpubr::theme_pubr()




ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# SOME OF THE GENE EXPRESSION DIFFERENCES THAT ARE DRIVING VARIATION ARE GENETIC ! ! ! ! ! !! ! ! ! !! 
genez = c("HLA-DRB3","HLA-DQA1", "HLA-DRB5", "HLA-A")
data |>
      filter(hgnc_symbol %in% genez) |>
      ggplot(aes(color = factor(donor), y = counts, x = drug)) +
      geom_point(position=position_jitter(w=0.3,h=0)) + 
      facet_wrap(~hgnc_symbol) +
      scale_y_log10() +
      labs(y = 'raw counts') +
      theme_minimal()

################################################
# END SANDBOX------
################################################

# CODE

data7 <- data |>
      filter(timepoint == timepointt)

degs <- get_degs(data7, 
                 deseq_formula = deseq_formula, 
                 contrast_vector = contrast_vector,
                 plot_volcano_filename = str_c('figures/degs/', filename, '.png'), timepointt = timepointt)


# need to add: plotDispEsts(dds)


enrichment <- get_enrichment(degs, pvalcuttoff = 0.5,
                             plot_filename=str_c('figures/gsea/', filename, '.png'), 
                             contrast_vector=contrast_vector, timepointt=timepointt)


res <- enrichment@result |>
      dplyr::select(ID, NES, p.adjust, core_enrichment) |>
      rename('nes_ctrl' = 'NES', 'padjust_ctrl' = 'p.adjust', 'core_ctrl' = 'core_enrichment')


res <- add_enrichments(res, data7, contrast_vector = c('drug', 'ctrl', 'dmso'))
res <- add_enrichments(res, data7, contrast_vector = c('drug', '2728', 'dmso'))
res <- add_enrichments(res, data7, contrast_vector = c('drug', '939', 'dmso'))
res <- add_enrichments(res, data7, contrast_vector = c('drug', '985', 'dmso'))
res <- add_enrichments(res, data7, contrast_vector = c('drug', 'sp', 'dmso'))
res <- add_enrichments(res, data7, contrast_vector = c('donor', '1519', '1520'))

#heatmap of normalized enrichment scores
for_heatmap <- res |>
      #drop_na() |>
      dplyr::select(ID, starts_with('nes')) |>
      as.data.frame() |>
      column_to_rownames(var = 'ID')

pheatmap::pheatmap(for_heatmap)

interesting_genes <- res |>
      drop_na() |> #play around w including this
      dplyr::select(starts_with('core')) |>
      pivot_longer(everything()) |>
      dplyr::select(value) |>
      separate_wider_delim(value, delim = '/', too_few = 'align_start', names_sep="-") |>
      pivot_longer(everything()) |>
      drop_na() |>
      dplyr::select(value) |>
      unique() 

vst_counts <- get_vst_counts(data, genes_to_select = all_of(interesting_genes$value))

rownames(vst_counts) <- vst_counts |> dplyr::select(donor, drug) |> transmute(dd_id = str_c(donor, drug, sep ="_")) |> pull(dd_id)
counts_4_pheat <- vst_counts |>
      rownames_to_column('ID') |>
      dplyr::select(!c(timepoint, drug, donor, drug_broad)) |>
      as.data.frame() |>
      column_to_rownames('ID') 

pheatmap::pheatmap(counts_4_pheat)


distance_matrix <- dist(counts_4_pheat) |> as.matrix()

####################################
####  ~ LMER INTERLUDE ~ 
pca_obj <- prcomp(counts_4_pheat, scale. = TRUE)

scores <- pca_obj$x |>
  as.data.frame() |>
  rownames_to_column(var = 'tdd_id') |>
  separate_wider_delim(cols = 'tdd_id', delim = '_', names = c('donor', 'drug')) |>
  mutate(donor = factor(donor),
         drug = factor(drug),
         drug_broad = case_when(
          drug == 'dmso' ~ 'dmso',
          drug == 'ctrl' ~ 'ctrl',
          TRUE ~ 'jnki'
         ),
        drug_broad = factor(drug_broad))

ggplot(scores, aes(x = PC2, y = PC4, color = drug_broad, shape = donor)) +
      geom_point() +
      theme_minimal() +
      facet_wrap(~donor)

lmem <- lme4::glmer(drug_broad ~ PC1+PC2+PC3+PC4+ (1|donor), 
                    data = scores,
                    family = binomial,
                    control=lme4::glmerControl(optimizer="bobyqa"))

summary(lmem)


#### ~ END LMER INTERLUDE ~
####################################

