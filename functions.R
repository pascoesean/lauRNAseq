# helper functions
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(msigdbr)

prep_data_ <- function(data, cpm_threshold = 1, num_samples_above = NA){
  # just a placeholder for code i would've used elsewhere
  count_matrix_new <- data |>
    dplyr::select(hgnc_symbol, counts, tdd_id) |>
    pivot_wider(values_from = counts, names_from = tdd_id) |>
    column_to_rownames('hgnc_symbol') 


  # Filter out lowly expressed genes
  # Retain those with greater than `cpm_threshold` counts in all samples
  # all is kinda strict. but so am i.
  cpm_counts=edgeR::cpm(count_matrix_new)
  cpm_more1<-apply(cpm_counts, 1, function(x) sum(x>cpm_threshold))
  idx_keep=which(cpm_more1>=max(cpm_more1))
  if (!is.na(num_samples_above)) {
    idx_keep=which(cpm_more1>=max(num_samples_above))
  }
  count_matrix_new_new <- count_matrix_new[idx_keep,]
  count_matrix_new_new
}

form_dds_obj_ <- function(count_matrix, data, deseq_formula){
    
  coldata <- data |>
    dplyr::select(tdd_id, drug_broad, timepoint, donor, drug) |>
    unique() |>
    mutate(across(everything(), as.factor)) 

    rownames(coldata) <- coldata$tdd_id

    dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = coldata,
                                  design = as.formula(deseq_formula))

  
  dds <- DESeq(dds)
}


get_degs <- function(data, 
                     deseq_formula = '~ drug', 
                     contrast_vector = c('drug', '3013', 'dmso'),
                     plot_volcano_filename = NA, timepointt = NA, print_plot=TRUE){

  # INPUT SHOULD BE THE DATA YOU WANT TO USE FOR DESEQ: IF U WANT A SUBSET OF THE DATA, SUBSET THE DATA.
  # I AM NOT YOUR MAID

  count_matrix_new_new <- prep_data_(data)
  dds <- form_dds_obj_(count_matrix_new_new, data, deseq_formula)

  res <- results(dds, contrast=contrast_vector, independentFiltering = FALSE)

  dmsoDEGs <- results(dds, contrast=contrast_vector) |>
    as.data.frame() |>
    rownames_to_column(var = 'gene') 

  plot <- EnhancedVolcano::EnhancedVolcano(dmsoDEGs,
      lab = dmsoDEGs$gene,
      x = 'log2FoldChange',
      title = str_c(contrast_vector[2], ' vs ', contrast_vector[3]),
      subtitle = str_c('timepoint ', timepointt),
      y = 'padj', pCutoff = 0.05)
    
  if (print_plot){
    print(plot)
  }
  if (!is.na(plot_volcano_filename)) { 
    ggsave(plot_volcano_filename, plot = plot)
  }
    
  # return
  dmsoDEGs
}


get_vst_counts <- function(data, genes_to_select = everything(), 
                           deseq_formula = '~ drug + donor') {

  count_matrix_new_new <- prep_data_(data)
  dds <- form_dds_obj_(count_matrix_new_new, data)
  vsd.fixed <- varianceStabilizingTransformation(dds, blind=FALSE)
  vst_expr <- assay(vsd.fixed)

  vst_data <- vst_expr |>
    as.data.frame() |>
    t() |>
    as.data.frame() |>
    dplyr::select(genes_to_select) |>
    rownames_to_column('tdd_id') |>
    separate_wider_delim(cols = 'tdd_id', delim = '_', names = c('timepoint', 'donor', 'drug')) |>
    mutate(drug_broad = case_when(
           drug == 'dmso' ~ 'dmso',
           drug == 'ctrl' ~ 'ctrl',
           TRUE ~ 'jnki'
         ), .after = drug) 
  
  vst_data
}

run_pca_explorer <- function(data) {
  # this is super cool
  # and is an easy way to show that we don't have drug specific variance in the top 1,000 hvgs
  # maybe try on 7day data only? see what happens
  count_matrix_new_new <- prep_data_(data)
  dds <- form_dds_obj_(count_matrix_new_new, data)
  pcaExplorer::pcaExplorer(dds= dds)
}

get_enrichment <- function(degs, collection = "H", pvalcuttoff = 0.05, plot_filename=NA, contrast_vector=NA, timepointt=NA, print_plot=TRUE) {

  degs <- degs |>
    dplyr::select(gene, log2FoldChange) |>
    arrange(desc(log2FoldChange)) 
  
  m_df <- msigdbr(species = "Homo sapiens", collection = collection)

  gene_list <- degs$log2FoldChange
  names(gene_list) <- degs$gene

  gsea_results <- GSEA(sort(gene_list, decreasing = TRUE),
                     TERM2GENE = m_df[, c("gs_name", "gene_symbol")],
                     TERM2NAME = m_df[, c("gs_name", "gs_name")],
                     pvalueCutoff = pvalcuttoff,
                     minGSSize = 10,
                     maxGSSize = 500)

  # return results
  # only return tibble if its not empty
  if (nrow(gsea_results@result) != 0L) {
  plot <- dotplot(gsea_results, x='NES', 
                    title = str_c(contrast_vector[2], ' vs ', contrast_vector[3], ', timepoint: ', timepointt))
  
  if (print_plot){
    print(plot)
  }

  if (!is.na(plot_filename)) {
    ggsave(plot_filename, plot = plot)
  }
    
  }
  else {
    gsea_results <- NA
  }
  gsea_results
}

add_enrichments <- function(enrichment_results, data, contrast_vector, deseq_formula = '~ drug + donor') {
      degs <- get_degs(data, 
                 deseq_formula = deseq_formula, 
                 contrast_vector = contrast_vector, print_plot=FALSE)

      enrichment <- get_enrichment(degs, pvalcuttoff = 0.5, 
                             contrast_vector=contrast_vector, timepointt=timepointt, print_plot=FALSE)
      
      new_res <- enrichment@result |>
            dplyr::select(ID, NES, p.adjust, core_enrichment) |>
            rename_with(~ str_c(., contrast_vector[2], sep= "_"), .cols = !ID)
      
      enrichment_results |>
            left_join(y = new_res)
}

format_scientific_notation <- function(x) {
  formatC(x, format = "e", digits = 3)
}

run_deg_enrichment <- function(contrast, data, timepointt, deseq_formula = '~donor + drug'){

  print(contrast)
  degs <- get_degs(data, 
                 deseq_formula = deseq_formula, 
                 contrast_vector = contrast,
                 plot_volcano_filename = NA, timepointt = timepointt)
  

  degs_ft <- degs |>
    arrange(padj) |>
    head(n=25) |>
    flextable::flextable() |>
    flextable::set_formatter(pvalue = format_scientific_notation, padj = format_scientific_notation) |>
    flextable::theme_vanilla()

  plot(degs_ft, scaling = 'full')
  # need to add: plotDispEsts(dds)

  enrichment <- get_enrichment(degs, pvalcuttoff = 0.5,
                             plot_filename=NA, 
                             contrast_vector=contrast, timepointt=timepointt)
  
  if (!is.na(enrichment)){
  enrich_print <- enrichment |>
    #arrange(padj) |>
    head(n=25) |>
    dplyr::select(ID, setSize, NES, qvalue) |>
    flextable::flextable() |>
    flextable::set_formatter(qvalue = format_scientific_notation) |>
    flextable::theme_vanilla()

  plot(enrich_print, scaling = 'full')
  }
  else {
    print('no pathways with p < 0.5')
  }

  }


get_tf_activities <- function(contrast_vector = c("drug", "3013", "dmso"), data, net, subtitlee, print_plot = TRUE){


  count_matrix_new_new <- prep_data_(data)
  dds <- form_dds_obj_(count_matrix_new_new, data, deseq_formula)

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

  contrast_acts$p_adjust <- stats::p.adjust(contrast_acts$p_value, method = 'BH')

  plot <- contrast_acts |>
    filter(p_adjust < 0.01) |>
    ggplot(aes(x = fct_reorder(factor(source), score), y = score, fill = score)) +
    geom_col() +
    labs(title = str_c(contrast_vector[2], " vs " , contrast_vector[3]), subtitle = subtitlee) +
    theme_minimal()

  if (print_plot){
    print(plot)
  }

  contrast_acts |>
    rename_with(~ str_c(.x, "_", contrast_vector[2], "_vs_" , contrast_vector[3]))

}
