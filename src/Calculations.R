
#' Create DESeqDataSet object
#'
#' @param counts Count matrix
#' @param coldata Phenotype data
#' @return DESeqDataSet object
create_deseq_dataset <- function(counts, coldata) {
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = "NR")
  
  return(dds)
}

#' Filter dataset
#'
#' @param dds DESeqDataSet object
#' @param min_count Minimum count threshold
#' @param min_samples Minimum samples threshold
#' @return Filtered DESeqDataSet object
filter_dataset <- function(dds, min_count = 10, min_samples = 15) {
  # Filter based on minimum count and samples
  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  dds <- dds[keep,]
  
  return(dds)
}

#' Run differential expression analysis
#'
#' @param dds DESeqDataSet object
#' @return DESeqDataSet object with results
run_differential_expression <- function(dds) {
  # Run DESeq analysis
  dds <- DESeq(dds, fitType = "parametric")
  
  return(dds)
}

#' Get results for specific contrast
#'
#' @param dds DESeqDataSet object
#' @param contrast Contrast specification
#' @return DESeq2 results object
get_contrast_results <- function(dds, contrast) {
  # Get results for specified contrast
  res <- results(dds, contrast = contrast)
  
  return(res)
}

#' Get significant results
#'
#' @param res DESeq2 results object
#' @param alpha Significance threshold
#' @return Data frame with significant results
get_significant_results <- function(res, alpha = 0.05) {
  # Get significant results
  signres <- results(dds, contrast = contrast, alpha = alpha)
  
  return(signres)
}

#' Get sorted results
#'
#' @param res DESeq2 results object
#' @return Data frame with sorted results
get_sorted_results <- function(res) {
  # Sort results by padj
  res_sorted <- as.data.frame(res[order(res$padj), ])
  
  return(res_sorted)
}

#' Get unique miRNAs for each group
#'
#' @param venn_list List of miRNAs for each group
#' @return List of unique miRNAs for each group
get_unique_miRNAs <- function(venn_list) {
  # Get unique miRNAs for each group
  unique_no_c <- setdiff(venn_list$Non, c(venn_list$Сellular, venn_list$AMR, venn_list$CAV))
  unique_cell <- setdiff(venn_list$Сellular, c(venn_list$Non, venn_list$AMR, venn_list$CAV))
  unique_hum <- setdiff(venn_list$AMR, c(venn_list$Non, venn_list$Сellular, venn_list$CAV))
  unique_CAV <- setdiff(venn_list$CAV, c(venn_list$Non, venn_list$Сellular, venn_list$AMR))
  
  return(list(
    NR = unique_no_c,
    ACR = unique_cell,
    AMR = unique_hum,
    CAV = unique_CAV
  ))
}

#' Convert miRNA versions
#'
#' @param miRNA_list List of miRNAs
#' @return Data frame with converted miRNAs
convert_miRNA_versions <- function(miRNA_list) {
  # Convert miRNA versions
  converted <- miRNAVersionConvert(miRNA_list)
  
  return(converted)
}

#' Get multiMiR targets
#'
#' @param mirna List of miRNAs
#' @param org Organism
#' @param table Table type
#' @return Data frame with targets
get_multimir_targets <- function(mirna, org = "hsa", table = "validated") {
  # Get multiMiR targets
  targets <- unique(get_multimir(org = org, mirna = mirna, table = table)@data$target_symbol)
  
  return(targets)
}

#' Perform GO enrichment analysis
#'
#' @param gene_list List of genes
#' @param OrgDb Organism database
#' @param keyType Key type
#' @param ont Ontology type
#' @param pAdjustMethod Method for p-value adjustment
#' @param qvalueCutoff Cutoff for q-value
#' @return GO enrichment result
perform_go_enrichment <- function(gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                               ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05) {
  # Perform GO enrichment
  GO_enrich <- enrichGO(
    gene = gene_list,  
    OrgDb = OrgDb,
    keyType = keyType,
    ont = ont, 
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff
  )
  
  return(GO_enrich)
}

#' Perform KEGG enrichment analysis
#'
#' @param gene_list List of genes
#' @param species Species
#' @return KEGG enrichment result
perform_kegg_enrichment <- function(gene_list, species = "Homo sapiens") {
  # Get KEGG gene sets
  msig_go_bp <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG")
  
  # Perform enrichment
  GO_enrich <- enricher(gene = gene_list, TERM2GENE = msig_go_bp[, c("gs_name", "gene_symbol")])
  
  return(GO_enrich)
}

#' Process mapped data
#'
#' @param mapped_file Path to mapped data file
#' @return Processed data frame
process_mapped_data <- function(mapped_file = "data/mapped.csv") {
  # Load mapped data
  all_df <- read.csv(mapped_file, header = TRUE, sep = ",")
  all_df <- all_df[, -c(1, 2)]
  
  # Convert empty strings to NA
  cols <- c("exact.miRNA", "hairpin.miRNA", "mature.tRNA", "primary.tRNA", "snoRNA", "rRNA", "ncrna.others", "mRNA", "isomiR.miRNA")
  all_df[cols] <- lapply(all_df[cols], function(x) ifelse(x == "", NA, x))
  
  # Create merged column
  all_df$merged_col <- apply(all_df[, c("exact.miRNA", "hairpin.miRNA", "mature.tRNA", "primary.tRNA", "snoRNA", "rRNA", "ncrna.others", "mRNA", "isomiR.miRNA")], 1, function(x) na.omit(x)[1])
  all_df <- all_df[, -c(1:9)]
  
  # Collapse data
  collapsed_df <- all_df %>%
    group_by(merged_col) %>%
    summarise(across(everything(), sum, na.rm = FALSE)) %>%
    as.data.frame()  
  
  # Set row names
  rownames(collapsed_df) <- collapsed_df$merged_col
  collapsed_df$merged_col <- NULL 
  
  # Process column names
  colnames(collapsed_df) <- gsub("^X", "", colnames(collapsed_df))
  rownames(collapsed_df) <- collapsed_df$X
  collapsed_df$X <- NULL
  
  # Filter samples
  common_samples <- intersect(colnames(collapsed_df), coldata$sample)
  collapsed_df <- rownames_to_column(collapsed_df, var = "ncRNA")
  collapsed_df <- collapsed_df[, c("ncRNA", common_samples)] 
  rownames(collapsed_df) <- collapsed_df$ncRNA
  collapsed_df$ncRNA <- NULL
  
  return(collapsed_df)
}