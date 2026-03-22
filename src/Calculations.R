# Get unique Venn RNA

get_unique_venn_elements <- function(venn_list) {
  unique_list <- lapply(names(venn_list), function(g) {
    others <- setdiff(names(venn_list), g)
    setdiff(venn_list[[g]], unlist(venn_list[others]))
  })
  names(unique_list) <- names(venn_list)
  return(unique_list)
}


# Create DESeqDataSet object

create_dds <- function(counts, coldata, condition_col,
                       ref_level, check_rank = FALSE) {
  
  # Convert column name to formula
  design_formula <- as.formula(paste("~", condition_col))
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = design_formula
  )
  
  # Set reference level
  dds[[condition_col]] <- relevel(dds[[condition_col]], ref = ref_level)
  
  # Optional: check model matrix rank
  if (check_rank) {
    modelMatrix <- model.matrix(design_formula, coldata)
    cat("Rank of model matrix:", qr(modelMatrix)$rank, "\n")
    cat("Number of columns:", ncol(modelMatrix), "\n")
  }
  
  return(dds)
}

# Transformation of dds
auto_transform_dds <- function(dds, blind = FALSE, plot = TRUE) {
  
  n_samples <- ncol(dds)
  n_genes <- nrow(dds)
  
  cat("Samples:", n_samples, "\n")
  cat("Genes:", n_genes, "\n")
  
  # Choose transformation
  if (n_samples > 30) {
    cat("Using variance stabilizing transformation (vst)\n")
    transformed <- varianceStabilizingTransformation(dds, blind = blind)
    method <- "vst"
  } else {
    cat("Using rlog transformation\n")
    transformed <- rlog(dds, blind = blind)
    method <- "rlog"
  }
  
  # Optional mean-SD plot
  if (plot) {
    vsn::meanSdPlot(assay(transformed))
  }
  
  return(list(
    transformed = transformed,
    method = method
  ))
}


# Get significant results

get_significant_results <- function(dds, contrast, output_res = NULL, output_sign_res = NULL, alpha = 0.05, lfc_threshold = 1) {
  
  # Get DESeq2 results
  res <- results(dds, contrast = contrast, alpha = alpha)
  
  # Filter significant genes
  res_sign_df <- as.data.frame(subset(res, padj < alpha & !is.na(padj) & abs(log2FoldChange) > lfc_threshold))
  res_df <- as.data.frame(res)
  
  # Save if output_file specified
  if (!is.null(output_res)) {
    write.csv(res_df, file = output_res, row.names = TRUE)
    cat("Significant results saved to:", output_res, "\n")
  }
  
  if (!is.null(output_sign_res)) {
    write.csv(res_sign_df, file = output_sign_res, row.names = TRUE)
    cat("Significant results saved to:", output_sign_res, "\n")
  }
  
  return(list(
    res_df = res_df,
    res_sign_df = res_sign_df
  ))
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