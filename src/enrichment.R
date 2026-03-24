#' Convert miRNA versions

convert_miRNA_versions <- function(res_sign_df,
                                   contrast_name = NULL,
                                   org = "hsa",
                                   table = "validated",
                                   min_miRNA_per_target = 2) {
  
  # Stop early if no significant miRNAs
  if (is.null(res_sign_df) || nrow(res_sign_df) == 0) {
    message(sprintf("No significant miRNAs found in %s. Skipping.", contrast_name))
    return(NULL)
  }
  
  # Split up/down regulated miRNAs
  up <- res_sign_df %>%
    as.data.frame() %>%
    dplyr::filter(log2FoldChange > 0)
  
  down <- res_sign_df %>%
    as.data.frame() %>%
    dplyr::filter(log2FoldChange < 0)
  
  if (nrow(up) == 0 && nrow(down) == 0) {
    message(sprintf("No up- or down-regulated miRNAs found in %s.", contrast_name))
    return(NULL)
  }
  
  targets_up <- character(0)
  targets_down <- character(0)
  
  if (nrow(up) > 0) {
    mir_up <- rownames(up)
    converted_up <- miRNAVersionConvert(mir_up)
    targets_up <- get_multimir(org = org, mirna = converted_up, table = table)@data %>%
      dplyr::distinct(target_symbol)
  }
  
  if (nrow(down) > 0) {
    mir_down <- rownames(down)
    converted_down <- miRNAVersionConvert(mir_down)
    targets_down <- get_multimir(org = org, mirna = converted_down, table = table)@data%>%
        dplyr::distinct(target_symbol)
  }
  
  return(list(
    down = targets_up,
    up = targets_down
  ))
}


#' Perform GO enrichment analysis

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

perform_kegg_enrichment <- function(gene_list, contrast_name = "", direction = "", species = "Homo sapiens") {
  if(length(gene_list) == 0) {
    message("Empty gene list, skipping KEGG enrichment.")
    return(NULL)
  }
  gene_list <- gene_list$target_symbol
  valid_genes <- gene_list[gene_list %in% keys(org.Hs.eg.db, keytype = "SYMBOL")]
  entrez_ids <- bitr(valid_genes,  fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Hs.eg.db
  )
  
  if(nrow(entrez_ids) == 0) {
    message("No valid SYMBOLs for KEGG conversion. Skipping...")
    return(NULL)
  }
  
  kegg_enrich <- enrichKEGG(
    gene         = entrez_ids$ENTREZID,
    organism     = "hsa",
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 3
  )
  
  if(nrow(kegg_enrich@result) == 0){
    message("No enriched KEGG pathways found.")
    return(NULL)
  }
  
  plt <- dotplot(kegg_enrich, showCategory = 20) +
    ggtitle(sprintf("KEGG enrichment: %s (%s)", contrast_name, direction)) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 11))
  
  return(plt)
}

#' Process mapped data

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

#' Create GO enrichment emapplot

create_go_enrichment_emapplot <- function(go_enrichment, title, output_file = NULL) {
  # Create pairwise termsim
  GO_enrich_BP <- enrichplot::pairwise_termsim(go_enrichment, method = "JC")
  
  # Create plot
  plt <- emapplot(GO_enrich_BP, 
                  repel = TRUE,
                  showCategory = 20) +
    ggtitle(title) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 3)
    )    
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 16, height = 10, dpi = 300)
  }
  
  return(plt)
}