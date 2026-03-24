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