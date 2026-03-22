# src/calculations.R

#' Calculate differential expression statistics
#'
#' @param expr_matrix Expression matrix (genes x samples)
#' @param phenotable Phenotype data with 'condition' column
#' @return data.frame with DE statistics
calculate_differential_expression <- function(expr_matrix, phenotable) {
  
  # Validate inputs
  if (!all(colnames(expr_matrix) == phenotable$sample)) {
    stop("Sample names mismatch between expression matrix and phenotable")
  }
  
  if (!"condition" %in% colnames(phenotable)) {
    stop("phenotable must contain 'condition' column")
  }
  
  # Get condition groups
  conditions <- unique(phenotable$condition)
  if (length(conditions) != 2) {
    stop("Expected exactly 2 conditions, got: ", length(conditions))
  }
  
  group1_samples <- phenotable$sample[phenotable$condition == conditions[1]]
  group2_samples <- phenotable$sample[phenotable$condition == conditions[2]]
  
  # Calculate statistics for each gene
  results <- lapply(1:nrow(expr_matrix), function(i) {
    gene_expr <- expr_matrix[i, ]
    
    group1_values <- gene_expr[group1_samples]
    group2_values <- gene_expr[group2_samples]
    
    # t-test
    test_result <- t.test(group1_values, group2_values)
    
    # Calculate metrics
    data.frame(
      gene = rownames(expr_matrix)[i],
      mean_group1 = mean(group1_values),
      mean_group2 = mean(group2_values),
      log2FC = log2(mean(group2_values) / mean(group1_values)),
      pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results)
  
  # Adjust p-values for multiple testing
  results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
  
  # Add significance flag
  results_df$significant <- results_df$padj < 0.05 & abs(results_df$log2FC) > 1
  
  # Sort by p-value
  results_df <- results_df[order(results_df$pvalue), ]
  rownames(results_df) <- NULL
  
  return(results_df)
}


#' Get top differentially expressed genes
#'
#' @param de_results data.frame from calculate_differential_expression
#' @param n Number of genes to return
#' @param by Sort criterion: "pvalue", "padj", or "log2FC"
#' @return data.frame with top DE genes
get_top_genes <- function(de_results, n = 20, by = "padj") {
  
  if (!by %in% c("pvalue", "padj", "log2FC")) {
    stop("by must be one of: pvalue, padj, log2FC")
  }
  
  if (by == "log2FC") {
    # Sort by absolute log2FC
    de_results <- de_results[order(abs(de_results$log2FC), decreasing = TRUE), ]
  } else {
    de_results <- de_results[order(de_results[[by]]), ]
  }
  
  head(de_results, n)
}


#' Calculate summary statistics
#'
#' @param de_results data.frame from calculate_differential_expression
#' @return list with summary metrics
calculate_summary_stats <- function(de_results) {
  
  list(
    total_genes = nrow(de_results),
    significant_genes = sum(de_results$significant),
    upregulated = sum(de_results$log2FC > 1 & de_results$significant),
    downregulated = sum(de_results$log2FC < -1 & de_results$significant),
    percent_significant = 100 * sum(de_results$significant) / nrow(de_results)
  )
}


#' Filter genes by criteria
#'
#' @param de_results data.frame from calculate_differential_expression
#' @param padj_threshold Adjusted p-value threshold
#' @param log2fc_threshold Log2 fold change threshold
#' @return Filtered data.frame
filter_genes <- function(de_results, 
                         padj_threshold = 0.05, 
                         log2fc_threshold = 1) {
  
  de_results[
    de_results$padj < padj_threshold & 
      abs(de_results$log2FC) > log2fc_threshold,
  ]
}
