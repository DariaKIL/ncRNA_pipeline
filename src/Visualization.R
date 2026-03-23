# Create Venn diagram for miRNA expression

create_venn_diagram <- function(
    counts,
    phenotable,
    condition_col = "condition",
    count_threshold = HIGH_EXPRESSION_THRESHOLD,
    output_file = NULL
) {
  
  # Convert counts to sample-level dataframe
  df <- counts %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
  df$sample <- rownames(df)
  
  # Add condition
  df <- df %>%
    left_join(phenotable[, c("sample", condition_col)], by = "sample")
  
  # Summarise counts per condition
  summed <- df %>%
    group_by(.data[[condition_col]]) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")
  
  # Long format
  long_df <- summed %>%
    pivot_longer(
      cols = -all_of(condition_col),
      names_to = "feature",
      values_to = "counts"
    ) %>%
    filter(counts > count_threshold)
  
  # Build list automatically
  venn_list <- split(long_df$feature, long_df[[condition_col]])
  
  # Plot
  plt <- ggvenn(
    venn_list,
    fill_alpha = 0.5,
    fill_color = GROUP_COLORS,
    stroke_size = 0,
    set_name_size = 5,
    text_size = 4
  )
  
  # Save
  if (!is.null(output_file)) {
    ggsave(output_file, plt, width = 8, height = 6, dpi = 300, bg = "white")
    cat("   ✓ Saved to:", output_file, "\n")
  }
  
  return(list(
    plot = plt,
    venn_list = venn_list
  ))
}


# Create bar plot for all dataset
create_barplot_all_dataset <- function(data, output_file = NULL) {
  # Create long format data
  anno_long <- data %>%
    pivot_longer(cols = -sample, names_to = "RNA_Type", values_to = "Count")
  
  # Create plot
  plt <- ggplot(anno_long, aes(x = sample, y = Count, fill = RNA_Type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Sample", y = "Read Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#Create normalized bar plot for all dataset
create_barplot_all_dataset_normalized <- function(data, output_file = NULL) {
  # Create normalized data
  anno_long <- data %>%
    rowwise() %>%
    mutate(across(-sample, ~ . / sum(c_across(-sample)))) %>% 
    ungroup() %>%
    pivot_longer(cols = -sample, names_to = "RNA_Type", values_to = "Proportion")
  
  # Create plot
  plt <- ggplot(anno_long, aes(x = sample, y = Proportion, fill = RNA_Type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Sample", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

# Create bar plot for dataset by groups

create_grouped_barplots <- function(data, phenotable, condition_col = "condition", output_prefix = NULL) {
  
  # Add condition information from phenotable
  data <- data %>%
    left_join(phenotable[, c("sample", condition_col)], by = "sample")
  
  # Identify columns with numeric count data
  count_cols <- setdiff(colnames(data), c("sample", condition_col))
  
  ## 1. Plot with raw counts
  data_long <- data %>%
    pivot_longer(
      cols = all_of(count_cols),
      names_to = "RNA_Type",
      values_to = "Count"
    )
  
  plt_counts <- ggplot(data_long, aes(x = .data[[condition_col]], y = Count, fill = RNA_Type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Condition", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
  
  # Save raw counts plot if output prefix is provided
  if (!is.null(output_prefix)) {
    ggsave(paste0(output_prefix, "_counts.png"), plt_counts, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  ## 2. Plot with normalized proportions per group
  # First, sum counts within each group
  data_norm_group <- data %>%
    group_by(.data[[condition_col]]) %>%
    summarise(across(all_of(count_cols), sum), .groups = "drop") %>%
    # Normalize counts to proportions within each group
    rowwise() %>%
    mutate(across(all_of(count_cols), ~ . / sum(c_across(all_of(count_cols))))) %>%
    ungroup() %>%
    # Convert to long format for plotting
    pivot_longer(
      cols = all_of(count_cols),
      names_to = "RNA_Type",
      values_to = "Proportion"
    )
  
  plt_prop <- ggplot(data_norm_group, aes(x = .data[[condition_col]], y = Proportion, fill = RNA_Type)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    labs(x = "Condition", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
  
  # Save proportion plot if output prefix is provided
  if (!is.null(output_prefix)) {
    ggsave(paste0(output_prefix, "_proportion.png"), plt_prop, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  # Return both plots in a list
  return(list(counts = plt_counts, proportion = plt_prop))
}



#' Create PCA plot

create_pca_plot <- function(rlt, condition_col = "condition", output_file = NULL) {
  
  # Extract PCA data
  pcaData <- plotPCA(rlt, intgroup = condition_col, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData$sample <- gsub("_.*", "", rownames(pcaData))
  
  # Create plot
  plt <- ggplot(pcaData, aes(PC1, PC2, color = .data[[condition_col]])) +
    geom_point(size = 3) +
    geom_text_repel(
      aes(label = pcaData$sample),
      size = 3,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.4,
      max.overlaps = Inf
    ) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    coord_fixed() +
    theme_bw() +
    scale_color_brewer(palette = "Set2")
  
  # Save if needed
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

# Create heatmap of most variable genes
create_heatmap_most_expressed <- function(rlt, top_n = 50, group_col = "condition", output_file = NULL) {
  
  # Select most variable genes
  vars <- apply(assay(rlt), 1, var)
  select <- order(vars, decreasing = TRUE)[1:top_n]
  
  # Create sample labels: SampleID_Group
  colnames(rlt) <- paste(
    gsub("_.*", "", colnames(rlt)),
    colData(rlt)[[group_col]],
    sep = "_"
  )
  
  # Annotation for heatmap
  df <- as.data.frame(colData(rlt)[group_col])
  
  # Heatmap
  plt <- pheatmap::pheatmap(
    assay(rlt)[select, ],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    annotation_col = df,
    fontsize_row = 6
    #scale = "row"
  )
  
  # Save if needed
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300)
  }
  
  return(plt)
}

#' Create sample distance heatmap

create_sample_distance_heatmap <- function(rlt, group_col = "condition", output_file = NULL) {
  
  # Build sample labels: SampleName_Group
  sample_labels <- paste(gsub("_.*", "", colnames(rlt)), colData(rlt)[[group_col]], sep = "_")
  
  # Calculate Euclidean distances between samples
  sampleDists <- dist(t(assay(rlt)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- sample_labels
  colnames(sampleDistMatrix) <- sample_labels
  
  # Annotation for heatmap
  annotation_df <- data.frame(Group = colData(rlt)[[group_col]])
  rownames(annotation_df) <- sample_labels
  
  # Colors
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  # Plot
  plt <- pheatmap(sampleDistMatrix,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  annotation_col = annotation_df,
                  annotation_row = annotation_df,
                  color = colors)
  
  # Save if requested
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create MA plot
plot_ma_from_res <- function(res_df, contrast_name = NULL, lfc_threshold = 1, alpha = 0.05, output_file = NULL) {
  
  # Set plot title based on contrast name
  plot_title <- ifelse(!is.null(contrast_name), paste0("MA plot: ", contrast_name), "MA plot")
  
  # Filter out genes with baseMean <= 0
  res_df <- res_df[res_df$baseMean > 0, ]
  
  # Determine significance
  res_df$significance <- "not_significant"
  res_df$significance[!is.na(res_df$pvalue) &
                        res_df$pvalue < alpha &
                        abs(res_df$log2FoldChange) > lfc_threshold] <- "pvalue"
  res_df$significance[!is.na(res_df$padj) &
                        res_df$padj < alpha &
                        abs(res_df$log2FoldChange) > lfc_threshold] <- "padj"
  
  # Create plot
  plt <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_x_log10() +
    scale_color_manual(values = c(
      "not_significant" = "gray70",
      "pvalue" = "blue",
      "padj" = "red"
    )) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = c(-lfc_threshold, lfc_threshold), linetype = "dotted") +
    labs(title = plot_title) +   
    theme_bw()
  
  # Save plot if requested
  if(!is.null(output_file)) {
    ggsave(output_file, plt, width = 7, height = 5, dpi = 300)
  }
  
  return(plt)
}

#' Create volcano plot

create_volcano_plot <- function(res_df, res_sign_df, contrast_name = NULL, output_file = NULL) {
  
  # Если контраст указан, используем его как заголовок
  plot_title <- ifelse(!is.null(contrast_name), paste0("Volcano plot: ", contrast_name), NULL)
  
  plt <- EnhancedVolcano(
    res_df,
    lab = rownames(res_df),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    labSize = 3.0,
    boxedLabels = FALSE,
    col = c('black', '#CBD5E8', '#B3E2CD', 'red'),
    colAlpha = 1,
    title = plot_title,
    selectLab = rownames(res_sign_df),
    drawConnectors = TRUE
  )
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create boxplot for specific gene

create_gene_boxplot <- function(dds, gene, output_file = NULL) {
  # Get gene index if gene name provided
  if (is.character(gene)) {
    gene_idx <- which(rownames(dds) == gene)
    gene_name <- gene
  } else {
    gene_idx <- gene
    gene_name <- rownames(dds)[gene_idx]
  }
  
  # Get count data
  count_data <- plotCounts(dds, gene=gene_idx, intgroup="condition", returnData=TRUE)
  count_data$sample <- gsub("_.*", "", rownames(count_data))  
  
  # Create plot
  plt <- ggplot(count_data, aes(x=condition, y=count, fill=condition)) +
    geom_boxplot(alpha=0.5) +
    geom_text(aes(label=sample), size=3, vjust=-1) + 
    geom_jitter(width=0.2, size=2, alpha=0.7) +
    scale_y_log10() +
    ggtitle(paste("Expression of", gene_name)) +
    theme_minimal()
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 16, height = 10, dpi = 300)
  }
  
  return(plt)
}

#' Create heatmap of differentially expressed genes

create_heatmap_diff_expressed <- function(res_sign_df, rlt, coldata, condition_col = "condition",
                                          conditions = NULL, output_file = NULL, scale_rows = TRUE) {
  
  # Check if there are enough significant genes
  if (is.null(res_sign_df) || nrow(res_sign_df) < 2) {
    message("No significant genes to plot for this contrast. Skipping...")
    return(NULL)
  }
  
  # Get expression matrix for significant genes
  de_mat <- assay(rlt)[rownames(res_sign_df), , drop = FALSE]
  
  # Filter samples by selected conditions (or use all if NULL)
  if (!is.null(conditions)) {
    coldata_filtered <- coldata[coldata[[condition_col]] %in% conditions, , drop = FALSE]
  } else {
    coldata_filtered <- coldata
  }
  
  # Keep only samples that exist in the expression matrix
  coldata_filtered <- coldata_filtered[coldata_filtered$sample %in% colnames(de_mat), , drop = FALSE]
  
  # Drop NA in condition column
  coldata_filtered <- coldata_filtered[!is.na(coldata_filtered[[condition_col]]), , drop = FALSE]
  
  # Subset expression matrix
  de_mat_filtered <- de_mat[, coldata_filtered$sample, drop = FALSE]
  
  if (ncol(de_mat_filtered) < 2) {
    message("Not enough samples to plot heatmap. Skipping...")
    return(NULL)
  }
  
  # Scale by row if requested
  datamatrix <- if(scale_rows) t(scale(t(de_mat_filtered))) else de_mat_filtered
  
  # Create annotation dataframe for columns with nice sample names
  annotation_col <- data.frame(condition = coldata_filtered[[condition_col]])
  rownames(annotation_col) <- paste0(gsub("_.*", "", coldata_filtered$sample), "_", coldata_filtered[[condition_col]])
  colnames(datamatrix) <- rownames(annotation_col)
  
  # Only pick colors for levels present in filtered data
  cond_levels <- unique(annotation_col$condition)
  palette_colors <- RColorBrewer::brewer.pal(min(8, length(cond_levels)), "Set2")
  annotation_colors <- list(condition = setNames(palette_colors, cond_levels))
  
  # Create heatmap
  plt <- pheatmap(datamatrix,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  show_rownames = TRUE,
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  fontsize = 12)
  
  # Save plot if needed
  if(!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create GO enrichment dotplot
#'
#' @param go_enrichment GO enrichment result
#' @param title Plot title
#' @param output_file Output file path
#' @return ggplot object
create_go_enrichment_dotplot <- function(go_enrichment, title, output_file = NULL) {
  # Create plot
  p1 <- dotplot(go_enrichment, showCategory = 20) +
    ggtitle(title) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.text.y = element_text(size = 20)
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40))
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p1, width = 16, height = 10, dpi = 300)
  }
  
  return(p1)
}

#' Create GO enrichment emapplot
#'
#' @param go_enrichment GO enrichment result
#' @param title Plot title
#' @param output_file Output file path
#' @return ggplot object
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

#' Create t-SNE plot
#'
#' @param data Expression matrix
#' @param coldata Phenotype data
#' @param output_file Output file path
#' @return ggplot object
create_tsne_plot <- function(data, coldata, output_file = NULL) {
  # Prepare data for t-SNE
  tsne_input <- t(as.matrix(data[, -1]))  # exclude first column
  
  # Run t-SNE
  set.seed(41)  # for reproducibility
  tsne_result <- Rtsne(tsne_input, dims = 2, perplexity = 5, verbose = TRUE, max_iter = 500)
  
  # Create data frame
  tsne_df <- data.frame(
    Sample = rownames(tsne_input),
    tSNE1 = tsne_result$Y[,1],
    tSNE2 = tsne_result$Y[,2]
  )
  
  # Merge with phenotype data
  tsne_df <- merge(tsne_df, coldata, by.x = "Sample", by.y = "sample")
  tsne_df$Sample <- gsub("_.*", "", tsne_df$Sample)
  
  # Create plot
  plt <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = condition)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label=Sample), size=3, vjust=1.5) +
    theme_minimal() +
    labs(title = "t-SNE plot", x = "tSNE 1", y = "tSNE 2")
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}