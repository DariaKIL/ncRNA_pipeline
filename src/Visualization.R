#' Create Venn diagram for miRNA expression
#'
#' @param data Data frame with miRNA expression data
#' @param condition_col Column name for condition
#' @param count_threshold Minimum count threshold
#' @param output_file Output file path
#' @return ggplot object
create_venn_diagram <- function(data, condition_col = "condition", count_threshold = 150, output_file = NULL) {
  # Process data for Venn diagram
  summed_venn <- data %>%
    group_by(!!sym(condition_col)) %>%
    summarise(across(where(is.numeric), sum)) %>%
    pivot_longer(cols = -!!sym(condition_col), names_to = "miRNA", values_to = "counts") %>%
    filter(counts > count_threshold)
  
  # Extract miRNAs for each condition
  venn_no_c <- summed_venn[summed_venn[[condition_col]] == "NR", ]$miRNA 
  venn_cell <- summed_venn[summed_venn[[condition_col]] == "ACR", ]$miRNA
  venn_hum <- summed_venn[summed_venn[[condition_col]] == "AMR", ]$miRNA 
  venn_CAV <- summed_venn[summed_venn[[condition_col]] == "CAV", ]$miRNA 
  
  # Create Venn list
  venn_list <- list(
    Non = venn_no_c,
    Сellular = venn_cell,
    Humoral = venn_hum,
    CAV = venn_CAV
  )
  
  # Create plot
  plt <- ggvenn(
    venn_list,
    fill_color = c("#4CAF50", "#F44336", "#2196F3", "#FF9800"),
    fill_alpha = 0.5,
    stroke_size = 0,
    set_name_size = 5,
    text_size = 4
  )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create bar plot for all dataset
#'
#' @param data Data frame with annotation data
#' @param output_file Output file path
#' @return ggplot object
create_barplot_all_dataset <- function(data, output_file = NULL) {
  # Create long format data
  anno_long <- data %>%
    pivot_longer(cols = -Sample.name.s., names_to = "RNA_Type", values_to = "Count")
  
  # Create plot
  plt <- ggplot(anno_long, aes(x = Sample.name.s., y = Count, fill = RNA_Type)) +
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

#' Create normalized bar plot for all dataset
#'
#' @param data Data frame with annotation data
#' @param output_file Output file path
#' @return ggplot object
create_barplot_all_dataset_normalized <- function(data, output_file = NULL) {
  # Create normalized data
  anno_long <- data %>%
    rowwise() %>%
    mutate(across(-Sample.name.s., ~ . / sum(c_across(-Sample.name.s.)))) %>% 
    ungroup() %>%
    pivot_longer(cols = -Sample.name.s., names_to = "RNA_Type", values_to = "Proportion")
  
  # Create plot
  plt <- ggplot(anno_long, aes(x = Sample.name.s., y = Proportion, fill = RNA_Type)) +
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

#' Create bar plot for dataset by groups
#'
#' @param data Data frame with summed annotation data
#' @param output_file Output file path
#' @return ggplot object
create_barplot_by_groups <- function(data, output_file = NULL) {
  # Create long format data
  summed_anno_plt <- data %>%
    pivot_longer(cols = -condition, names_to = "RNA_Type", values_to = "Count")
  
  # Create plot
  plt <- ggplot(summed_anno_plt, aes(x = condition, y = Count, fill = RNA_Type)) +
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

#' Create normalization plot
#'
#' @param raw_counts Raw counts matrix
#' @param normalized_counts Normalized counts matrix
#' @param output_file Output file path
#' @return ggplot object
create_normalization_plot <- function(raw_counts, normalized_counts, output_file = NULL) {
  # Create data frame
  df <- data.frame(
    Sample = rep(colnames(raw_counts), 2),
    Counts = c(colSums(raw_counts), colSums(normalized_counts)),
    Type = rep(c("Raw", "Normalized"), each = ncol(raw_counts))
  )
  
  # Create plot
  plt <- ggplot(df, aes(x = Sample, y = Counts, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Counts before and after normalization", x = "Sample", y = "Total Counts") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create PCA plot
#'
#' @param rlt rlog transformed data
#' @param coldata Phenotype data
#' @param output_file Output file path
#' @return ggplot object
create_pca_plot <- function(rlt, coldata, output_file = NULL) {
  # Create PCA data
  pcaData <- plotPCA(rlt, intgroup=c("condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Add sample names
  pcaData$sample <- gsub("_.*", "", coldata$sample)
  
  # Create plot
  plt <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_text(aes(label=sample), size=3, vjust=1.5) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) + 
    coord_fixed() +
    theme_bw() +
    scale_color_brewer(palette = "Set2")
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create PCA plot with multiple components
#'
#' @param rlt rlog transformed data
#' @param coldata Phenotype data
#' @param pcs_count Number of principal components
#' @param output_file Output file path
#' @return ggplot object
create_pca_plot_multiple <- function(rlt, coldata, pcs_count = 4, output_file = NULL) {
  # Create PCA result
  pca_result <- prcomp(t(assay(rlt)), center = TRUE, scale. = TRUE)
  
  # Get explained variance
  explained_variance_ratio <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  
  # Create data frame
  miRNA_pca_df <- as.data.frame(pca_result$x[, 1:pcs_count])
  colnames(miRNA_pca_df) <- paste0("PC", 1:pcs_count, " (", round(explained_variance_ratio[1:pcs_count] * 100, 2), "%)")
  miRNA_pca_df$Group <- paste(coldata$condition)
  
  # Create plot
  plt <- ggpairs(
    miRNA_pca_df, aes(color = Group, alpha = 0.7),
    upper = list(continuous = wrap("cor", size = 3))
  )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create heatmap of most expressed genes
#'
#' @param rlt rlog transformed data
#' @param dds DESeqDataSet object
#' @param output_file Output file path
#' @return pheatmap object
create_heatmap_most_expressed <- function(rlt, dds, output_file = NULL) {
  # Select top 50 genes
  select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:50]
  
  # Create annotation data
  df <- as.data.frame(colData(dds)$condition)
  colnames(df) <- "condition"
  rownames(df) <- colnames(counts(dds))
  
  # Create heatmap
  plt <- pheatmap(assay(rlt)[select,], 
           cluster_rows = TRUE, 
           show_rownames = TRUE, 
           cluster_cols = TRUE, 
           annotation_col = df,
           fontsize_row = 6)
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create sample distance heatmap
#'
#' @param rlt rlog transformed data
#' @param output_file Output file path
#' @return pheatmap object
create_sample_distance_heatmap <- function(rlt, output_file = NULL) {
  # Calculate distances
  sampleDists <- dist(t(assay(rlt)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rlt$condition)
  colnames(sampleDistMatrix) <- paste(rlt$condition)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
  
  # Create heatmap
  plt <- pheatmap(sampleDistMatrix,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           color = colors)
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create MA plot
#'
#' @param res DESeq2 results object
#' @param alpha Significance threshold
#' @param ylim Y-axis limits
#' @param output_file Output file path
#' @return ggplot object
create_ma_plot <- function(res, alpha = 0.05, ylim = c(-8, 8), output_file = NULL) {
  # Create plot
  plt <- plotMA(res, alpha = alpha, ylim = ylim)
  
  # Save if output file specified
  if (!is.null(output_file)) {
    tiff(output_file, width = 8, height = 6, units = "in", res = 300, bg = "white")
    print(plt)
    dev.off()
  }
  
  return(plt)
}

#' Create custom MA plot
#'
#' @param res DESeq2 results object
#' @param output_file Output file path
#' @return ggplot object
create_custom_ma_plot <- function(res, output_file = NULL) {
  # Convert to data frame
  res_df <- res %>%
    as.data.frame() %>%
    mutate(color = case_when( 
      padj < 0.05  ~ "padj < 0.05",   
      pvalue < 0.05  ~ "pvalue < 0.05", 
      TRUE ~ "All"
    ))
  
  # Create plot
  plt <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = color)) +
    geom_point(alpha = 0.7, size = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray40", size = 1.5) +
    scale_color_manual(values = c("All" = "gray70", 
                                  "pvalue < 0.05" = "blue", 
                                  "padj < 0.05" = "red")) +
    scale_x_log10(labels = scales::scientific) + 
    theme_minimal() +
    labs(x = "mean of normalized counts", 
         y = "log fold change", 
         color = NULL)
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create volcano plot
#'
#' @param res DESeq2 results object
#' @param output_file Output file path
#' @return EnhancedVolcano object
create_volcano_plot <- function(res, output_file = NULL) {
  # Create plot
  plt <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    labSize = 3.0,
    boxedLabels = FALSE,
    col = c('black', '#CBD5E8', '#B3E2CD', 'red'),
    colAlpha = 1,
    title = NULL,
    selectLab = rownames(subset(res, padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 1.0)),  # all labels will be displayed
    drawConnectors = TRUE               # arrows if labels overlap
  )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, plot = plt, width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  return(plt)
}

#' Create boxplot for specific gene
#'
#' @param dds DESeqDataSet object
#' @param gene Gene name or index
#' @param output_file Output file path
#' @return ggplot object
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
#'
#' @param res DESeq2 results object
#' @param rlt rlog transformed data
#' @param coldata Phenotype data
#' @param output_file Output file path
#' @return pheatmap object
create_heatmap_diff_expressed <- function(res, rlt, coldata, output_file = NULL) {
  # Filter significant genes
  res_sign <- subset(res, padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 1.0)
  res_sign <- res_sign[order(res_sign$log2FoldChange, decreasing = TRUE), ]
  sig_genes <- rownames(res_sign)  
  
  # Get expression matrix for significant genes
  de_mat <- assay(rlt)[sig_genes, ] 
  
  # Filter samples for specific conditions
  coldata_filtered <- coldata[coldata$condition %in% c("humoral", "no_complications"), ]
  de_mat_filtered <- de_mat[, coldata_filtered$sample]
  datamatrix <- t(scale(t(de_mat_filtered)))
  
  # Create annotation data
  annotation_col <- data.frame(condition = coldata_filtered$condition)
  rownames(annotation_col) <- colnames(datamatrix)
  
  # Create annotation colors
  annotation_colors <- list(
    condition = c("no_complications" = "#FFCC00", "humoral" = "#3399FF"))
  
  # Create heatmap
  plt <- pheatmap(datamatrix,
           cluster_rows = TRUE, 
           show_rownames = TRUE, 
           cluster_cols = TRUE, 
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           display_numbers = FALSE,
           legend = TRUE,
           fontsize = 15)  
  
  # Save if output file specified
  if (!is.null(output_file)) {
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