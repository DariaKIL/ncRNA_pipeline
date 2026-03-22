# config.R


COUNTS_MIR_FILE <- here("data", "raw", "miR.Counts.csv")
COUNTS_TRF_FILE <- here("data", "raw", "tRF.Counts.csv")
ANNOTATION_FILE <- here("data", "raw", "annotation.report.csv")

# FILTERING PARAMETERS
# Minimum expression level
MIN_COUNTS <- 10

# Significance threshold
PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1.0

# Threshold for "high" expression (for Venn diagram)
HIGH_EXPRESSION_THRESHOLD <- 150


# ENRICHMENT ANALYSIS PARAMETERS

# P-value cutoff for enrichment
ENRICH_PVAL_CUTOFF <- 0.05

# Minimum gene set size
MIN_GENESET_SIZE <- 10

# Maximum gene set size
MAX_GENESET_SIZE <- 500

set.seed(42)

# REQUIRED LIBRARIES

required_packages <- c(
  "DESeq2", "ggplot2", "pheatmap", "RColorBrewer", 
  "vsn", "hexbin", "tidyverse", "ggrepel", 
  "EnhancedVolcano", "ggvenn", "clusterProfiler",
  "org.Hs.eg.db", "multiMiR", "miRBaseConverter", "here", 
  "dplyr", "tidyr", "biomaRt", "GenomicRanges", "msigdbr",
  "enrichplot", "Rtsne", "GGally", "rvest", "patchwork"
)

# Function to install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% c("DESeq2", "clusterProfiler", "org.Hs.eg.db", 
                     "EnhancedVolcano", "vsn", "hexbin")) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
}

# Load all libraries
suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
})

cat("✓ Configuration loaded successfully\n")