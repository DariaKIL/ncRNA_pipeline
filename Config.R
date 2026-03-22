# config.R

DATA_DIR <- "data"

COUNTS_MIR_FILE <- file.path(DATA_DIR, "/raw/miR.Counts.csv")
COUNTS_TRF_FILE <- file.path(DATA_DIR, "/raw/tRF.Counts.csv")
ANNOTATION_FILE <- file.path(DATA_DIR, "/raw/annotation.report.csv")


# FILTERING PARAMETERS
# Minimum expression level
MIN_COUNTS <- 10

# Significance threshold
PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1.0

# Threshold for "high" expression (for Venn diagram)
HIGH_EXPRESSION_THRESHOLD <- 150

# NORMALIZATION PARAMETERS

# Transformation method: "rlog" or "vst"
TRANSFORM_METHOD <- "rlog"

# VISUALIZATION PARAMETERS

# Number of top genes for heatmap
TOP_GENES_HEATMAP <- 50

# Base text size for plots
BASE_TEXT_SIZE <- 12


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
  "org.Hs.eg.db", "multiMiR", "miRBaseConverter"
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