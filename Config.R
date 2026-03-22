# config.R

DATA_DIR <- "data"

PHENOTYPE_FILE <- file.path(DATA_DIR, "phenotable.tsv")
COUNTS_FILE <- file.path(DATA_DIR, "miR.Counts.csv")

OUTPUT_DIR <- "results"
FIGURES_DIR <- "figures"

# Create directories if they do not exist
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)

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

# Color palette for groups
GROUP_COLORS <- c(
  "no_complications" = "#4DAF4A",
  "humoral" = "#E41A1C",
  "cellular" = "#377EB8",
  "TCAD" = "#FF7F00"
)

# ENRICHMENT ANALYSIS PARAMETERS

# P-value cutoff for enrichment
ENRICH_PVAL_CUTOFF <- 0.05

# Minimum gene set size
MIN_GENESET_SIZE <- 10

# Maximum gene set size
MAX_GENESET_SIZE <- 500


CONTRASTS <- list(
  humoral = c("condition", "AMR", "NR"),
  cellular = c("condition", "ACR", "NR"),
  TCAD = c("condition", "CAV", "NR")
)

set.seed(42)

# REQUIRED LIBRARIES

required_packages <- c(
  "DESeq2", "ggplot2", "pheatmap", "RColorBrewer", 
  "vsn", "hexbin", "tidyverse", "ggrepel", 
  "EnhancedVolcano", "VennDiagram", "clusterProfiler",
  "org.Hs.eg.db", "multiMiR", "miRBaseConverter", "ggVennDiagram"
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