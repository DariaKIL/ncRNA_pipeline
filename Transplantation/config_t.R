# Transplantation/config_t.R

PHENOTYPE_FILE <- file.path("Transplantation/data/phenotable.tsv")

OUTPUT_DIR <- "results"
FIGURES_DIR <- "figures"

# Create directories if they do not exist
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)


# Color palette for groups
GROUP_COLORS <- c(
  "NR" = "#4DAF4A",
  "AMR" = "#E41A1C",
  "ACR" = "#377EB8",
  "CAV" = "#FF7F00"
)


CONTRASTS <- list(
  humoral = c("condition", "AMR", "NR"),
  cellular = c("condition", "ACR", "NR"),
  CAV = c("condition", "CAV", "NR")
)

set.seed(42)

cat("✓ Configuration for transplantation loaded successfully\n")