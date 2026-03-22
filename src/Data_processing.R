# src/data_processing.R
# Data loading and preprocessing functions

load_and_process_data <- function() {
  
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("LOADING AND PROCESSING DATA\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  cat("Loading phenotable data...\n")
  cat("   File:", PHENOTYPE_FILE, "\n")
  
  if (!file.exists(PHENOTYPE_FILE)) {
    stop("❌ Phenotype file not found: ", PHENOTYPE_FILE)
  }
  
  phenotable <- read.delim(
    PHENOTYPE_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Validate required columns
  if (!"sample" %in% colnames(phenotable)) {
    stop("❌ phenotable file must have 'sample' column")
  }
  
  cat(sprintf("   ✓ Loaded %d samples\n", nrow(phenotable)))
  cat(sprintf("   ✓ Groups: %s\n", 
              paste(names(table(phenotable$condition)), collapse = ", ")))
  
  cat("\n Loading counts data...\n")
  cat("   File:", COUNTS_MIR_FILE, "\n")
  
  if (!file.exists(COUNTS_MIR_FILE)) {
    stop("❌ Counts file not found: ", COUNTS_MIR_FILE)
  }
  
  # Load counts
  counts <- read.table(
    COUNTS_MIR_FILE,
    header = TRUE,
    sep = ",",
    row.names = 1,
    check.names = FALSE
  )
  
  cat(sprintf("   ✓ Loaded %d features × %d samples\n", 
              nrow(counts), ncol(counts)))
  
  
  cat("\n Cleaning sample names...\n")
  
  original_names <- colnames(counts)
  
  # Remove leading X (R adds this to numeric column names)
  colnames(counts) <- gsub("^X", "", colnames(counts))
  
  n_changed <- sum(original_names != colnames(counts))
  if (n_changed > 0) {
    cat(sprintf("   ✓ Cleaned %d sample names\n", n_changed))
    cat("   Example: ", original_names[1], " → ", colnames(counts)[1], "\n")
  } else {
    cat("   ✓ No cleaning needed\n")
  }
  
  cat("\n🔗 Matching samples between counts and phenotable...\n")
  
  counts_samples <- colnames(counts)
  phenotable_samples <- phenotable$sample
  
  common_samples <- intersect(counts_samples, phenotable_samples)
  
  cat(sprintf("   Counts:    %d samples\n", length(counts_samples)))
  cat(sprintf("   phenotable: %d samples\n", length(phenotable_samples)))
  cat(sprintf("   Common:    %d samples\n", length(common_samples)))
  
  if (length(common_samples) == 0) {
    cat("\n❌ ERROR: No matching samples found!\n\n")
    cat("Counts samples (first 10):\n")
    print(head(counts_samples, 10))
    cat("\nphenotable samples (first 10):\n")
    print(head(phenotable_samples, 10))
    stop("Cannot proceed without matching samples!")
  }
  
  if (length(common_samples) < length(counts_samples)) {
    cat(sprintf(" %d counts samples not in phenotable\n",
                length(counts_samples) - length(common_samples)))
  }
  
  if (length(common_samples) < length(phenotable_samples)) {
    cat(sprintf(" %d phenotable samples not in counts\n",
                length(phenotable_samples) - length(common_samples)))
  }
  
  # Filter to common samples
  colnames(counts) <- trimws(colnames(counts))
  phenotable$sample <- trimws(phenotable$sample)
  counts <- counts[, common_samples, drop = FALSE]
  phenotable <- phenotable[phenotable$sample %in% common_samples, ]
  
  # Reorder phenotable to match counts
  phenotable <- phenotable[match(colnames(counts), phenotable$sample), ]
  
  # Validate order
  if (!all(colnames(counts) == phenotable$sample)) {
    stop("❌ ERROR: Sample order mismatch after matching!")
  }
  
  cat("   ✓ Sample order verified\n")
  
  cat("\n Filtering low-expressed features...\n")
  
  n_before <- nrow(counts)
  min_group <- min(table(phenotable$condition))
  keep <- rowSums(counts >= MIN_COUNTS) >= min_group
  counts <- counts[keep, , drop = FALSE]
  
  n_after <- nrow(counts)
  n_removed <- n_before - n_after
  pct_kept <- 100 * n_after / n_before
  
  cat(sprintf("   Smallest group size: %d samples\n", min_group))
  cat(sprintf("   Threshold: >= %d counts in >= %d samples\n",
              MIN_COUNTS, min_group))
  cat(sprintf("   ✓ Kept %d / %d features (%.1f%%)\n",
              n_after, n_before, pct_kept))
  cat(sprintf("   ✓ Removed %d low-expressed features\n", n_removed))
  
  # FINAL SUMMARY
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("✅ DATA PROCESSING COMPLETE\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  cat(sprintf("Final dataset: %d features × %d samples\n",
              nrow(counts), ncol(counts)))
  cat(sprintf("Groups: %s\n",
              paste(names(table(phenotable$condition)), collapse = ", ")))
  cat(sprintf("Samples per group:\n"))
  print(table(phenotable$condition))
  cat("\n")
  
  
  metadata <- list(
    n_features_original = n_before,
    n_features_filtered = n_after,
    n_samples = ncol(counts),
    groups = unique(phenotable$condition),
    samples_per_group = as.list(table(phenotable$condition)),
    filter_params = list(
      min_counts = MIN_COUNTS,
      min_samples = min_group
    )
  )
  
  return(list(
    counts = counts,
    phenotable = phenotable,
    metadata = metadata
  ))
}

#' Load and align sample annotation
#'
#' @param annotation_file Path to annotation.report.csv
#' @param phenotable Sample phenotype table
#' @return Annotation table aligned with phenotable
load_sample_annotation <- function(annotation_file, phenotable) {
  
  cat("📋 Loading sample annotation...\n")
  cat("   File:", annotation_file, "\n")
  
  if (!file.exists(annotation_file)) {
    stop("❌ Annotation file not found: ", annotation_file)
  }
  
  anno <- read.csv(annotation_file, header = TRUE)
  
  cat(sprintf("   ✓ Loaded %d samples × %d columns\n", nrow(anno), ncol(anno)))
  
  # Rename sample column for consistency
  colnames(anno)[colnames(anno) == "Sample.name.s."] <- "sample"
  
  # Remove unnecessary columns
  anno <- anno[, -c(2:5, 7, 15)]
  
  cat(sprintf("   ✓ Retained %d columns after filtering\n", ncol(anno)))
  
  cat("\n🔗 Matching annotation with phenotype...\n")
  
  # Clean sample names
  anno$sample <- trimws(anno$sample)
  phenotable$sample <- trimws(phenotable$sample)
  
  common_samples <- intersect(anno$sample, phenotable$sample)
  
  cat(sprintf("   Annotation samples: %d\n", nrow(anno)))
  cat(sprintf("   Phenotype samples:  %d\n", nrow(phenotable)))
  cat(sprintf("   Common samples:     %d\n", length(common_samples)))
  
  if (length(common_samples) == 0) {
    stop("❌ No matching samples between annotation and phenotype")
  }
  
  # Filter annotation
  anno <- anno[anno$sample %in% common_samples, ]
  
  # Reorder to match phenotable
  anno <- anno[match(phenotable$sample, anno$sample), ]
  
  # Verify order
  if (!all(phenotable$sample == anno$sample)) {
    stop("❌ Sample order mismatch after annotation alignment")
  }
  
  cat("   ✓ Annotation successfully aligned\n\n")
  
  return(anno)
}

