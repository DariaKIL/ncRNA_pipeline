# src/data_processing.R
# Data loading and preprocessing functions

#' Load and process all data
#'
#' This is the MAIN function - it handles everything:
#' - Loading counts and phenotype
#' - Cleaning sample names
#' - Matching samples
#' - Filtering low counts
#' - Validating consistency
#'
#' @return List with:
#'   - counts: matrix (features x samples)
#'   - phenotype: data.frame with sample and condition columns
#'   - metadata: list with processing stats
load_and_process_data <- function() {
  
  cat("\n" , rep("=", 60), "\n", sep = "")
  cat("LOADING AND PROCESSING DATA\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  # ============================================================================
  # 1. LOAD PHENOTYPE
  # ============================================================================
  
  cat("📋 Loading phenotype data...\n")
  cat("   File:", PHENOTYPE_FILE, "\n")
  
  if (!file.exists(PHENOTYPE_FILE)) {
    stop("❌ Phenotype file not found: ", PHENOTYPE_FILE)
  }
  
  phenotype <- read.delim(
    PHENOTYPE_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Validate required columns
  if (!"sample" %in% colnames(phenotype)) {
    stop("❌ Phenotype file must have 'sample' column")
  }
  if (!"condition" %in% colnames(phenotype)) {
    stop("❌ Phenotype file must have 'condition' column")
  }
  
  cat(sprintf("   ✓ Loaded %d samples\n", nrow(phenotype)))
  cat(sprintf("   ✓ Groups: %s\n", 
              paste(names(table(phenotype$condition)), collapse = ", ")))
  
  # ============================================================================
  # 2. LOAD COUNTS
  # ============================================================================
  
  cat("\n📊 Loading counts data...\n")
  cat("   File:", COUNTS_FILE, "\n")
  
  if (!file.exists(COUNTS_FILE)) {
    stop("❌ Counts file not found: ", COUNTS_FILE)
  }
  
  # Load counts
  counts <- read.table(
    COUNTS_FILE,
    header = TRUE,
    sep = ",",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
  
  cat(sprintf("   ✓ Loaded %d features × %d samples\n", 
              nrow(counts), ncol(counts)))
  
  # ============================================================================
  # 3. CLEAN SAMPLE NAMES
  # ============================================================================
  
  cat("\n🧹 Cleaning sample names...\n")
  
  original_names <- colnames(counts)
  
  # Remove leading X (R adds this to numeric column names)
  colnames(counts) <- gsub("^X", "", colnames(counts))
  
  # Remove any other common prefixes/suffixes if needed
  # colnames(counts) <- gsub("\\..*$", "", colnames(counts))  # Remove after dot
  
  n_changed <- sum(original_names != colnames(counts))
  if (n_changed > 0) {
    cat(sprintf("   ✓ Cleaned %d sample names\n", n_changed))
    cat("   Example: ", original_names[1], " → ", colnames(counts)[1], "\n")
  } else {
    cat("   ✓ No cleaning needed\n")
  }
  
  # ============================================================================
  # 4. MATCH SAMPLES
  # ============================================================================
  
  cat("\n🔗 Matching samples between counts and phenotype...\n")
  
  counts_samples <- colnames(counts)
  phenotype_samples <- phenotype$sample
  
  common_samples <- intersect(counts_samples, phenotype_samples)
  
  cat(sprintf("   Counts:    %d samples\n", length(counts_samples)))
  cat(sprintf("   Phenotype: %d samples\n", length(phenotype_samples)))
  cat(sprintf("   Common:    %d samples\n", length(common_samples)))
  
  if (length(common_samples) == 0) {
    cat("\n❌ ERROR: No matching samples found!\n\n")
    cat("Counts samples (first 10):\n")
    print(head(counts_samples, 10))
    cat("\nPhenotype samples (first 10):\n")
    print(head(phenotype_samples, 10))
    stop("Cannot proceed without matching samples!")
  }
  
  if (length(common_samples) < length(counts_samples)) {
    cat(sprintf("   ⚠️  Warning: %d counts samples not in phenotype\n",
                length(counts_samples) - length(common_samples)))
  }
  
  if (length(common_samples) < length(phenotype_samples)) {
    cat(sprintf("   ⚠️  Warning: %d phenotype samples not in counts\n",
                length(phenotype_samples) - length(common_samples)))
  }
  
  # Filter to common samples
  counts <- counts[, common_samples, drop = FALSE]
  phenotype <- phenotype[phenotype$sample %in% common_samples, ]
  
  # Reorder phenotype to match counts
  phenotype <- phenotype[match(colnames(counts), phenotype$sample), ]
  
  # Validate order
  if (!all(colnames(counts) == phenotype$sample)) {
    stop("❌ ERROR: Sample order mismatch after matching!")
  }
  
  cat("   ✓ Sample order verified\n")
  
  # ============================================================================
  # 5. FILTER LOW COUNTS
  # ============================================================================
  
  cat("\n🔬 Filtering low-expressed features...\n")
  cat(sprintf("   Threshold: >= %d counts in >= %d samples\n",
              MIN_COUNTS, MIN_SAMPLES))
  
  n_before <- nrow(counts)
  
  keep <- rowSums(counts >= MIN_COUNTS) >= MIN_SAMPLES
  counts <- counts[keep, , drop = FALSE]
  
  n_after <- nrow(counts)
  n_removed <- n_before - n_after
  pct_kept <- 100 * n_after / n_before
  
  cat(sprintf("   ✓ Kept %d / %d features (%.1f%%)\n",
              n_after, n_before, pct_kept))
  cat(sprintf("   ✓ Removed %d low-expressed features\n", n_removed))
  
  # ============================================================================
  # 6. FINAL SUMMARY
  # ============================================================================
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("✅ DATA PROCESSING COMPLETE\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  cat(sprintf("Final dataset: %d features × %d samples\n",
              nrow(counts), ncol(counts)))
  cat(sprintf("Groups: %s\n",
              paste(names(table(phenotype$condition)), collapse = ", ")))
  cat(sprintf("Samples per group:\n"))
  print(table(phenotype$condition))
  cat("\n")
  
  # ============================================================================
  # 7. RETURN
  # ============================================================================
  
  metadata <- list(
    n_features_original = n_before,
    n_features_filtered = n_after,
    n_samples = ncol(counts),
    groups = unique(phenotype$condition),
    samples_per_group = as.list(table(phenotype$condition)),
    filter_params = list(
      min_counts = MIN_COUNTS,
      min_samples = MIN_SAMPLES
    )
  )
  
  return(list(
    counts = counts,
    phenotype = phenotype,
    metadata = metadata
  ))
}


#' Annotate RNA types based on feature names
#'
#' @param counts Matrix or data.frame with features as rows
#' @return Vector of RNA types
annotate_rna_types <- function(counts) {
  
  cat("🏷️  Annotating RNA types...\n")
  
  types <- case_when(
    grepl("^hsa-miR", rownames(counts)) ~ "miRNA",
    grepl("^hsa-let", rownames(counts)) ~ "let-7",
    grepl("piR", rownames(counts), ignore.case = TRUE) ~ "piRNA",
    grepl("snoRNA", rownames(counts), ignore.case = TRUE) ~ "snoRNA",
    grepl("tRNA", rownames(counts), ignore.case = TRUE) ~ "tRNA",
    grepl("rRNA", rownames(counts), ignore.case = TRUE) ~ "rRNA",
    grepl("snRNA", rownames(counts), ignore.case = TRUE) ~ "snRNA",
    grepl("misc_RNA", rownames(counts), ignore.case = TRUE) ~ "misc_RNA",
    grepl("Y_RNA", rownames(counts), ignore.case = TRUE) ~ "Y_RNA",
    TRUE ~ "other"
  )
  
  cat("   ✓ Annotation complete:\n")
  type_counts <- table(types)
  for (t in names(type_counts)) {
    cat(sprintf("      %s: %d\n", t, type_counts[t]))
  }
  
  return(types)
}


#' Get highly expressed features per group
#'
#' @param counts Matrix of counts
#' @param phenotype Data frame with sample and condition columns
#' @param threshold Expression threshold (default from config)
#' @return Named list of feature vectors per group
get_highly_expressed_by_group <- function(counts, phenotype,
                                          threshold = HIGH_EXPRESSION_THRESHOLD) {
  
  cat(sprintf("🔝 Finding highly expressed features (threshold: %d)...\n", 
              threshold))
  
  groups <- unique(phenotype$condition)
  result <- list()
  
  for (group in groups) {
    group_samples <- phenotype$sample[phenotype$condition == group]
    group_counts <- counts[, group_samples, drop = FALSE]
    
    # Features with expression > threshold in at least one sample
    highly_expressed <- rownames(group_counts)[
      apply(group_counts, 1, function(x) any(x > threshold))
    ]
    
    result[[group]] <- highly_expressed
    
    cat(sprintf("   %s: %d features\n", group, length(highly_expressed)))
  }
  
  return(result)
}


# ============================================================================
# LEGACY FUNCTIONS (deprecated - use load_and_process_data instead)
# ============================================================================

#' @deprecated Use load_and_process_data() instead
load_phenotypes <- function(file_path) {
  .Deprecated("load_and_process_data")
  message("⚠️  This function is deprecated. Use load_and_process_data() instead.")
}

#' @deprecated Use load_and_process_data() instead
load_counts <- function(file_path) {
  .Deprecated("load_and_process_data")
  message("⚠️  This function is deprecated. Use load_and_process_data() instead.")
}

#' @deprecated Use load_and_process_data() instead
match_samples <- function(counts, phenotable) {
  .Deprecated("load_and_process_data")
  message("⚠️  This function is deprecated. Use load_and_process_data() instead.")
}

#' @deprecated Use load_and_process_data() instead
filter_low_counts <- function(counts, min_counts, min_samples) {
  .Deprecated("load_and_process_data")
  message("⚠️  This function is deprecated. Use load_and_process_data() instead.")
}


cat("✅ Data processing functions loaded\n")

