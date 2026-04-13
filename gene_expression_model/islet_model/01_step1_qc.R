###############################################################################
# STEP 1
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

cat("=== PHASE 1 — STEP 1 STARTED ===\n")

###############################################################################
# 1) Set project root
###############################################################################
project_root <- "C:/Users/talha/OneDrive/Desktop/gene expression work/gene Expression Phse 1 "
if (!dir.exists(project_root)) project_root <- getwd()
setwd(project_root)
cat("Project root:", project_root, "\n\n")

###############################################################################
# 2) Ensure directories exist
###############################################################################
dirs <- c("meta", "data", "qc", "logs", "processed")
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 3) Load metadata
###############################################################################
meta_path <- "meta/samples_master.csv"
if (!file.exists(meta_path)) stop("Metadata file missing: meta/samples_master.csv")

cat("Loading metadata...\n")
meta <- fread(meta_path, data.table = FALSE)

required_cols <- c("GSM_ID", "SRR_ID", "Disease_Status")
missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Metadata loaded. Samples:", nrow(meta), "\n\n")

###############################################################################
# 4) Load expression matrix
###############################################################################
expr_path <- "data/GSE50244_Exons_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz"
if (!file.exists(expr_path)) {
  stop("Expression file missing in data/: ", expr_path)
}

cat("Loading expression matrix (fast fread)...\n")
dt <- fread(expr_path, data.table = TRUE)

# First column contains gene IDs
gene_ids <- dt[[1]]
dt[[1]] <- NULL

expr <- as.matrix(dt)
rownames(expr) <- gene_ids

rm(dt); gc()

cat("Expression loaded. Dimensions (genes × samples): ", paste(dim(expr), collapse = " × "), "\n\n")

###############################################################################
# 5) Fix sample names (SRR → GSM) — sequential mapping
###############################################################################

if (ncol(expr) != nrow(meta)) {
  stop("Number of expression columns DOES NOT match metadata rows.\nCheck metadata or expression file.")
}

cat("Assigning SRR_IDs to expression columns based on ordering...\n")
colnames(expr) <- meta$SRR_ID

cat("Mapping SRR → GSM IDs...\n")
colnames(expr) <- meta$GSM_ID[match(colnames(expr), meta$SRR_ID)]

# Reorder metadata to match expression
meta <- meta[match(colnames(expr), meta$GSM_ID), ]

cat("Sample name mapping complete. Matched samples:", ncol(expr), "\n\n")

###############################################################################
#fix NA values
###############################################################################

cat("Converting expression matrix to numeric...\n")
mode(expr) <- "numeric"
expr <- as.matrix(expr)

cat("Checking for NA values...\n")
if (anyNA(expr)) {
  num_na <- sum(is.na(expr))
  cat("Found", num_na, "NA values — replacing with 0.\n")
  expr[is.na(expr)] <- 0
} else {
  cat("No NA values found.\n")
}
cat("NA cleaning complete.\n\n")

###############################################################################
# 7) Remove genes with zero counts across ALL samples
###############################################################################

cat("Identifying genes with zero counts across all samples...\n")
zero_genes <- rowSums(expr == 0) == ncol(expr)
zero_genes[is.na(zero_genes)] <- FALSE  

if (any(zero_genes)) {
  cat("Removing", sum(zero_genes), "genes with zero total counts.\n")
  expr <- expr[!zero_genes, , drop = FALSE]
} else {
  cat("No zero-only genes detected.\n")
}

cat("Remaining genes after zero-filter:", nrow(expr), "\n\n")



save(expr, meta, file = "processed/step1_loaded_expression.RData")
cat("Saved: processed/step1_loaded_expression.RData\n")

cat("=== STEP 1 COMPLETED SUCCESSFULLY ===\n")
