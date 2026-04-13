###############################################################################
# BLOOD MODEL BUILDING — STEP 1
# QC, Filtering, and Preparation
#
# INPUT:  GSE184050 DESeq2-normalised counts + metadata
# OUTPUT: log_blood_expr_filtered.csv  (genes x samples, top 5000 by variance)
#         sample_qc_metrics.csv
#         plots: library_size, genes_detected, variance_distribution, PCA
#
# KEY DECISIONS vs old pipeline:
#   - NO cell adjustment before WGCNA (it destroys co-expression structure)
#   - Filter to top 5000 variable genes (WGCNA needs this with n=116)
#   - Keep Age + Sex as covariates for model stage (Step 5), not regression here
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(ggplot2)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT    <- "C:/Users/talha/OneDrive/Desktop/gene expression work"
OUT_DIR <- file.path(ROOT, "blood model building")
dir.create(OUT_DIR,                            recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "processed"),    showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "qc"),           showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "plots"),        showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "meta"),         showWarnings=FALSE)

META_FILE <- file.path(ROOT, "gene_expression_phase_3", "meta",
                       "blood_samples_master.csv")
EXPR_FILE <- file.path(ROOT, "gene_expression_phase_3", "raw",
                       "GSE184050_CCHC_T2D_DEseq2_normalized_count.txt.gz")

stopifnot(file.exists(META_FILE))
stopifnot(file.exists(EXPR_FILE))

cat("=======================================================================\n")
cat("  BLOOD MODEL — STEP 1: QC AND FILTERING\n")
cat("=======================================================================\n\n")

###############################################################################
# SECTION 1: LOAD METADATA
###############################################################################

cat("--- Section 1: Loading metadata ---\n")

meta <- fread(META_FILE)
stopifnot("Sample_ID"      %in% names(meta))
stopifnot("Disease_Status" %in% names(meta))
stopifnot(!anyDuplicated(meta$Sample_ID))

meta$Disease_Status <- factor(meta$Disease_Status,
                              levels = c("Control","T2D"))

cat(sprintf("Metadata: %d samples (%d T2D, %d Control)\n",
            nrow(meta),
            sum(meta$Disease_Status=="T2D"),
            sum(meta$Disease_Status=="Control")))
cat(sprintf("Age range: %d-%d, mean=%.1f\n",
            min(meta$Age), max(meta$Age), mean(meta$Age)))
cat(sprintf("Sex: %d female, %d male\n\n",
            sum(meta$Sex=="female"), sum(meta$Sex=="male")))

# Save copy to new folder
fwrite(meta, file.path(OUT_DIR, "meta", "blood_samples_master.csv"))

###############################################################################
# SECTION 2: LOAD EXPRESSION
###############################################################################

cat("--- Section 2: Loading expression matrix ---\n")

expr      <- fread(EXPR_FILE)
gene_ids  <- expr[[1]]
expr_mat  <- as.matrix(expr[, -1])
rownames(expr_mat) <- gene_ids
mode(expr_mat) <- "numeric"

cat(sprintf("Expression matrix: %d genes x %d samples\n\n",
            nrow(expr_mat), ncol(expr_mat)))

###############################################################################
# SECTION 3: ALIGN SAMPLES
###############################################################################

cat("--- Section 3: Aligning samples ---\n")

stopifnot(ncol(expr_mat) == nrow(meta))
colnames(expr_mat) <- meta$Sample_ID
stopifnot(identical(colnames(expr_mat), meta$Sample_ID))

cat(sprintf("Samples aligned: %d\n\n", ncol(expr_mat)))

###############################################################################
# SECTION 4: SAMPLE QC
###############################################################################

cat("--- Section 4: Sample-level QC ---\n")

qc_df <- data.frame(
  Sample_ID      = colnames(expr_mat),
  Disease_Status = meta$Disease_Status,
  LibrarySize    = colSums(expr_mat),
  GenesDetected  = colSums(expr_mat > 0)
)

fwrite(qc_df, file.path(OUT_DIR, "qc", "sample_qc_metrics.csv"))

# Flag outliers: library size > 3 SD from mean
lib_mean <- mean(qc_df$LibrarySize)
lib_sd   <- sd(qc_df$LibrarySize)
outliers <- qc_df$Sample_ID[abs(qc_df$LibrarySize - lib_mean) > 3*lib_sd]

if (length(outliers) > 0) {
  cat(sprintf("WARNING: %d library size outliers detected: %s\n",
              length(outliers), paste(outliers, collapse=", ")))
} else {
  cat("No library size outliers detected (within 3 SD)\n")
}

cat(sprintf("Library size: median=%.0f, range=%.0f-%.0f\n",
            median(qc_df$LibrarySize),
            min(qc_df$LibrarySize),
            max(qc_df$LibrarySize)))
cat(sprintf("Genes detected: median=%.0f\n\n",
            median(qc_df$GenesDetected)))

# Plot 1: Library size
p1 <- ggplot(qc_df, aes(x=LibrarySize, fill=Disease_Status)) +
  geom_histogram(bins=30, alpha=0.7, position="identity") +
  scale_fill_manual(values=c(Control="#2980B9", T2D="#E74C3C")) +
  theme_minimal(base_size=12) +
  labs(title="Library Size Distribution",
       subtitle=sprintf("n=%d samples, median=%.0f",
                        nrow(qc_df), median(qc_df$LibrarySize)),
       x="Library Size (DESeq2-normalised)", y="Count") +
  theme(legend.position="top")
ggsave(file.path(OUT_DIR, "plots", "step1_library_size.png"),
       p1, width=7, height=5, dpi=150)

# Plot 2: Genes detected
p2 <- ggplot(qc_df, aes(x=GenesDetected, fill=Disease_Status)) +
  geom_histogram(bins=30, alpha=0.7, position="identity") +
  scale_fill_manual(values=c(Control="#2980B9", T2D="#E74C3C")) +
  theme_minimal(base_size=12) +
  labs(title="Genes Detected per Sample",
       x="Number of Genes Detected (count > 0)", y="Count") +
  theme(legend.position="top")
ggsave(file.path(OUT_DIR, "plots", "step1_genes_detected.png"),
       p2, width=7, height=5, dpi=150)

###############################################################################
# SECTION 5: GENE FILTERING
###############################################################################

cat("--- Section 5: Gene filtering ---\n")

# Keep genes expressed (>1) in at least 15% of samples
min_samples <- ceiling(0.15 * ncol(expr_mat))
keep_expr   <- rowSums(expr_mat > 1) >= min_samples
cat(sprintf("Genes before filtering: %d\n", nrow(expr_mat)))
cat(sprintf("Genes after expression filter (>1 in ≥15%% samples): %d\n",
            sum(keep_expr)))

expr_filt <- expr_mat[keep_expr, ]

###############################################################################
# SECTION 6: LOG2 TRANSFORMATION
###############################################################################

cat("\n--- Section 6: Log2 transformation ---\n")

log_expr <- log2(expr_filt + 1)
cat("Log2(x+1) transformation applied\n")

###############################################################################
# SECTION 7: SELECT TOP 5000 VARIABLE GENES FOR WGCNA
###############################################################################

cat("\n--- Section 7: Selecting top 5000 variable genes ---\n")

# WHY 5000: WGCNA needs n_samples >> n_genes/module_size
# With 116 samples, >5000 genes leads to network sparsity (as seen before)
# Top 5000 by MAD (more robust than variance for RNA-seq)
gene_mad  <- rowMads(log_expr)
top5000   <- order(gene_mad, decreasing=TRUE)[1:5000]
log_expr_top <- log_expr[top5000, ]

cat(sprintf("Genes after top-5000 MAD filter: %d\n", nrow(log_expr_top)))
cat(sprintf("MAD threshold: %.4f\n", min(gene_mad[top5000])))

# Plot 3: Variance/MAD distribution
mad_df <- data.frame(MAD = gene_mad)
p3 <- ggplot(mad_df, aes(x=MAD)) +
  geom_histogram(bins=50, fill="#85C1E9", color="white") +
  geom_vline(xintercept=min(gene_mad[top5000]),
             color="#E74C3C", linetype="dashed", linewidth=1) +
  annotate("text", x=min(gene_mad[top5000])+0.05,
           y=Inf, vjust=2, hjust=0,
           label=sprintf("Top 5000 cutoff\nMAD=%.3f", min(gene_mad[top5000])),
           color="#E74C3C", size=3.5) +
  theme_minimal(base_size=12) +
  labs(title="Gene Variability Distribution (MAD)",
       subtitle=sprintf("Red line = top 5000 threshold (%d total genes)",
                        nrow(log_expr)),
       x="Median Absolute Deviation (MAD)", y="Number of Genes")
ggsave(file.path(OUT_DIR, "plots", "step1_gene_mad_distribution.png"),
       p3, width=7, height=5, dpi=150)

###############################################################################
# SECTION 8: PCA FOR QC
###############################################################################

cat("\n--- Section 8: PCA quality check ---\n")

pca     <- prcomp(t(log_expr_top), scale.=TRUE)
pca_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(
  PC1            = pca$x[,1],
  PC2            = pca$x[,2],
  PC3            = pca$x[,3],
  Disease_Status = meta$Disease_Status,
  Age            = meta$Age,
  Sex            = meta$Sex
)

# Plot 4a: PCA by disease status
p4a <- ggplot(pca_df, aes(PC1, PC2, color=Disease_Status, shape=Sex)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_manual(values=c(Control="#2980B9", T2D="#E74C3C")) +
  theme_minimal(base_size=12) +
  labs(title="PCA — Blood Expression (Top 5000 genes)",
       subtitle="Coloured by Disease Status, shaped by Sex",
       x=sprintf("PC1 (%.1f%%)", pca_var[1]),
       y=sprintf("PC2 (%.1f%%)", pca_var[2]))
ggsave(file.path(OUT_DIR, "plots", "step1_pca_disease.png"),
       p4a, width=8, height=6, dpi=150)

# Plot 4b: PCA coloured by Age
p4b <- ggplot(pca_df, aes(PC1, PC2, color=Age)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_gradient(low="#AED6F1", high="#1A5276") +
  theme_minimal(base_size=12) +
  labs(title="PCA — Coloured by Age",
       x=sprintf("PC1 (%.1f%%)", pca_var[1]),
       y=sprintf("PC2 (%.1f%%)", pca_var[2]))
ggsave(file.path(OUT_DIR, "plots", "step1_pca_age.png"),
       p4b, width=7, height=5, dpi=150)

cat(sprintf("PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%\n\n",
            pca_var[1], pca_var[2], pca_var[3]))

###############################################################################
# SECTION 9: SAVE OUTPUT
###############################################################################

cat("--- Section 9: Saving filtered expression matrix ---\n")

out_file <- file.path(OUT_DIR, "processed", "log_blood_expr_top5000.csv")
fwrite(data.table(Gene=rownames(log_expr_top), log_expr_top), out_file)

cat(sprintf("Saved: %s\n", out_file))
cat(sprintf("Dimensions: %d genes x %d samples\n\n",
            nrow(log_expr_top), ncol(log_expr_top)))

###############################################################################
# FINAL SUMMARY
###############################################################################

cat("=======================================================================\n")
cat("  STEP 1 SUMMARY\n")
cat("=======================================================================\n")
cat(sprintf("  Input:    %d genes x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("  After expression filter: %d genes\n", sum(keep_expr)))
cat(sprintf("  After top-5000 MAD:      5000 genes\n"))
cat(sprintf("  Output:   log_blood_expr_top5000.csv\n"))
cat(sprintf("  Plots:    step1_library_size.png\n"))
cat(sprintf("            step1_genes_detected.png\n"))
cat(sprintf("            step1_gene_mad_distribution.png\n"))
cat(sprintf("            step1_pca_disease.png\n"))
cat(sprintf("            step1_pca_age.png\n"))
cat("=======================================================================\n")
cat("=== STEP 1 COMPLETED ===\n")