###############################################################################
#  STEP 2
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(edgeR)
  library(ggplot2)
  library(tidyverse)
})

cat("===  STEP 2 STARTED ===\n")

###############################################################################
# 0) LOAD STEP 1 OUTPUT
###############################################################################
load("processed/step1_loaded_expression.RData")


cat("Loaded Step 1 objects.\n")
cat("Matrix:", nrow(expr), "genes ×", ncol(expr), "samples\n\n")

###############################################################################
# 1) LIBRARY SIZE QC 
###############################################################################
cat("Computing library sizes...\n")

library_sizes <- colSums2(expr)

lib_df <- tibble(
  GSM_ID = colnames(expr),
  Library_Size = library_sizes,
  Disease_Status = meta$Disease_Status
)

ggsave(
  filename = "qc/library_sizes.png",
  ggplot(lib_df, aes(x = reorder(GSM_ID, Library_Size), y = Library_Size, fill = Disease_Status)) +
    geom_col() +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    labs(title="Library Sizes", x="Samples", y="Total Counts"),
  width=10, height=6
)

cat("Library size plot saved.\n\n")
###############################################################################
# FINALIZE FILTERED MATRICES 
###############################################################################

cat("Finalizing matrices...\n")

# Assign final filtered matrices
expr_wgcna_final <- expr_wgcna_final2   
expr_ml_final    <- expr_ml_final2      

cat("Final WGCNA dimensions:", nrow(expr_wgcna_final), "x", ncol(expr_wgcna_final), "\n")
cat("Final ML dimensions:", nrow(expr_ml_final), "x", ncol(expr_ml_final), "\n\n")

###############################################################################
# 4) VST / LOG2 TRANSFORMATION 
###############################################################################

cat("Computing VST/log2 matrices...\n")

vst_wgcna <- log2(expr_wgcna_final + 1)
vst_ml    <- log2(expr_ml_final + 1)

dir.create("processed/vst", recursive = TRUE, showWarnings = FALSE)

write.csv(vst_wgcna, "processed/vst/vst_wgcna.csv")
write.csv(vst_ml,     "processed/vst/vst_ml.csv")

cat("Saved VST matrices.\n\n")

###############################################################################
# 5) PCA ON FINAL WGCNA MATRIX
###############################################################################

cat("Running PCA on final WGCNA matrix...\n")


vst_for_pca <- vst_wgcna[rowVars(vst_wgcna) > 0, , drop = FALSE]

pca_res <- prcomp(t(vst_for_pca), scale. = TRUE)
var_exp <- summary(pca_res)$importance[2, 1:2] * 100

pca_df <- tibble(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  GSM_ID = rownames(pca_res$x)
) %>% left_join(meta, by="GSM_ID")

# Save PCA plot
ggsave(
  "qc/pca_disease.png",
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Disease_Status)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = "PCA - Disease Status (Final Filtered Matrix)",
      x = sprintf("PC1 (%.2f%%)", var_exp[1]),
      y = sprintf("PC2 (%.2f%%)", var_exp[2])
    ) +
    scale_color_brewer(palette = "Set1"),
  width = 8, height = 6
)

cat("PCA saved.\n\n")

###############################################################################
# 6) SAVE ALL CORRECT FINAL OUTPUTS
###############################################################################

save(
  expr_wgcna_final,
  expr_ml_final,
  vst_wgcna,
  vst_ml,
  meta,
  pca_res,
  file = "processed/step2_outputs_optionC.RData"
)

cat("=== STEP 2 COMPLETED SUCCESSFULLY ===\n")
cat("Saved final output: processed/step2_outputs_optionC.RData\n")
