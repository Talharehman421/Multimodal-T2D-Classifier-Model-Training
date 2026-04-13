###############################################################################
# PHASE 1 — STEP 3B
# Module stability analysis via subsampling
###############################################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

allowWGCNAThreads()
options(stringsAsFactors = FALSE)

cat("=== PHASE 1 — STEP 3B: MODULE STABILITY ANALYSIS ===\n")

###############################################################################
# 0) Load WGCNA results
###############################################################################
load("wgcna/step3_wgcna_results.RData")
# Objects: net, moduleColors, MEs, traits

# Load expression matrix
load("processed/step2_outputs_optionC.RData")
datExpr <- t(vst_wgcna)   
softPower <- 8   

gene_names <- colnames(datExpr)

###############################################################################
# 1) Parameters
###############################################################################
n_iter <- 50
subsample_frac <- 0.80
minModuleSize <- 30

set.seed(2025)

dir.create("wgcna/stability", recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2) Reference module assignment
###############################################################################
ref_modules <- moduleColors
names(ref_modules) <- gene_names
ref_levels <- setdiff(unique(ref_modules), "grey")

###############################################################################
# 3) Stability storage
###############################################################################
stability_list <- vector("list", length(ref_levels))
names(stability_list) <- ref_levels

for (m in ref_levels) {
  stability_list[[m]] <- numeric(n_iter)
}

###############################################################################
# 4) Subsampling loop
###############################################################################
cat("Running", n_iter, "subsampling iterations...\n")

for (i in seq_len(n_iter)) {
  
  cat("Iteration", i, "of", n_iter, "\n")
  
  # Subsample samples
  keep_samples <- sample(
    rownames(datExpr),
    size = floor(subsample_frac * nrow(datExpr)),
    replace = FALSE
  )
  
  dat_sub <- datExpr[keep_samples, ]
  
  # Run WGCNA on subset
  net_sub <- blockwiseModules(
    dat_sub,
    power =softPower ,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = minModuleSize,
    deepSplit = 2,
    mergeCutHeight = 0.25,
    numericLabels = FALSE,
    pamRespectsDendro = TRUE,
    saveTOMs = FALSE,
    verbose = 0
  )
  
  sub_modules <- net_sub$colors
  names(sub_modules) <- colnames(dat_sub)
  
  # Compare to reference
  for (m in ref_levels) {
    
    ref_genes <- names(ref_modules)[ref_modules == m]
    sub_genes <- names(sub_modules)[sub_modules == m]
    
    if (length(ref_genes) == 0 || length(sub_genes) == 0) {
      stability_list[[m]][i] <- 0
    } else {
      stability_list[[m]][i] <-
        length(intersect(ref_genes, sub_genes)) /
        length(union(ref_genes, sub_genes))
    }
  }
}

###############################################################################
# 5) Summarize stability
###############################################################################
stability_summary <- tibble(
  Module = ref_levels,
  Mean_Jaccard = sapply(stability_list, mean),
  SD_Jaccard   = sapply(stability_list, sd)
)

write.csv(
  stability_summary,
  "wgcna/stability/module_stability_summary.csv",
  row.names = FALSE
)

###############################################################################
# 6) Identify stable modules
###############################################################################
stable_modules <- stability_summary %>%
  filter(Mean_Jaccard >= 0.70) %>%
  arrange(desc(Mean_Jaccard))

unstable_modules <- stability_summary %>%
  filter(Mean_Jaccard < 0.70)

write.csv(
  stable_modules,
  "wgcna/stability/stable_modules.csv",
  row.names = FALSE
)

write.csv(
  unstable_modules,
  "wgcna/stability/unstable_modules.csv",
  row.names = FALSE
)

###############################################################################
# 7) Save objects
###############################################################################
save(
  stability_list,
  stability_summary,
  stable_modules,
  unstable_modules,
  file = "wgcna/stability/step3b_module_stability.RData"
)

cat("=== STEP 3B COMPLETED SUCCESSFULLY ===\n")
cat("Stable modules:", nrow(stable_modules), "\n")
cat("Unstable modules:", nrow(unstable_modules), "\n")
