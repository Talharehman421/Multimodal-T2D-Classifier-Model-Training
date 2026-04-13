###############################################################################
# PHASE 1 — STEP 4
# Module Eigengene (ME) analysis + robustness assessment
###############################################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

options(stringsAsFactors = FALSE)

cat("=== PHASE 1 — STEP 4: ME ANALYSIS STARTED ===\n")

###############################################################################
# 0) Load required objects
###############################################################################
load("wgcna/step3_wgcna_results.RData")   
load("processed/step2_outputs_optionC.RData")  

###############################################################################
# 1) Prepare trait matrix
###############################################################################
traits <- meta %>%
  transmute(
    Disease = as.numeric(factor(Disease_Status)),
    Age = Age,
    BMI = BMI,
    Sex = ifelse(Sex == "male", 1, 0)
  )

rownames(traits) <- meta$GSM_ID

###############################################################################
# 2) Align MEs and traits
###############################################################################
MEs_ordered <- orderMEs(MEs)

common_samples <- intersect(rownames(MEs_ordered), rownames(traits))
MEs_use <- MEs_ordered[common_samples, ]
traits_use <- traits[common_samples, ]

###############################################################################
# 3) ME–trait correlations 
###############################################################################
ME_trait_cor <- cor(MEs_use, traits_use, use = "p")
ME_trait_p   <- corPvalueStudent(ME_trait_cor, nSamples = nrow(MEs_use))


ME_trait_fdr <- apply(ME_trait_p, 2, p.adjust, method = "fdr")

###############################################################################
# 4) Save ME–trait tables
###############################################################################
dir.create("wgcna/ME_analysis", recursive = TRUE, showWarnings = FALSE)

write.csv(ME_trait_cor,
          "wgcna/ME_analysis/ME_trait_correlations.csv")
write.csv(ME_trait_p,
          "wgcna/ME_analysis/ME_trait_pvalues.csv")
write.csv(ME_trait_fdr,
          "wgcna/ME_analysis/ME_trait_FDR.csv")

###############################################################################
# 5) Robustness: Leave-One-Out (LOO) sensitivity
###############################################################################
cat("Running leave-one-out robustness analysis...\n")

loo_results <- list()

for (i in seq_len(nrow(MEs_use))) {
  
  idx <- setdiff(seq_len(nrow(MEs_use)), i)
  
  cor_loo <- cor(
    MEs_use[idx, ],
    traits_use[idx, ],
    use = "p"
  )
  
  loo_results[[i]] <- cor_loo
}

loo_array <- simplify2array(loo_results)


loo_sd <- apply(loo_array, c(1,2), sd, na.rm = TRUE)

write.csv(
  loo_sd,
  "wgcna/ME_analysis/ME_trait_LOO_SD.csv"
)

###############################################################################
# 6) Save objects
###############################################################################
save(
  MEs_use,
  traits_use,
  ME_trait_cor,
  ME_trait_p,
  ME_trait_fdr,
  loo_sd,
  file = "wgcna/ME_analysis/step4_ME_trait_results.RData"
)

cat("=== STEP 4 COMPLETED SUCCESSFULLY ===\n")
