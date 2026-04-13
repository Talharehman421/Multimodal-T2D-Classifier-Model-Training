###############################################################################
# BLOOD MODEL BUILDING — STEP 4
# Build Model Input Object (Step 1 equivalent for Phase 5)
#
# INPUT:  WGCNA_blood_network.RData  (from Step 2)
#         blood_hub_genes.csv        (from Step 3)
#         blood_samples_master.csv   (metadata)
#
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

ROOT      <- "C:/Users/talha/OneDrive/Desktop/gene expression work"
BM_DIR    <- file.path(ROOT, "blood model building")
WGCNA_DIR <- file.path(BM_DIR, "wgcna")
HUB_DIR   <- file.path(BM_DIR, "hubs")
OBJ_DIR   <- file.path(BM_DIR, "model_input")
dir.create(OBJ_DIR,                        showWarnings=FALSE)
dir.create(file.path(OBJ_DIR, "plots"),    showWarnings=FALSE)

cat("=======================================================================\n")
cat("  BLOOD MODEL — STEP 4: BUILD MODEL INPUT OBJECT\n")
cat("=======================================================================\n\n")

###############################################################################
# SECTION 1: LOAD ALL INPUTS
###############################################################################

cat("--- Section 1: Loading inputs ---\n")

load(file.path(WGCNA_DIR, "WGCNA_blood_network.RData"))
# loads: net, datExpr, softPower, MEs, module_colors

meta <- fread(file.path(BM_DIR, "meta", "blood_samples_master.csv"))
meta$Disease_Status <- factor(meta$Disease_Status,
                              levels=c("Control","T2D"))
meta$Sex_bin <- ifelse(meta$Sex == "female", 0L, 1L)  

hub_genes <- fread(file.path(HUB_DIR, "blood_hub_genes.csv"))

cat(sprintf("Samples: %d\n", nrow(meta)))
cat(sprintf("MEs available: %d\n", ncol(MEs)))
cat(sprintf("Hub genes: %d across %d modules\n\n",
            nrow(hub_genes), length(unique(hub_genes$Module))))

###############################################################################
# SECTION 2: ALIGN METADATA WITH WGCNA SAMPLES
###############################################################################

cat("--- Section 2: Aligning samples ---\n")

wgcna_samples <- rownames(datExpr)
stopifnot(all(wgcna_samples %in% meta$Sample_ID))

meta_aligned <- meta[match(wgcna_samples, meta$Sample_ID)]
stopifnot(identical(meta_aligned$Sample_ID, wgcna_samples))

cat(sprintf("Samples aligned: %d\n\n", nrow(meta_aligned)))

###############################################################################
# SECTION 3: BUILD X MATRIX (MODULE EIGENGENES ONLY, NO GREY)
###############################################################################

cat("--- Section 3: Building feature matrix X ---\n")

# Remove grey module eigengene if present
ME_cols  <- colnames(MEs)
ME_keep  <- ME_cols[!grepl("grey", ME_cols, ignore.case=TRUE)]
X_me     <- as.matrix(MEs[, ME_keep, drop=FALSE])
rownames(X_me) <- rownames(MEs)

cat(sprintf("MEs used (grey excluded): %d\n", ncol(X_me)))
cat(sprintf("ME names: %s\n\n", paste(ME_keep, collapse=", ")))

# Sanity check: variance should not be near-zero
me_vars <- apply(X_me, 2, var)
cat("ME variances:\n")
print(round(me_vars, 4))

if (any(me_vars < 0.001)) {
  cat("\nWARNING: Some ME variances near zero:")
  cat(names(me_vars)[me_vars < 0.001], "\n")
  cat("These may not be informative for modelling.\n\n")
} else {
  cat("\nAll ME variances healthy (>0.001). Good.\n\n")
}

###############################################################################
# SECTION 4: OUTCOME VECTOR
###############################################################################

y <- meta_aligned$Disease_Status
names(y) <- meta_aligned$Sample_ID

cat(sprintf("y: %d Control, %d T2D\n\n",
            sum(y=="Control"), sum(y=="T2D")))

###############################################################################
# SECTION 5: ME-TRAIT CORRELATIONS
###############################################################################

cat("--- Section 5: ME-trait correlations ---\n")

y_int    <- ifelse(y == "T2D", 1L, 0L)
age_vec  <- meta_aligned$Age
sex_vec  <- meta_aligned$Sex_bin

me_cor_disease <- cor(X_me, y_int,   use="complete.obs")
me_cor_age     <- cor(X_me, age_vec, use="complete.obs")
me_cor_sex     <- cor(X_me, sex_vec, use="complete.obs")

trait_cor <- data.table(
  ME             = ME_keep,
  Cor_Disease    = round(as.numeric(me_cor_disease), 4),
  Cor_Age        = round(as.numeric(me_cor_age),     4),
  Cor_Sex        = round(as.numeric(me_cor_sex),     4)
)

cat("ME correlations with traits:\n")
print(trait_cor)

fwrite(trait_cor, file.path(OBJ_DIR, "ME_trait_correlations.csv"))


n_disease_me <- sum(abs(trait_cor$Cor_Disease) >= 0.10)
cat(sprintf("\nMEs with |r|>=0.10 with disease: %d / %d\n\n",
            n_disease_me, nrow(trait_cor)))

###############################################################################
# SECTION 6: PLOTS
###############################################################################

cat("--- Section 6: Generating plots ---\n")

me_long <- as.data.frame(X_me)
me_long$Sample_ID      <- rownames(me_long)
me_long$Disease_Status <- meta_aligned$Disease_Status

me_melt <- reshape(me_long,
                   varying   = ME_keep,
                   v.names   = "Eigengene",
                   timevar   = "ME",
                   times     = ME_keep,
                   direction = "long")

# Plot 1: ME boxplots by disease status
p1 <- ggplot(me_melt, aes(x=Disease_Status, y=Eigengene,
                          fill=Disease_Status)) +
  geom_boxplot(outlier.size=0.8, alpha=0.8) +
  geom_jitter(width=0.15, size=0.5, alpha=0.4) +
  facet_wrap(~ME, scales="free_y") +
  scale_fill_manual(values=c(Control="#2980B9", T2D="#E74C3C")) +
  theme_minimal(base_size=11) +
  theme(legend.position="top",
        axis.text.x=element_text(angle=30, hjust=1)) +
  labs(title="Module Eigengenes by Disease Status",
       x="", y="Eigengene Value")
ggsave(file.path(OBJ_DIR, "plots", "step4_ME_boxplots.png"),
       p1, width=max(10, ncol(X_me)*2), height=6, dpi=150)

# Plot 2: Heatmap of ME-trait correlations
trait_long <- melt(trait_cor, id.vars="ME",
                   variable.name="Trait", value.name="Correlation")

p2 <- ggplot(trait_long, aes(x=Trait, y=ME, fill=Correlation)) +
  geom_tile(color="white") +
  geom_text(aes(label=sprintf("%.2f", Correlation)),
            size=3.5, color="black") +
  scale_fill_gradient2(low="#2980B9", mid="white", high="#E74C3C",
                       midpoint=0, limits=c(-1,1)) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1)) +
  labs(title="ME–Trait Correlation Heatmap",
       subtitle="Disease = T2D (1) vs Control (0)",
       x="Trait", y="Module Eigengene")
ggsave(file.path(OBJ_DIR, "plots", "step4_ME_trait_heatmap.png"),
       p2, width=7, height=max(5, nrow(trait_cor)*0.5+2), dpi=150)

# Plot 3: Age and Sex distribution by disease (covariate check)
cov_df <- data.frame(
  Sample_ID      = meta_aligned$Sample_ID,
  Disease_Status = meta_aligned$Disease_Status,
  Age            = meta_aligned$Age,
  Sex            = meta_aligned$Sex
)

p3a <- ggplot(cov_df, aes(x=Disease_Status, y=Age, fill=Disease_Status)) +
  geom_boxplot(alpha=0.8) +
  geom_jitter(width=0.1, size=1, alpha=0.5) +
  scale_fill_manual(values=c(Control="#2980B9", T2D="#E74C3C")) +
  theme_minimal(base_size=12) +
  labs(title="Age by Disease Status", x="", y="Age (years)")

p3b <- ggplot(cov_df, aes(x=Disease_Status, fill=Sex)) +
  geom_bar(position="fill", alpha=0.85) +
  scale_fill_manual(values=c(female="#F1948A", male="#85C1E9")) +
  theme_minimal(base_size=12) +
  labs(title="Sex by Disease Status", x="", y="Proportion")

library(gridExtra)
p3_combined <- arrangeGrob(p3a, p3b, ncol=2)
ggsave(file.path(OBJ_DIR, "plots", "step4_covariate_check.png"),
       p3_combined, width=10, height=5, dpi=150)

###############################################################################
# SECTION 7: SAVE RDS OBJECT
###############################################################################

cat("--- Section 7: Saving model input object ---\n")

model_obj <- list(
  X             = X_me,                    
  y             = y,                        
  metadata      = meta_aligned,           
  feature_names = ME_keep,                 
  hub_genes     = hub_genes,              
  trait_cors    = trait_cor,               
  wgcna_summary = list(
    softPower       = softPower,
    n_modules       = length(ME_keep),
    minModuleSize   = 20,
    mergeCutHeight  = 0.30
  )
)

saveRDS(model_obj,
        file.path(OBJ_DIR, "blood_model_step1_object.rds"))

cat(sprintf("Saved: blood_model_step1_object.rds\n"))
cat(sprintf("  X: %d samples x %d features\n", nrow(X_me), ncol(X_me)))
cat(sprintf("  y: %d Control, %d T2D\n\n",
            sum(y=="Control"), sum(y=="T2D")))

###############################################################################
# FINAL SUMMARY
###############################################################################

cat("=======================================================================\n")
cat("  STEP 4 SUMMARY\n")
cat("=======================================================================\n")
cat(sprintf("  Samples:          %d\n", nrow(X_me)))
cat(sprintf("  Features (MEs):   %d\n", ncol(X_me)))
cat(sprintf("  T2D:              %d\n", sum(y=="T2D")))
cat(sprintf("  Control:          %d\n", sum(y=="Control")))
cat(sprintf("  Covariates stored: Age, Sex (used in Step 5 model)\n"))
cat(sprintf("  MEs with |r|>0.10 with disease: %d\n", n_disease_me))
cat("  \n  Next step: Run BM_step5_model.R\n")
cat("=======================================================================\n")
cat("=== STEP 4 COMPLETED ===\n")