###############################################################################
# STEP 6


suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(themis)
  library(tidyverse)
})

set.seed(2025)


# SECTION 1: LOAD DATA
###############################################################################

cat("--- Section 1: Loading data ---\n")


load("wgcna/ME_analysis/step4_ME_trait_results.RData")

# Load authoritative metadata
meta_full <- read.csv("meta/samples_master.csv",
                      stringsAsFactors = FALSE)

cat("MEs_use dimensions:", nrow(MEs_use), "samples x", ncol(MEs_use), "MEs\n")
cat("ME names:", paste(colnames(MEs_use), collapse = ", "), "\n\n")


###############################################################################

cat("--- Section 2: Filtering to T2D vs Normal (exact match to step 5) ---\n")

rownames(meta_full) <- meta_full$GSM_ID
common_ids    <- intersect(rownames(MEs_use), rownames(meta_full))
meta_aligned  <- meta_full[common_ids, ]
MEs_aligned   <- MEs_use[common_ids, ]


keep_idx      <- meta_aligned$Disease_Status %in% c("T2D", "Normal")
MEs_62        <- as.matrix(MEs_aligned[keep_idx, ])
mode(MEs_62)  <- "numeric"

y_labels      <- meta_aligned$Disease_Status[keep_idx]
y             <- ifelse(y_labels == "T2D", 1L, 0L)

cat("Samples after filtering:\n")
print(table(y, dnn = "y (0=Normal, 1=T2D)"))
cat("\n")

stopifnot(
  sum(y == 1) == 11,   
  sum(y == 0) == 51,   
  nrow(MEs_62) == 62
)
cat("Sample counts verified: 11 T2D, 51 Normal, 62 total\n\n")

###############################################################################
# SECTION 3: DEFINE BEST HYPERPARAMETERS FROM LOO
###############################################################################

cat("--- Section 3: Hyperparameters from LOO ---\n")

BEST_ALPHA  <- 0     
BEST_LAMBDA <- 0.0271  

cat(sprintf("Best alpha:  %g  (Ridge — selected in ALL 62/62 folds)\n", BEST_ALPHA))
cat(sprintf("Best lambda: %.4f  (median across 62 LOO folds)\n\n", BEST_LAMBDA))

###############################################################################
# SECTION 4: APPLY SAME PREPROCESSING 
###############################################################################

cat("--- Section 4: ME filtering on full 62-sample dataset ---\n")

# Step A: Filter MEs by correlation with y (|r| >= 0.10)
cors <- apply(MEs_62, 2, function(col) {
  if (var(col) == 0) return(0)
  cor(col, y, use = "complete.obs")
})

keep_cor <- which(abs(cors) >= 0.10)

if (length(keep_cor) == 0) {
  keep_cor <- order(abs(cors), decreasing = TRUE)[1:2]
}

cat("MEs passing correlation filter (|r| >= 0.10):\n")
cor_df <- data.frame(ME = names(cors[keep_cor]),
                     Correlation = round(cors[keep_cor], 4))
cor_df <- cor_df[order(abs(cor_df$Correlation), decreasing = TRUE), ]
print(cor_df)
cat("\n")

MEs_filtered <- MEs_62[, keep_cor, drop = FALSE]

# Step B: Remove collinear MEs (|r| > 0.85)
if (ncol(MEs_filtered) > 1) {
  cor_mat <- cor(MEs_filtered, use = "complete.obs")
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- 0
  drop_set <- c()
  for (i in seq_len(nrow(cor_mat))) {
    if (i %in% drop_set) next
    high_j <- which(abs(cor_mat[i, ]) > 0.85)
    drop_set <- union(drop_set, high_j)
  }
  keep_col <- setdiff(seq_len(ncol(MEs_filtered)), drop_set)
  if (length(keep_col) == 0) keep_col <- 1
  MEs_filtered <- MEs_filtered[, keep_col, drop = FALSE]
}

cat("MEs after collinearity removal:\n")
cat(paste(colnames(MEs_filtered), collapse = ", "), "\n")
cat("Total MEs in final model:", ncol(MEs_filtered), "\n\n")

###############################################################################
# SECTION 5: SCALE — FIT ON ALL 62 SAMPLES
# These scaling parameters MUST be saved and used for external validation
###############################################################################

cat("--- Section 5: Scaling (parameters saved for external use) ---\n")

# Fit scaler on all 62 training samples
scaling_means <- colMeans(MEs_filtered)
scaling_sds   <- apply(MEs_filtered, 2, sd)
scaling_sds[scaling_sds == 0] <- 1  

MEs_scaled <- scale(MEs_filtered,
                    center = scaling_means,
                    scale  = scaling_sds)

cat("Scaling parameters (means):\n")
print(round(scaling_means, 4))
cat("\nScaling parameters (SDs):\n")
print(round(scaling_sds, 4))
cat("\n")

###############################################################################
# SECTION 6: APPLY SMOTE ON ALL 62 SAMPLES
# Same parameters as step 5: K=5, over_ratio=0.5
###############################################################################

cat("--- Section 6: SMOTE on full 62-sample training set ---\n")

df_smote     <- as.data.frame(MEs_scaled)
df_smote$.outcome <- factor(ifelse(y == 1, "T2D", "Normal"),
                            levels = c("Normal", "T2D"))

set.seed(2025)
df_smoted <- tryCatch(
  themis::smote(df_smote, var = ".outcome", k = 5, over_ratio = 0.5),
  error = function(e) {
    cat("SMOTE failed:", conditionMessage(e), "— using original data\n")
    df_smote
  }
)

y_sm    <- ifelse(df_smoted$.outcome == "T2D", 1L, 0L)
X_sm    <- as.matrix(df_smoted[, colnames(MEs_scaled), drop = FALSE])
mode(X_sm) <- "numeric"

# Class weights
n1  <- sum(y_sm == 1)
n0  <- sum(y_sm == 0)
w_sm <- ifelse(y_sm == 1, n0 / n1, 1.0)

cat(sprintf("After SMOTE: %d T2D, %d Normal (class weight T2D = %.2fx)\n\n",
            n1, n0, n0 / n1))

###############################################################################
# SECTION 7: FIT FROZEN MODEL
###############################################################################

cat("--- Section 7: Fitting frozen model ---\n")
cat(sprintf("  alpha  = %g  (Ridge)\n", BEST_ALPHA))
cat(sprintf("  lambda = %.4f\n\n", BEST_LAMBDA))

frozen_model <- glmnet(
  x           = X_sm,
  y           = y_sm,
  family      = "binomial",
  alpha       = BEST_ALPHA,
  lambda      = BEST_LAMBDA,
  weights     = w_sm,
  standardize = FALSE  
)

# Extract coefficients
coef_mat  <- as.matrix(coef(frozen_model))
coef_df   <- data.frame(
  Feature     = rownames(coef_mat),
  Coefficient = round(coef_mat[, 1], 6),
  stringsAsFactors = FALSE
)
coef_df   <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ]

cat("Model coefficients (sorted by |magnitude|):\n")
print(coef_df)
cat("\n")

###############################################################################
# SECTION 8: VERIFY — IN-SAMPLE PREDICTIONS ON ORIGINAL 62 SAMPLES
###############################################################################

cat("--- Section 8: In-sample verification (optimistic — for check only) ---\n")

preds_insample <- as.numeric(
  predict(frozen_model, MEs_scaled, type = "response", s = BEST_LAMBDA)
)

roc_insample <- pROC::roc(y, preds_insample,
                          levels = c(0, 1), direction = "<", quiet = TRUE)
auc_insample <- as.numeric(pROC::auc(roc_insample))

cat(sprintf("In-sample AUC: %.3f  (OPTIMISTIC — report LOO AUC = 0.836 instead)\n\n",
            auc_insample))

###############################################################################
# SECTION 9: SAVE ALL OUTPUTS
###############################################################################

cat("--- Section 9: Saving outputs ---\n")

dir.create("ml_smote", showWarnings = FALSE)

# 1. Frozen model
saveRDS(frozen_model, "ml_smote/frozen_model.rds")
cat("Saved: ml_smote/frozen_model.rds\n")

# 2. Scaling parameters — CRITICAL for external validation
scaling_params <- list(
  means        = scaling_means,
  sds          = scaling_sds,
  feature_names = colnames(MEs_filtered),
  alpha        = BEST_ALPHA,
  lambda       = BEST_LAMBDA,
  n_train      = 62,
  n_t2d        = sum(y == 1),
  n_normal     = sum(y == 0),
  loo_auc      = 0.836
)
saveRDS(scaling_params, "ml_smote/frozen_scaling_params.rds")
cat("Saved: ml_smote/frozen_scaling_params.rds\n")

# 3. Feature list with coefficients
write.csv(coef_df,
          "ml_smote/frozen_model_features.csv",
          row.names = FALSE)
cat("Saved: ml_smote/frozen_model_features.csv\n")

# 4. Full save with all objects
save(frozen_model, scaling_params, coef_df,
     MEs_filtered, MEs_scaled, y, y_labels,
     BEST_ALPHA, BEST_LAMBDA,
     file = "ml_smote/step5b_frozen_model.RData")
cat("Saved: ml_smote/step5b_frozen_model.RData\n\n")

###############################################################################
# FINAL SUMMARY
###############################################################################

cat("=======================================================================\n")
cat("  FROZEN MODEL SUMMARY\n")
cat("=======================================================================\n")
cat(sprintf("  Training samples:  62 (11 T2D, 51 Normal)\n"))
cat(sprintf("  Features (MEs):    %d (after filtering)\n", ncol(MEs_filtered))  )
cat(sprintf("  Algorithm:         Ridge Regression (alpha = 0)\n"))
cat(sprintf("  Lambda:            %.4f\n", BEST_LAMBDA))
cat(sprintf("  LOO-CV AUC:        0.836  (honest, unbiased estimate)\n"))
cat(sprintf("  In-sample AUC:     %.3f  (optimistic — do not report)\n\n",
            auc_insample))
cat("  Next step: run 07_phase2_external_validation.R\n")
cat("=======================================================================\n")
cat("=== STEP 5b COMPLETED SUCCESSFULLY ===\n")