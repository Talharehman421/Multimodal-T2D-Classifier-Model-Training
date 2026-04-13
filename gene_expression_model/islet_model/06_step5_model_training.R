###############################################################################
# STEP 5 
# Elastic Net with SMOTE + Class Weights + Proper Nested CV
###############################################################################

suppressPackageStartupMessages({
  library(glmnet)   # elastic net
  library(pROC)     # ROC-AUC
  library(PRROC)    # PR-AUC
  library(themis)   # SMOTE
  library(tidyverse)
})

set.seed(2025)

cat("=======================================================================\n")
cat("  PHASE 1 — STEP 5 (REVISED): Elastic Net + SMOTE + Class Weights\n")
cat("=======================================================================\n\n")

###############################################################################
# SECTION 1: LOAD DATA
###############################################################################
cat("--- Section 1: Loading data ---\n")

load("wgcna/ME_analysis/step4_ME_trait_results.RData")

stopifnot(
  exists("MEs_use"),
  exists("traits_use")
)

meta_full <- read.csv("meta/samples_master.csv",
                      stringsAsFactors = FALSE)

cat("Samples in MEs_use:", nrow(MEs_use), "\n")
cat("Columns in MEs_use:", paste(colnames(MEs_use), collapse = ", "), "\n\n")

###############################################################################
# SECTION 2: CORRECT LABEL ENCODING
###############################################################################

cat("--- Section 2: Label encoding ---\n")

# Align metadata to MEs_use row order using GSM_ID
rownames(meta_full) <- meta_full$GSM_ID

common_ids <- intersect(rownames(MEs_use), rownames(meta_full))
meta_aligned <- meta_full[common_ids, ]
MEs_aligned  <- MEs_use[common_ids, ]

cat("Original Disease_Status distribution:\n")
print(table(meta_aligned$Disease_Status, useNA = "always"))
cat("\n")

# Keep ONLY T2D and Normal — exclude IGT and NA
keep_idx <- meta_aligned$Disease_Status %in% c("T2D", "Normal")

cat("Samples kept (T2D + Normal):", sum(keep_idx), "\n")
cat("Samples excluded (IGT):",
    sum(meta_aligned$Disease_Status == "IGT", na.rm = TRUE), "\n")
cat("Samples excluded (NA):",
    sum(is.na(meta_aligned$Disease_Status)), "\n\n")

MEs_binary  <- as.matrix(MEs_aligned[keep_idx, ])
mode(MEs_binary) <- "numeric"

# Explicit encoding: T2D = 1 (positive/case), Normal = 0 (negative/control)
y_labels <- meta_aligned$Disease_Status[keep_idx]
y <- ifelse(y_labels == "T2D", 1L, 0L)

cat("Final label distribution:\n")
print(table(y, dnn = "y (0=Normal, 1=T2D)"))
cat("\n")

n_t2d   <- sum(y == 1)
n_norm  <- sum(y == 0)
n_total <- length(y)

stopifnot(
  n_t2d > 0,
  n_norm > 0,
  !anyNA(y),
  nrow(MEs_binary) == n_total
)

###############################################################################
# SECTION 3: HELPER FUNCTIONS
###############################################################################

cat("--- Section 3: Defining helper functions ---\n\n")


filter_MEs_by_correlation <- function(X_train, y_train,
                                      cor_threshold = 0.10) {
  cors <- apply(X_train, 2, function(col) {
    if (var(col) == 0) return(0)
    cor(col, y_train, use = "complete.obs")
  })
  keep <- abs(cors) >= cor_threshold
  if (sum(keep) == 0) {

    top2 <- order(abs(cors), decreasing = TRUE)[1:min(2, ncol(X_train))]
    keep[top2] <- TRUE
  }
  return(which(keep))
}

remove_collinear_MEs <- function(X_train, cor_cutoff = 0.85) {
  if (ncol(X_train) <= 1) return(1:ncol(X_train))
  cor_mat  <- cor(X_train, use = "complete.obs")
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- 0
  drop_set <- c()
  for (i in seq_len(nrow(cor_mat))) {
    if (i %in% drop_set) next
    high_j <- which(abs(cor_mat[i, ]) > cor_cutoff)
    drop_set <- union(drop_set, high_j)
  }
  keep <- setdiff(seq_len(ncol(X_train)), drop_set)
  if (length(keep) == 0) keep <- 1  # safety: always keep at least 1
  return(keep)
}


apply_smote <- function(X_train, y_train, k = 5, over_ratio = 0.5) {
  # themis::smote() requires a data.frame with the outcome as a factor
  df_train <- as.data.frame(X_train)
  df_train$.outcome <- factor(ifelse(y_train == 1, "T2D", "Normal"),
                              levels = c("Normal", "T2D"))
  
  n_minority <- sum(y_train == 1)
  k_use <- min(k, n_minority - 1)
  if (k_use < 1) {
    cat("    [SMOTE skipped: not enough minority samples]\n")
    df_out <- df_train
  } else {
    df_out <- tryCatch(
      themis::smote(df_train,
                    var      = ".outcome",
                    k        = k_use,
                    over_ratio = over_ratio),
      error = function(e) {
        cat("    [SMOTE error:", conditionMessage(e), "— skipping]\n")
        df_train
      }
    )
  }
  
  y_out <- ifelse(df_out$.outcome == "T2D", 1L, 0L)
  X_out <- as.matrix(df_out[, colnames(X_train), drop = FALSE])
  mode(X_out) <- "numeric"
  
  list(X = X_out, y = y_out)
}

# ------------------------------------------------------------------
# compute_class_weights()
# Returns per-sample weights vector for glmnet: upweights T2D samples
# Weight for T2D = n_normal / n_t2d  (inverse frequency weighting)
# ------------------------------------------------------------------
compute_class_weights <- function(y_vec) {
  n1 <- sum(y_vec == 1)  # T2D
  n0 <- sum(y_vec == 0)  # Normal
  w  <- ifelse(y_vec == 1, n0 / n1, 1.0)
  return(w)
}

# ------------------------------------------------------------------
# stratified_kfold_ids()
# Creates stratified k-fold indices preserving class proportions
# Returns list of length k, each element = test indices for that fold
# ------------------------------------------------------------------
stratified_kfold_ids <- function(y_vec, k = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx0 <- which(y_vec == 0)
  idx1 <- which(y_vec == 1)
  # Shuffle within each class
  idx0 <- sample(idx0)
  idx1 <- sample(idx1)
  # Assign fold IDs
  fold_id <- integer(length(y_vec))
  fold_id[idx0] <- ((seq_along(idx0) - 1) %% k) + 1
  fold_id[idx1] <- ((seq_along(idx1) - 1) %% k) + 1
  # Return test indices per fold
  lapply(1:k, function(f) which(fold_id == f))
}

# ------------------------------------------------------------------
# inner_cv_tune()
# Runs stratified 10-fold CV on training data to find best (alpha, lambda)
# Returns: list(best_alpha, best_lambda)
# ------------------------------------------------------------------
inner_cv_tune <- function(X_tr, y_tr,
                          alpha_grid = c(0, 0.25, 0.5, 0.75, 1.0),
                          k_inner    = 10,
                          smote_k    = 5,
                          smote_ratio = 0.5,
                          seed       = 42) {
  
  inner_folds <- stratified_kfold_ids(y_tr, k = k_inner, seed = seed)
  
  best_auc    <- -Inf
  best_alpha  <- 0.5
  best_lambda <- NULL
  
  for (alpha_val in alpha_grid) {
    
    # Collect predictions across inner folds for this alpha
    # We use the glmnet lambda path — find lambda maximising AUC
    lambda_aucs <- list()
    
    for (fold_i in seq_along(inner_folds)) {
      
      val_idx   <- inner_folds[[fold_i]]
      train_idx <- setdiff(seq_len(nrow(X_tr)), val_idx)
      
      X_inner_tr <- X_tr[train_idx, , drop = FALSE]
      y_inner_tr <- y_tr[train_idx]
      X_inner_val <- X_tr[val_idx, , drop = FALSE]
      y_inner_val  <- y_tr[val_idx]
      
      # Scale: fit on inner train, apply to inner val
      tr_mean <- colMeans(X_inner_tr)
      tr_sd   <- apply(X_inner_tr, 2, sd)
      tr_sd[tr_sd == 0] <- 1  # avoid division by zero
      
      X_inner_tr_sc  <- scale(X_inner_tr,  center = tr_mean, scale = tr_sd)
      X_inner_val_sc <- scale(X_inner_val, center = tr_mean, scale = tr_sd)
      
      # SMOTE on inner training fold only
      smoted <- apply_smote(X_inner_tr_sc, y_inner_tr,
                            k = smote_k, over_ratio = smote_ratio)
      X_sm <- smoted$X
      y_sm <- smoted$y
      
      # Class weights
      w_sm <- compute_class_weights(y_sm)
      
      # Fit elastic net on full lambda path
      # Suppress convergence warnings for speed
      fit <- tryCatch(
        suppressWarnings(
          glmnet(x       = X_sm,
                 y       = y_sm,
                 family  = "binomial",
                 alpha   = alpha_val,
                 weights = w_sm,
                 standardize = FALSE)  # already scaled
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      
      # Predict on inner validation fold
      preds <- predict(fit, X_inner_val_sc, type = "response")
      # preds: matrix rows=val samples, cols=lambda values
      
      # For each lambda, compute AUC (only if both classes present in val)
      if (length(unique(y_inner_val)) < 2) next
      
      lambda_auc_vec <- apply(preds, 2, function(p) {
        tryCatch(
          as.numeric(pROC::auc(
            pROC::roc(y_inner_val, p,
                      levels = c(0, 1), direction = "<",
                      quiet = TRUE)
          )),
          error = function(e) NA_real_
        )
      })
      
      lambda_aucs[[fold_i]] <- data.frame(
        lambda    = fit$lambda,
        auc       = lambda_auc_vec,
        fold      = fold_i,
        alpha     = alpha_val
      )
    }
    
    if (length(lambda_aucs) == 0) next
    
    # Average AUC across inner folds for each lambda
    ldf <- do.call(rbind, lambda_aucs)
    avg_auc_by_lambda <- ldf %>%
      group_by(lambda) %>%
      summarise(mean_auc = mean(auc, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_auc))
    
    top_lambda <- avg_auc_by_lambda$lambda[1]
    top_auc    <- avg_auc_by_lambda$mean_auc[1]
    
    if (!is.na(top_auc) && top_auc > best_auc) {
      best_auc    <- top_auc
      best_alpha  <- alpha_val
      best_lambda <- top_lambda
    }
  }
  
  # Fallback if tuning failed
  if (is.null(best_lambda)) {
    best_alpha  <- 0.5
    best_lambda <- 0.01
  }
  
  list(best_alpha  = best_alpha,
       best_lambda = best_lambda,
       best_inner_auc = best_auc)
}

###############################################################################
# SECTION 4: OUTER LOO LOOP — OUT-OF-FOLD PREDICTIONS
###############################################################################

cat("--- Section 4: Outer LOO loop ---\n")
cat("Total samples:", n_total, "| T2D:", n_t2d, "| Normal:", n_norm, "\n")
cat("Running", n_total, "outer LOO iterations...\n\n")

oof_probs  <- numeric(n_total)   # out-of-fold predicted probabilities
oof_alphas  <- numeric(n_total)  # best alpha selected per fold
oof_lambdas <- numeric(n_total)  # best lambda selected per fold

pb_step <- max(1, floor(n_total / 10))  # progress every 10%

for (i in seq_len(n_total)) {
  
  if (i %% pb_step == 0 || i == 1 || i == n_total) {
    cat(sprintf("  Outer fold %d / %d (T2D sample: %s)\n",
                i, n_total, ifelse(y[i] == 1, "YES", "no")))
  }
  
  # Split: test = sample i, train = all others
  test_idx  <- i
  train_idx <- setdiff(seq_len(n_total), i)
  
  X_train <- MEs_binary[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  X_test  <- MEs_binary[test_idx,  , drop = FALSE]
  
  # --- Step A: ME filtering (on training fold only) ---
  keep_cols_cor <- filter_MEs_by_correlation(X_train, y_train,
                                             cor_threshold = 0.10)
  X_train_f <- X_train[, keep_cols_cor, drop = FALSE]
  X_test_f  <- X_test[,  keep_cols_cor, drop = FALSE]
  
  # --- Step B: Remove collinear MEs (on training fold only) ---
  keep_cols_col <- remove_collinear_MEs(X_train_f, cor_cutoff = 0.85)
  X_train_f <- X_train_f[, keep_cols_col, drop = FALSE]
  X_test_f  <- X_test_f[,  keep_cols_col, drop = FALSE]
  
  # Ensure at least 1 column
  if (ncol(X_train_f) == 0) {
    X_train_f <- X_train
    X_test_f  <- X_test
  }
  
  # --- Step C: Inner CV to tune alpha + lambda ---
  tuned <- inner_cv_tune(
    X_tr        = X_train_f,
    y_tr        = y_train,
    alpha_grid  = c(0, 0.25, 0.5, 0.75, 1.0),
    k_inner     = 10,
    smote_k     = 5,
    smote_ratio = 0.5,
    seed        = i * 7  
  )
  
  best_alpha  <- tuned$best_alpha
  best_lambda <- tuned$best_lambda
  
  oof_alphas[i]  <- best_alpha
  oof_lambdas[i] <- best_lambda
  
  # --- Step D: Refit final model on FULL training fold ---
  # Scale: fit on full training fold
  tr_mean <- colMeans(X_train_f)
  tr_sd   <- apply(X_train_f, 2, sd)
  tr_sd[tr_sd == 0] <- 1
  
  X_train_sc <- scale(X_train_f, center = tr_mean, scale = tr_sd)
  X_test_sc  <- scale(X_test_f,  center = tr_mean, scale = tr_sd)
  
  # SMOTE on full training fold
  smoted_full <- apply_smote(X_train_sc, y_train,
                             k = 5, over_ratio = 0.5)
  X_tr_sm <- smoted_full$X
  y_tr_sm <- smoted_full$y
  
  # Class weights on SMOTE-augmented training set
  w_tr <- compute_class_weights(y_tr_sm)
  
  # Fit final model with best (alpha, lambda)
  final_fit <- tryCatch(
    suppressWarnings(
      glmnet(x           = X_tr_sm,
             y           = y_tr_sm,
             family      = "binomial",
             alpha       = best_alpha,
             lambda      = best_lambda,
             weights     = w_tr,
             standardize = FALSE)
    ),
    error = function(e) {
      # Fallback: ridge with no SMOTE
      glmnet(x           = X_train_sc,
             y           = y_train,
             family      = "binomial",
             alpha       = 0,
             lambda      = 0.01,
             standardize = FALSE)
    }
  )
  
  # --- Step E: Predict on held-out test sample ---
  oof_probs[i] <- as.numeric(
    predict(final_fit, X_test_sc, type = "response",
            s = best_lambda)
  )
}

cat("\nOuter LOO complete.\n\n")

###############################################################################
# SECTION 5: PERFORMANCE METRICS
###############################################################################

cat("--- Section 5: Computing performance metrics ---\n\n")

# ROC-AUC
roc_obj  <- pROC::roc(y, oof_probs,
                      levels = c(0, 1), direction = "<", quiet = TRUE)
roc_auc  <- as.numeric(pROC::auc(roc_obj))
roc_ci   <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95))

# PR-AUC (Precision-Recall)
pr_obj   <- PRROC::pr.curve(
  scores.class0 = oof_probs[y == 1],  # T2D probabilities
  scores.class1 = oof_probs[y == 0],  # Normal probabilities
  curve = TRUE
)
pr_auc   <- pr_obj$auc.integral

# Baseline PR-AUC = prevalence of positive class
baseline_pr_auc <- mean(y)

# Optimal threshold by Youden's J
coords_df  <- pROC::coords(roc_obj, "best", ret = c("threshold",
                                                    "sensitivity",
                                                    "specificity"))
best_thresh <- as.numeric(coords_df["threshold"])
pred_class  <- ifelse(oof_probs >= best_thresh, 1L, 0L)
conf_mat    <- table(Predicted = pred_class, Observed = y)

cat(sprintf("ROC-AUC:  %.3f (95%% CI: %.3f – %.3f)\n",
            roc_auc, roc_ci[1], roc_ci[3]))
cat(sprintf("PR-AUC:   %.3f (baseline = %.3f for random classifier)\n",
            pr_auc, baseline_pr_auc))
cat(sprintf("Best threshold (Youden's J): %.3f\n\n", best_thresh))
cat("Confusion matrix at best threshold:\n")
print(conf_mat)
cat("\n")

# Sensitivity and Specificity
TP <- conf_mat["1", "1"]
TN <- conf_mat["0", "0"]
FP <- conf_mat["1", "0"]
FN <- conf_mat["0", "1"]
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
cat(sprintf("Sensitivity (T2D recall): %.3f\n", sensitivity))
cat(sprintf("Specificity (Normal recall): %.3f\n\n", specificity))

###############################################################################
# SECTION 6: PERMUTATION TEST
###############################################################################

cat("--- Section 6: Permutation test (100 shuffles) ---\n")
cat("This tests whether the observed AUC is better than chance.\n")
cat("Null hypothesis: AUC achievable by randomly shuffling labels.\n\n")

n_perms <- 100
perm_aucs <- numeric(n_perms)

for (perm in seq_len(n_perms)) {
  
  set.seed(perm * 13)
  y_perm <- sample(y)  
  
  perm_probs <- numeric(n_total)
  
  for (i in seq_len(n_total)) {
    
    train_idx <- setdiff(seq_len(n_total), i)
    
    X_train_p <- MEs_binary[train_idx, , drop = FALSE]
    y_train_p <- y_perm[train_idx]
    X_test_p  <- MEs_binary[i,         , drop = FALSE]
    
    
    keep_cor <- filter_MEs_by_correlation(X_train_p, y_train_p, 0.10)
    X_tr_p   <- X_train_p[, keep_cor, drop = FALSE]
    X_te_p   <- X_test_p[,  keep_cor, drop = FALSE]
    
    tr_mean <- colMeans(X_tr_p)
    tr_sd   <- apply(X_tr_p, 2, sd); tr_sd[tr_sd == 0] <- 1
    X_tr_p_sc <- scale(X_tr_p, center = tr_mean, scale = tr_sd)
    X_te_p_sc <- scale(X_te_p, center = tr_mean, scale = tr_sd)
    
    w_p <- compute_class_weights(y_train_p)
    
    fit_p <- tryCatch(
      suppressWarnings(
        glmnet(X_tr_p_sc, y_train_p,
               family = "binomial", alpha = 0.5,
               lambda = 0.01, weights = w_p,
               standardize = FALSE)
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_p)) {
      perm_probs[i] <- 0.5
    } else {
      perm_probs[i] <- as.numeric(
        predict(fit_p, X_te_p_sc, type = "response", s = 0.01)
      )
    }
  }
  
  perm_aucs[perm] <- tryCatch(
    as.numeric(pROC::auc(pROC::roc(y_perm, perm_probs,
                                   levels = c(0, 1), direction = "<",
                                   quiet = TRUE))),
    error = function(e) 0.5
  )
  
  if (perm %% 20 == 0) cat(sprintf("  Permutation %d / %d done\n",
                                   perm, n_perms))
}

perm_p_value <- mean(perm_aucs >= roc_auc)

cat(sprintf("\nPermutation test results:\n"))
cat(sprintf("  Observed AUC:  %.3f\n", roc_auc))
cat(sprintf("  Null AUC mean: %.3f (SD: %.3f)\n",
            mean(perm_aucs), sd(perm_aucs)))
cat(sprintf("  Permutation p-value: %.4f\n", perm_p_value))
if (perm_p_value < 0.05) {
  cat("  ✓ Signal is statistically significant (p < 0.05)\n\n")
} else {
  cat("  ✗ Signal is NOT significant — model may not beat chance\n\n")
}

###############################################################################
# SECTION 7: HYPERPARAMETER SUMMARY
###############################################################################

cat("--- Section 7: Hyperparameter summary across LOO folds ---\n\n")

alpha_table <- table(oof_alphas)
cat("Alpha selected per fold:\n")
print(alpha_table)
cat(sprintf("Most common alpha: %.2f\n\n", as.numeric(names(which.max(alpha_table)))))
cat(sprintf("Lambda range: %.5f – %.5f\n",
            min(oof_lambdas), max(oof_lambdas)))
cat(sprintf("Lambda median: %.5f\n\n", median(oof_lambdas)))

###############################################################################
# SECTION 8: SAVE RESULTS
###############################################################################

cat("--- Section 8: Saving results ---\n")

dir.create("ml_smote", showWarnings = FALSE)

# Save all objects
save(
  oof_probs, y,
  roc_auc, roc_ci, pr_auc, baseline_pr_auc,
  roc_obj, pr_obj,
  conf_mat, sensitivity, specificity, best_thresh,
  perm_aucs, perm_p_value,
  oof_alphas, oof_lambdas,
  file = "ml_smote/step5_smote_results.RData"
)

# Summary CSV
summary_df <- data.frame(
  Metric         = c("ROC_AUC", "ROC_AUC_CI_lower", "ROC_AUC_CI_upper",
                     "PR_AUC", "Baseline_PR_AUC",
                     "Sensitivity", "Specificity",
                     "Best_Threshold",
                     "Permutation_p_value",
                     "N_total", "N_T2D", "N_Normal"),
  Value          = c(round(roc_auc, 4), round(roc_ci[1], 4), round(roc_ci[3], 4),
                     round(pr_auc, 4), round(baseline_pr_auc, 4),
                     round(sensitivity, 4), round(specificity, 4),
                     round(best_thresh, 4),
                     round(perm_p_value, 4),
                     n_total, n_t2d, n_norm)
)
write.csv(summary_df, "ml_smote/performance_summary.csv", row.names = FALSE)

# OOF predictions CSV
oof_df <- data.frame(
  Sample_ID    = rownames(MEs_binary),
  True_Label   = y,
  Disease      = ifelse(y == 1, "T2D", "Normal"),
  Pred_Prob    = round(oof_probs, 4),
  Pred_Class   = pred_class,
  Best_Alpha   = oof_alphas,
  Best_Lambda  = round(oof_lambdas, 6)
)
write.csv(oof_df, "ml_smote/oof_predictions.csv", row.names = FALSE)

cat("Saved: ml_smote/step5_smote_results.RData\n")
cat("Saved: ml_smote/performance_summary.csv\n")
cat("Saved: ml_smote/oof_predictions.csv\n\n")

###############################################################################
# SECTION 9: PLOTS (PDF)
###############################################################################

cat("--- Section 9: Generating plots ---\n")

pdf("ml_smote/step5_diagnostic_plots.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1)

# --- Plot 1: ROC Curve ---
plot(roc_obj,
     col  = "#E74C3C", lwd = 2.5,
     main = sprintf("ROC Curve (LOO-CV)\nAUC = %.3f [95%% CI: %.3f – %.3f]",
                    roc_auc, roc_ci[1], roc_ci[3]),
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)",
     legacy.axes = TRUE)
abline(a = 0, b = 1, lty = 2, col = "grey60")
legend("bottomright",
       legend = c(sprintf("Model (AUC = %.3f)", roc_auc), "Random"),
       col    = c("#E74C3C", "grey60"),
       lwd    = c(2.5, 1), lty = c(1, 2), bty = "n")

# --- Plot 2: PR Curve ---
plot(pr_obj,
     col  = "#2980B9", lwd = 2.5,
     main = sprintf("Precision-Recall Curve\nPR-AUC = %.3f (baseline = %.3f)",
                    pr_auc, baseline_pr_auc),
     xlab = "Recall (Sensitivity)",
     ylab = "Precision (PPV)",
     auc.main = FALSE)
abline(h = baseline_pr_auc, lty = 2, col = "grey60")
legend("topright",
       legend = c(sprintf("Model (PR-AUC = %.3f)", pr_auc),
                  sprintf("Baseline = %.3f", baseline_pr_auc)),
       col    = c("#2980B9", "grey60"),
       lwd    = c(2.5, 1), lty = c(1, 2), bty = "n")

# --- Plot 3: Permutation Test ---
hist(perm_aucs,
     breaks = 20,
     col    = "#BDC3C7",
     border = "white",
     main   = sprintf("Permutation Test (n=%d shuffles)\np = %.4f",
                      n_perms, perm_p_value),
     xlab   = "AUC under null (shuffled labels)",
     ylab   = "Frequency")
abline(v = roc_auc, col = "#E74C3C", lwd = 2.5, lty = 1)
abline(v = mean(perm_aucs), col = "#7F8C8D", lwd = 1.5, lty = 2)
legend("topright",
       legend = c(sprintf("Observed AUC = %.3f", roc_auc),
                  sprintf("Null mean = %.3f", mean(perm_aucs))),
       col    = c("#E74C3C", "#7F8C8D"),
       lwd    = c(2.5, 1.5), lty = c(1, 2), bty = "n")

# --- Plot 4: Predicted Probabilities by True Class ---
t2d_probs  <- oof_probs[y == 1]
norm_probs <- oof_probs[y == 0]

boxplot(list(Normal = norm_probs, T2D = t2d_probs),
        col    = c("#AED6F1", "#F1948A"),
        border = c("#2980B9", "#E74C3C"),
        main   = "Predicted T2D Probability by True Class\n(Out-of-Fold)",
        ylab   = "Predicted Probability (T2D)",
        xlab   = "True Class",
        ylim   = c(0, 1),
        lwd    = 1.5)
abline(h = best_thresh, lty = 2, col = "grey40")
text(x = 2.4, y = best_thresh + 0.02,
     labels = sprintf("Threshold = %.2f", best_thresh),
     col = "grey30", cex = 0.85)
# Overlay individual points
stripchart(list(Normal = norm_probs, T2D = t2d_probs),
           vertical = TRUE, method = "jitter",
           jitter   = 0.08, pch = 16, cex = 0.7,
           col      = c("#2980B950", "#E74C3C50"),
           add      = TRUE)

dev.off()
cat("Saved: ml_smote/step5_diagnostic_plots.pdf\n\n")

###############################################################################
# FINAL SUMMARY
###############################################################################

cat("=======================================================================\n")
cat("  FINAL RESULTS SUMMARY\n")
cat("=======================================================================\n")
cat(sprintf("  Dataset:      GSE50244 — T2D vs Normal (IGT excluded)\n"))
cat(sprintf("  Samples:      %d total (%d T2D, %d Normal)\n",
            n_total, n_t2d, n_norm))
cat(sprintf("  Outer CV:     Leave-One-Out (%d folds)\n", n_total))
cat(sprintf("  Inner CV:     Stratified 10-fold\n"))
cat(sprintf("  Imbalance:    SMOTE (K=5, ratio=0.5) + class weights\n"))
cat(sprintf("  Tuning:       alpha ∈ {0, 0.25, 0.5, 0.75, 1.0} + lambda\n\n"))
cat(sprintf("  ROC-AUC:      %.3f [95%% CI: %.3f – %.3f]\n",
            roc_auc, roc_ci[1], roc_ci[3]))
cat(sprintf("  PR-AUC:       %.3f (baseline = %.3f)\n",
            pr_auc, baseline_pr_auc))
cat(sprintf("  Sensitivity:  %.3f\n", sensitivity))
cat(sprintf("  Specificity:  %.3f\n", specificity))
cat(sprintf("  Perm p-value: %.4f\n\n", perm_p_value))
cat("  Output files in: ml_smote/\n")
cat("    step5_smote_results.RData\n")
cat("    performance_summary.csv\n")
cat("    oof_predictions.csv\n")
cat("    step5_diagnostic_plots.pdf\n")
cat("=======================================================================\n")
cat("=== STEP 5 COMPLETED SUCCESSFULLY ===\n")