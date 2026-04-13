###############################################################################
# BLOOD MODEL BUILDING — STEP 6
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(SHAPforxgboost)
  library(themis)
})

ROOT    <- "C:/Users/talha/OneDrive/Desktop/gene expression work"
BM_DIR  <- file.path(ROOT, "blood model building")
OBJ_DIR <- file.path(BM_DIR, "model_input")
RES_DIR <- file.path(BM_DIR, "model_results")
OUT_DIR <- file.path(BM_DIR, "step6_final")
dir.create(OUT_DIR,                        showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "plots"),    showWarnings=FALSE)

cat("=======================================================================\n")
cat("  BLOOD MODEL — STEP 6: FINAL XGBoost EVALUATION + SHAP\n")
cat("=======================================================================\n\n")

###############################################################################
# SECTION 1: LOAD DATA + OOF PREDICTIONS
###############################################################################

cat("--- Section 1: Loading data ---\n")

obj    <- readRDS(file.path(OBJ_DIR, "blood_model_step1_object.rds"))
X      <- as.matrix(obj$X)
y      <- obj$y
y_int  <- ifelse(y == "T2D", 1L, 0L)
n      <- length(y)

oof_df  <- fread(file.path(RES_DIR, "oof_predictions.csv"))
oof_xgb <- oof_df$P_T2D_XGB

cat(sprintf("Samples: %d (%d T2D, %d Control)\n",
            n, sum(y_int), sum(y_int==0)))
cat(sprintf("Features: %d MEs\n\n", ncol(X)))

###############################################################################
# SECTION 2: LOAD BEST HYPERPARAMETERS FROM LOO
###############################################################################

cat("--- Section 2: Loading best hyperparameters ---\n")

hp_df <- fread(file.path(RES_DIR, "hyperparameters_per_fold.csv"))

best_depth   <- as.integer(names(which.max(table(hp_df$XGB_depth))))
best_eta     <- as.numeric(names(which.max(table(hp_df$XGB_eta))))
best_nrounds <- as.integer(names(which.max(table(hp_df$XGB_nrounds))))

cat(sprintf("Median hyperparams from LOO:\n"))
cat(sprintf("  max_depth = %d\n", best_depth))
cat(sprintf("  eta       = %.3f\n", best_eta))
cat(sprintf("  nrounds   = %d\n\n", best_nrounds))

###############################################################################
# SECTION 3: OOF PERFORMANCE METRICS
###############################################################################

cat("--- Section 3: OOF performance metrics ---\n")

roc_obj <- pROC::roc(y_int, oof_xgb,
                     levels=c(0,1), direction="<", quiet=TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))
ci_val  <- as.numeric(pROC::ci.auc(roc_obj, conf.level=0.95))

pr_obj  <- PRROC::pr.curve(
  scores.class0 = oof_xgb[y_int==1],
  scores.class1 = oof_xgb[y_int==0],
  curve = TRUE)
pr_auc  <- pr_obj$auc.integral

coords  <- pROC::coords(roc_obj, "best",
                        ret=c("threshold","sensitivity","specificity"))
thresh  <- if(is.finite(coords$threshold[1])) coords$threshold[1] else 0.5

pred_c  <- ifelse(oof_xgb >= thresh, 1L, 0L)
cm      <- table(Predicted=pred_c, Observed=y_int)

TP <- if("1" %in% rownames(cm) & "1" %in% colnames(cm)) cm["1","1"] else 0L
TN <- if("0" %in% rownames(cm) & "0" %in% colnames(cm)) cm["0","0"] else 0L
FP <- if("1" %in% rownames(cm) & "0" %in% colnames(cm)) cm["1","0"] else 0L
FN <- if("0" %in% rownames(cm) & "1" %in% colnames(cm)) cm["0","1"] else 0L
sens <- TP / max(TP+FN, 1)
spec <- TN / max(TN+FP, 1)
ppv  <- TP / max(TP+FP, 1)
npv  <- TN / max(TN+FN, 1)
f1   <- 2*TP / max(2*TP+FP+FN, 1)

cat(sprintf("XGBoost LOO-CV Performance:\n"))
cat(sprintf("  AUC:         %.3f [%.3f-%.3f]\n", auc_val, ci_val[1], ci_val[3]))
cat(sprintf("  PR-AUC:      %.3f\n", pr_auc))
cat(sprintf("  Threshold:   %.4f\n", thresh))
cat(sprintf("  Sensitivity: %.3f\n", sens))
cat(sprintf("  Specificity: %.3f\n", spec))
cat(sprintf("  PPV:         %.3f\n", ppv))
cat(sprintf("  NPV:         %.3f\n", npv))
cat(sprintf("  F1 Score:    %.3f\n", f1))
cat(sprintf("  TP=%d TN=%d FP=%d FN=%d\n\n", TP, TN, FP, FN))

###############################################################################
# SECTION 4: CONFUSION MATRIX PLOT
###############################################################################

cat("--- Section 4: Confusion matrix plot ---\n")

cm_df <- data.frame(
  Predicted = factor(c("Control","Control","T2D","T2D"),
                     levels=c("Control","T2D")),
  Observed  = factor(c("Control","T2D","Control","T2D"),
                     levels=c("Control","T2D")),
  Count     = c(as.integer(TN), as.integer(FN),
                as.integer(FP), as.integer(TP)),
  Label     = c(
    sprintf("TN\n%d", as.integer(TN)),
    sprintf("FN\n%d", as.integer(FN)),
    sprintf("FP\n%d", as.integer(FP)),
    sprintf("TP\n%d", as.integer(TP))
  )
)

p_cm <- ggplot(cm_df, aes(x=Observed, y=Predicted, fill=Count)) +
  geom_tile(color="white", linewidth=1.5) +
  geom_text(aes(label=Label), size=7, fontface="bold", color="white") +
  scale_fill_gradient(low="#AED6F1", high="#1A5276") +
  scale_x_discrete(position="top") +
  theme_minimal(base_size=14) +
  theme(axis.text=element_text(size=13, face="bold"),
        axis.title=element_text(size=13, face="bold"),
        panel.grid=element_blank(),
        legend.position="right") +
  labs(title="XGBoost — Confusion Matrix (LOO-CV)",
       subtitle=sprintf("Threshold=%.3f | Sens=%.3f | Spec=%.3f | AUC=%.3f",
                        thresh, sens, spec, auc_val),
       x="Observed", y="Predicted", fill="Count")

ggsave(file.path(OUT_DIR, "plots", "step6_confusion_matrix.png"),
       p_cm, width=7, height=6, dpi=150)
cat("Saved: confusion matrix\n\n")

###############################################################################
# SECTION 5: ROC CURVE PLOT
###############################################################################

cat("--- Section 5: ROC curve plot ---\n")

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x=FPR, y=TPR)) +
  geom_line(color="#E74C3C", linewidth=1.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              color="grey60", linewidth=0.8) +
  geom_point(aes(x=1-spec, y=sens),
             color="#2980B9", size=4, shape=16) +
  annotate("text", x=1-spec+0.03, y=sens-0.04,
           label=sprintf("Best threshold\nSens=%.2f, Spec=%.2f",
                         sens, spec),
           size=3.5, color="#2980B9", hjust=0) +
  annotate("text", x=0.6, y=0.15,
           label=sprintf("AUC = %.3f\n[%.3f-%.3f]",
                         auc_val, ci_val[1], ci_val[3]),
           size=4.5, color="#E74C3C", fontface="bold") +
  theme_minimal(base_size=13) +
  labs(title="XGBoost — ROC Curve (LOO-CV)",
       x="1 - Specificity (FPR)", y="Sensitivity (TPR)")

ggsave(file.path(OUT_DIR, "plots", "step6_roc_curve.png"),
       p_roc, width=7, height=6, dpi=150)
cat("Saved: ROC curve\n\n")

###############################################################################
# SECTION 6: TRAIN FINAL FROZEN MODEL ON ALL 116 SAMPLES
###############################################################################

cat("--- Section 6: Training frozen model on all 116 samples ---\n")

# Scale (fit on all 116)
X_means <- colMeans(X)
X_sds   <- apply(X, 2, sd); X_sds[X_sds==0] <- 1
X_sc    <- scale(X, X_means, X_sds)

# SMOTE on full data
sm_full <- tryCatch({
  df <- as.data.frame(X_sc)
  df$.outcome <- factor(ifelse(y_int==1,"T2D","Control"),
                        levels=c("Control","T2D"))
  out <- themis::smote(df, var=".outcome", k=3, over_ratio=0.5)
  y_out <- ifelse(out$.outcome=="T2D", 1L, 0L)
  X_out <- as.matrix(out[, colnames(X_sc), drop=FALSE])
  mode(X_out) <- "numeric"
  list(X=X_out, y=y_out)
}, error=function(e) list(X=X_sc, y=y_int))

scale_pos <- sum(sm_full$y==0) / sum(sm_full$y==1)
dtrain    <- xgboost::xgb.DMatrix(sm_full$X, label=sm_full$y)

final_xgb <- xgboost::xgb.train(
  params = list(
    objective         = "binary:logistic",
    eval_metric       = "auc",
    max_depth         = best_depth,
    eta               = best_eta,
    scale_pos_weight  = scale_pos,
    subsample         = 0.8,
    colsample_bytree  = 0.8
  ),
  data    = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

cat("Frozen model trained.\n\n")

