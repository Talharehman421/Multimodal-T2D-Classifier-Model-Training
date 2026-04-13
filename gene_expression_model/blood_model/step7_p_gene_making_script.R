###############################################################################
# BLOOD GENE EXPRESSION MODEL — p_gene EXPORT FOR LATE FUSION
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(ggplot2)
  library(PRROC)
})

set.seed(2025)

cat("=======================================================================\n")
cat("  GENE MODEL — p_gene EXPORT WITH PLATT CALIBRATION\n")
cat("=======================================================================\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 0: CONFIGURATION — UPDATE THESE PATHS
# ─────────────────────────────────────────────────────────────────────────────

ROOT    <- "C:/Users/talha/OneDrive/Desktop/gene expression work"
RES_DIR <- file.path(ROOT, "blood model building", "model_results")
META_DIR <- file.path(ROOT, "blood model building","meta")   # folder with master CSV
OUT_DIR <- file.path(ROOT, "blood model building", "p_gene_output")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OOF_PATH    <- file.path(RES_DIR,  "oof_predictions.csv")
MASTER_PATH <- file.path(META_DIR, "blood_samples_master.csv")

stopifnot(file.exists(OOF_PATH))
stopifnot(file.exists(MASTER_PATH))

###############################################################################
# SECTION 1: LOAD DATA
###############################################################################

cat("--- Section 1: Loading OOF predictions and master file ---\n")

oof_df     <- fread(OOF_PATH)
master_df  <- fread(MASTER_PATH)

cat(sprintf("OOF predictions loaded:  %d rows x %d cols\n",
            nrow(oof_df), ncol(oof_df)))
cat(sprintf("Master file loaded:      %d rows x %d cols\n",
            nrow(master_df), ncol(master_df)))

# Verify row counts match
stopifnot(nrow(oof_df) == nrow(master_df))
cat(sprintf("Row count match:         %d samples ✓\n\n", nrow(oof_df)))

###############################################################################
# SECTION 2: ALIGN AND VERIFY LABELS
###############################################################################

cat("--- Section 2: Aligning labels and verifying consistency ---\n")


y_label <- oof_df$True_Label          
y_int   <- oof_df$True_Int            

# Cross-check against master file Disease_Status column
master_labels <- master_df$Disease_Status
label_match   <- all(y_label == master_labels)

if (!label_match) {
  mismatches <- which(y_label != master_labels)
  cat(sprintf("WARNING: %d label mismatches found at rows: %s\n",
              length(mismatches), paste(mismatches, collapse=",")))
  stop("Label mismatch between oof_predictions.csv and master file. Cannot proceed.")
} else {
  cat(sprintf("Label alignment:         PASS — all %d labels match ✓\n", length(y_label)))
}

n_t2d  <- sum(y_int == 1)
n_ctrl <- sum(y_int == 0)
cat(sprintf("T2D samples:             %d\n", n_t2d))
cat(sprintf("Control samples:         %d\n\n", n_ctrl))


p_raw <- oof_df$P_T2D_XGB

cat(sprintf("p_gene_raw summary:\n"))
cat(sprintf("  T2D pool   — mean=%.4f  sd=%.4f  min=%.4f  max=%.4f\n",
            mean(p_raw[y_int==1]), sd(p_raw[y_int==1]),
            min(p_raw[y_int==1]), max(p_raw[y_int==1])))
cat(sprintf("  Control pool — mean=%.4f  sd=%.4f  min=%.4f  max=%.4f\n\n",
            mean(p_raw[y_int==0]), sd(p_raw[y_int==0]),
            min(p_raw[y_int==0]), max(p_raw[y_int==0])))

###############################################################################
# SECTION 3: PLATT SCALING CALIBRATION
###############################################################################
###############################################################################

cat("--- Section 3: Platt scaling calibration ---\n")

# Fit logistic regression: logit(p_calibrated) = a + b * p_raw
platt_df  <- data.frame(p_raw = p_raw, y = y_int)
platt_fit <- glm(y ~ p_raw, data = platt_df, family = binomial(link = "logit"))

cat("Platt scaling model fitted.\n")
cat(sprintf("  Intercept (a): %.4f\n", coef(platt_fit)[1]))
cat(sprintf("  Slope     (b): %.4f\n", coef(platt_fit)[2]))

# Sanity check: slope should be positive (higher raw prob → higher calibrated prob)
if (coef(platt_fit)[2] <= 0) {
  warning("Platt slope is non-positive — calibration may be unreliable. Check OOF probabilities.")
} else {
  cat("  Slope sign:    POSITIVE ✓ (calibration direction is correct)\n")
}

# Apply calibration to get p_gene_calibrated
p_calibrated <- predict(platt_fit, newdata = platt_df, type = "response")

cat(sprintf("\np_gene_calibrated summary:\n"))
cat(sprintf("  T2D pool   — mean=%.4f  sd=%.4f  min=%.4f  max=%.4f\n",
            mean(p_calibrated[y_int==1]), sd(p_calibrated[y_int==1]),
            min(p_calibrated[y_int==1]), max(p_calibrated[y_int==1])))
cat(sprintf("  Control pool — mean=%.4f  sd=%.4f  min=%.4f  max=%.4f\n",
            mean(p_calibrated[y_int==0]), sd(p_calibrated[y_int==0]),
            min(p_calibrated[y_int==0]), max(p_calibrated[y_int==0])))
cat(sprintf("  Overall mean: %.4f  (actual positive rate: %.4f)\n\n",
            mean(p_calibrated), mean(y_int)))

###############################################################################
# SECTION 4: CALIBRATION QUALITY CHECK
###############################################################################

cat("--- Section 4: Calibration quality check ---\n")

# AUC before and after calibration — should be identical
roc_raw <- pROC::roc(y_int, p_raw,
                     levels=c(0,1), direction="<", quiet=TRUE)
roc_cal <- pROC::roc(y_int, p_calibrated,
                     levels=c(0,1), direction="<", quiet=TRUE)

auc_raw <- as.numeric(pROC::auc(roc_raw))
auc_cal <- as.numeric(pROC::auc(roc_cal))

cat(sprintf("AUC (raw):        %.4f\n", auc_raw))
cat(sprintf("AUC (calibrated): %.4f\n", auc_cal))
cat(sprintf("AUC difference:   %.6f  (should be ~0)\n\n", abs(auc_raw - auc_cal)))

# Brier score — lower is better, calibrated should be <= raw
brier_raw <- mean((p_raw - y_int)^2)
brier_cal <- mean((p_calibrated - y_int)^2)
cat(sprintf("Brier score (raw):        %.4f\n", brier_raw))
cat(sprintf("Brier score (calibrated): %.4f\n", brier_cal))
if (brier_cal <= brier_raw) {
  cat("Brier improvement:        YES ✓\n\n")
} else {
  cat("Brier note: calibrated Brier slightly higher — acceptable for Platt on small n\n\n")
}

# Calibration in bins — expected vs observed
n_bins   <- 5   
bin_cuts <- quantile(p_calibrated, probs = seq(0, 1, length.out = n_bins + 1))
bin_id   <- cut(p_calibrated, breaks = bin_cuts, include.lowest = TRUE,
                labels = FALSE)

cat("Calibration check (5 quantile bins):\n")
cat(sprintf("  %-12s  %-14s  %-14s  %-10s\n",
            "Bin", "Mean predicted", "Observed rate", "N samples"))
for (b in 1:n_bins) {
  idx <- which(bin_id == b)
  if (length(idx) == 0) next
  cat(sprintf("  Bin %-8d  %-14.4f  %-14.4f  %-10d\n",
              b, mean(p_calibrated[idx]), mean(y_int[idx]), length(idx)))
}
cat("\n")

###############################################################################
# SECTION 5: SAVE p_gene.csv
###############################################################################

cat("--- Section 5: Saving p_gene.csv ---\n")

p_gene_df <- data.frame(
  sample_index    = oof_df$Sample_Index,
  sample_id       = master_df$Sample_ID,
  true_label      = y_label,           # "T2D" or "Control"
  true_int        = y_int,             # 1 or 0
  p_gene_raw      = round(p_raw,         6),
  p_gene          = round(p_calibrated,  6)   # use this column in fusion
)

p_gene_path <- file.path(OUT_DIR, "p_gene.csv")
fwrite(p_gene_df, p_gene_path)

cat(sprintf("Saved: %s\n", p_gene_path))
cat(sprintf("  Rows:    %d\n", nrow(p_gene_df)))
cat(sprintf("  Columns: %s\n\n", paste(names(p_gene_df), collapse=", ")))

# Print first 5 rows for visual check
cat("First 5 rows of p_gene.csv:\n")
print(head(p_gene_df, 5))
cat("\nLast 5 rows of p_gene.csv:\n")
print(tail(p_gene_df, 5))

###############################################################################
# SECTION 6: SAVE PLATT SCALER OBJECT
###############################################################################

cat("\n--- Section 6: Saving Platt scaler ---\n")

# Save as a list containing everything needed to apply calibration
platt_scaler <- list(
  model          = platt_fit,          
  intercept      = coef(platt_fit)[1],
  slope          = coef(platt_fit)[2],
  # metadata for documentation
  fitted_on      = "LOO-CV OOF XGBoost probabilities (n=116)",
  calibration    = "Platt scaling (logistic regression)",
  auc_raw        = auc_raw,
  auc_calibrated = auc_cal,
  brier_raw      = brier_raw,
  brier_cal      = brier_cal,
  # usage instructions stored inside the object
  usage = paste(
    "To calibrate new XGBoost predictions from the gene model:",
    "  new_p_calibrated <- predict(platt_scaler$model,",
    "                               newdata = data.frame(p_raw = new_p_raw),",
    "                               type = 'response')",
    sep = "\n"
  )
)

scaler_path <- file.path(OUT_DIR, "platt_scaler.rds")
saveRDS(platt_scaler, scaler_path)
cat(sprintf("Saved: %s\n", scaler_path))
cat(sprintf("  Intercept: %.4f\n", platt_scaler$intercept))
cat(sprintf("  Slope:     %.4f\n\n", platt_scaler$slope))

###############################################################################
# SECTION 7: DIAGNOSTIC PLOTS
###############################################################################

cat("--- Section 7: Generating diagnostic plots ---\n")

# ── Plot 1: Calibration curve 
n_cal_bins   <- 8
bin_breaks   <- quantile(p_raw, probs = seq(0, 1, length.out = n_cal_bins + 1))
bin_id_raw   <- cut(p_raw, breaks = bin_breaks,
                    include.lowest = TRUE, labels = FALSE)

cal_plot_df <- do.call(rbind, lapply(1:n_cal_bins, function(b) {
  idx <- which(bin_id_raw == b)
  if (length(idx) < 3) return(NULL)
  data.frame(
    mean_raw        = mean(p_raw[idx]),
    mean_calibrated = mean(p_calibrated[idx]),
    observed_rate   = mean(y_int[idx]),
    n               = length(idx)
  )
}))
cal_plot_df <- cal_plot_df[!is.null(cal_plot_df), ]

p_cal <- ggplot(cal_plot_df) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              color="grey60", linewidth=0.8) +
  geom_line(aes(x=mean_raw, y=observed_rate),
            color="#E74C3C", linewidth=1.2) +
  geom_point(aes(x=mean_raw, y=observed_rate),
             color="#E74C3C", size=3) +
  geom_line(aes(x=mean_calibrated, y=observed_rate),
            color="#27AE60", linewidth=1.2) +
  geom_point(aes(x=mean_calibrated, y=observed_rate),
             color="#27AE60", size=3) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  theme_minimal(base_size=13) +
  labs(
    title    = "Calibration Curve — Gene Model (XGBoost LOO-CV)",
    subtitle = sprintf("Red = raw (AUC=%.3f, Brier=%.4f)  |  Green = Platt calibrated (Brier=%.4f)",
                       auc_raw, brier_raw, brier_cal),
    x        = "Mean Predicted Probability",
    y        = "Observed Fraction of T2D Cases"
  ) +
  annotate("text", x=0.75, y=0.15,
           label="Perfect calibration\n(dashed diagonal)",
           color="grey50", size=3.5)

ggsave(file.path(OUT_DIR, "p_gene_calibration_plot.png"),
       p_cal, width=7, height=6, dpi=150)
cat("Saved: p_gene_calibration_plot.png\n")

# ── Plot 2: Probability distributions by class (raw vs calibrated) ──────────
dist_df <- data.frame(
  probability = c(p_raw, p_calibrated),
  type        = rep(c("Raw", "Platt Calibrated"), each=116),
  class       = rep(ifelse(y_int==1, "T2D", "Control"), 2)
)
dist_df$type  <- factor(dist_df$type,  levels=c("Raw","Platt Calibrated"))
dist_df$class <- factor(dist_df$class, levels=c("Control","T2D"))

p_dist <- ggplot(dist_df, aes(x=probability, fill=class)) +
  geom_histogram(binwidth=0.05, alpha=0.7, position="identity",
                 color="white", linewidth=0.3) +
  scale_fill_manual(values=c(Control="#2980B9", T2D="#E74C3C")) +
  facet_wrap(~type, ncol=2) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom") +
  labs(
    title    = "p_gene Distribution by Class — Raw vs Platt Calibrated",
    subtitle = sprintf("T2D: raw mean=%.3f → cal mean=%.3f  |  Control: raw mean=%.3f → cal mean=%.3f",
                       mean(p_raw[y_int==1]),        mean(p_calibrated[y_int==1]),
                       mean(p_raw[y_int==0]),        mean(p_calibrated[y_int==0])),
    x        = "Predicted Probability of T2D",
    y        = "Count",
    fill     = "True Class"
  )

ggsave(file.path(OUT_DIR, "p_gene_distribution_plot.png"),
       p_dist, width=10, height=5, dpi=150)
cat("Saved: p_gene_distribution_plot.png\n\n")

###############################################################################
# SECTION 8: POOLS SUMMARY FOR LATE FUSION
###############################################################################

cat("--- Section 8: Pools summary for late fusion ---\n\n")

t2d_pool  <- p_gene_df$p_gene[p_gene_df$true_int == 1]
ctrl_pool <- p_gene_df$p_gene[p_gene_df$true_int == 0]

cat("These are the two pools Python will sample from during\n")
cat("synthetic meta-training:\n\n")
cat(sprintf("P_gene_pos (T2D pool, n=%d):\n", length(t2d_pool)))
cat(sprintf("  min=%.4f  Q1=%.4f  median=%.4f  mean=%.4f  Q3=%.4f  max=%.4f\n\n",
            min(t2d_pool),
            quantile(t2d_pool, 0.25),
            median(t2d_pool),
            mean(t2d_pool),
            quantile(t2d_pool, 0.75),
            max(t2d_pool)))

cat(sprintf("P_gene_neg (Control pool, n=%d):\n", length(ctrl_pool)))
cat(sprintf("  min=%.4f  Q1=%.4f  median=%.4f  mean=%.4f  Q3=%.4f  max=%.4f\n\n",
            min(ctrl_pool),
            quantile(ctrl_pool, 0.25),
            median(ctrl_pool),
            mean(ctrl_pool),
            quantile(ctrl_pool, 0.75),
            max(ctrl_pool)))

cat(sprintf("Pool separation (mean T2D - mean Control): %.4f\n\n",
            mean(t2d_pool) - mean(ctrl_pool)))

###############################################################################
# FINAL SUMMARY
###############################################################################

cat("=======================================================================\n")
cat("  p_gene EXPORT COMPLETE\n")
cat("=======================================================================\n")
cat(sprintf("  Samples:           %d (%d T2D, %d Control)\n",
            nrow(p_gene_df), n_t2d, n_ctrl))
cat(sprintf("  Calibration:       Platt scaling (logistic regression)\n"))
cat(sprintf("  AUC (preserved):   %.4f (raw) → %.4f (calibrated)\n",
            auc_raw, auc_cal))
cat(sprintf("  Brier score:       %.4f (raw) → %.4f (calibrated)\n",
            brier_raw, brier_cal))
cat(sprintf("\n  Outputs saved to:  %s\n", OUT_DIR))
cat(sprintf("    p_gene.csv                   — use 'p_gene' column in fusion\n"))
cat(sprintf("    platt_scaler.rds             — save for web app deployment\n"))
cat(sprintf("    p_gene_calibration_plot.png  — calibration quality check\n"))
cat(sprintf("    p_gene_distribution_plot.png — probability distributions\n"))
cat("\n  NEXT STEP: Copy p_gene.csv to your Python late fusion folder\n")
cat("  alongside p_clinical.csv and p_lifestyle.csv\n")
cat("=======================================================================\n")