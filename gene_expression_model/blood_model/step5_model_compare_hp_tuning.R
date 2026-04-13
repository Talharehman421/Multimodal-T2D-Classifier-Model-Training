###############################################################################
#  STEP 5
# Blood ME Model: Multi-Model Nested LOO-CV + Ensemble
#
# MODELS:
#   1. Elastic Net (EN)   — alpha/lambda tuned
#   2. Random Forest (RF) — mtry/ntree tuned
#   3. XGBoost (XGB)      — depth/eta/nrounds tuned
#   4. SVM RBF            — cost/gamma tuned
#   5. Ensemble           — average of EN + RF + XGB + SVM probabilities
#
# PIPELINE PER LOO FOLD (train=115, test=1):
#   Feature filter (|r|>0.10) → collinearity removal (|r|<0.85)
#   → Scale (fit on train only) → SMOTE (K=3, ratio=0.5) on train only
#   → Stratified 10-fold inner CV → tune hyperparams → refit on full train
#   → predict 1 test sample

suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(e1071)
  library(pROC)
  library(PRROC)
  library(themis)
  library(tidyverse)
  library(data.table)
  library(doParallel)
  library(parallel)
  library(foreach)
})

set.seed(2025)

cat("=======================================================================\n")
cat("  PHASE 5 STEP 2B (REVISED): Multi-Model Blood ME Nested LOO-CV\n")
cat("=======================================================================\n\n")

###############################################################################
# SECTION 0: CONFIGURATION
###############################################################################

PROJECT_ROOT <- "C:/Users/talha/OneDrive/Desktop/gene expression work"
PHASE5_DIR   <- file.path(PROJECT_ROOT, "blood model building")

INPUT_RDS <- file.path(ROOT, "blood model building", "model_input",
                       "blood_model_step1_object.rds")
OUT_DIR <- file.path(ROOT, "blood model building", "model_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


USE_PARALLEL <- TRUE
N_CORES      <- max(1, detectCores() - 1)

# Tuning grids
ALPHA_GRID   <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
LAMBDA_GRID  <- exp(seq(-6, 2, length.out = 50))
RF_MTRY_GRID <- c(2, 3, 5, 7, 10)   
RF_NTREE     <- 500
XGB_GRID <- expand.grid(
  max_depth = c(2, 3, 4),
  eta       = c(0.01, 0.05, 0.1),
  nrounds   = c(50, 100, 200)
)
SVM_COST_GRID  <- c(0.01, 0.1, 1, 10, 100)
SVM_GAMMA_GRID <- c(0.001, 0.01, 0.1, 1)

# Inner CV folds (stratified, manual)
K_INNER <- 10

# SMOTE parameters
SMOTE_K     <- 3    # K=3 safer than K=5 for n_minority~44 in inner train
SMOTE_RATIO <- 0.5  # T2D:Control = 1:2 after SMOTE (not full balance)

# Permutation test
N_PERMS <- 100

stopifnot(file.exists(INPUT_RDS))

###############################################################################
# SECTION 1: LOAD DATA
###############################################################################

cat("--- Section 1: Loading Phase 5 Step 1 object ---\n")

obj <- readRDS(INPUT_RDS)
X   <- as.matrix(obj$X)
y   <- obj$y   # factor: levels = c("Control","T2D")

stopifnot(nrow(X) == length(y))
stopifnot(all(levels(y) == c("Control", "T2D")))

# Convert to binary integer for glmnet/xgboost
y_int <- ifelse(y == "T2D", 1L, 0L)

n_total  <- length(y)
n_t2d    <- sum(y_int == 1)
n_ctrl   <- sum(y_int == 0)
n_feat   <- ncol(X)

cat(sprintf("Samples: %d total (%d T2D, %d Control)\n",
            n_total, n_t2d, n_ctrl))
cat(sprintf("Features: %d MEs\n", n_feat))
cat(sprintf("Imbalance ratio: 1:%.2f (T2D:Control)\n\n",
            n_ctrl / n_t2d))

# Adjust RF mtry grid to not exceed n_feat
RF_MTRY_GRID <- RF_MTRY_GRID[RF_MTRY_GRID <= n_feat]
if (length(RF_MTRY_GRID) == 0) RF_MTRY_GRID <- max(1, floor(sqrt(n_feat)))

###############################################################################
# SECTION 2: HELPER FUNCTIONS
###############################################################################

cat("--- Section 2: Defining helper functions ---\n\n")

# ── Feature filtering ────────────────────────────────────────────────────────
filter_features <- function(X_train, y_train, cor_thr = 0.10,
                            collin_thr = 0.85) {
  # Step 1: correlation with outcome
  cors <- apply(X_train, 2, function(col) {
    if (var(col) < 1e-10) return(0)
    cor(col, y_train, use = "complete.obs")
  })
  keep <- which(abs(cors) >= cor_thr)
  if (length(keep) == 0)
    keep <- order(abs(cors), decreasing = TRUE)[1:min(3, ncol(X_train))]
  
  X_f <- X_train[, keep, drop = FALSE]
  
  # Step 2: remove collinear features
  if (ncol(X_f) > 1) {
    cm <- cor(X_f, use = "complete.obs")
    cm[lower.tri(cm, diag = TRUE)] <- 0
    drop_set <- c()
    for (j in seq_len(nrow(cm))) {
      if (j %in% drop_set) next
      high <- which(abs(cm[j, ]) > collin_thr)
      drop_set <- union(drop_set, high)
    }
    keep2 <- setdiff(seq_len(ncol(X_f)), drop_set)
    if (length(keep2) == 0) keep2 <- 1L
    X_f <- X_f[, keep2, drop = FALSE]
  }
  
  list(X = X_f, cols = colnames(X_f))
}

# ── SMOTE wrapper ────────────────────────────────────────────────────────────
apply_smote <- function(X_train, y_train, k = 3, over_ratio = 0.5) {
  df <- as.data.frame(X_train)
  df$.outcome <- factor(ifelse(y_train == 1, "T2D", "Control"),
                        levels = c("Control", "T2D"))
  n_min <- sum(y_train == 1)
  k_use <- min(k, max(1L, n_min - 1L))
  
  out <- tryCatch(
    themis::smote(df, var = ".outcome", k = k_use,
                  over_ratio = over_ratio),
    error = function(e) df
  )
  y_out <- ifelse(out$.outcome == "T2D", 1L, 0L)
  X_out <- as.matrix(out[, colnames(X_train), drop = FALSE])
  mode(X_out) <- "numeric"
  list(X = X_out, y = y_out)
}

# ── Class weights ─────────────────────────────────────────────────────────────
make_weights <- function(y_vec) {
  n1 <- sum(y_vec == 1); n0 <- sum(y_vec == 0)
  ifelse(y_vec == 1, n0 / n1, 1.0)
}

# ── Stratified k-fold ─────────────────────────────────────────────────────────
strat_kfold <- function(y_vec, k = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx0 <- sample(which(y_vec == 0))
  idx1 <- sample(which(y_vec == 1))
  fold_id <- integer(length(y_vec))
  fold_id[idx0] <- ((seq_along(idx0) - 1) %% k) + 1
  fold_id[idx1] <- ((seq_along(idx1) - 1) %% k) + 1
  lapply(1:k, function(f) which(fold_id == f))
}

# ── Tune Elastic Net ──────────────────────────────────────────────────────────
tune_en <- function(X_tr, y_tr, alpha_grid, k = 10, seed = 42) {
  folds    <- strat_kfold(y_tr, k = k, seed = seed)
  best_auc <- -Inf; best_a <- 0.5; best_l <- 0.01
  
  for (a in alpha_grid) {
    ldf_list <- lapply(seq_along(folds), function(fi) {
      vi <- folds[[fi]]; ti <- setdiff(seq_len(nrow(X_tr)), vi)
      Xti <- X_tr[ti, , drop=F]; yti <- y_tr[ti]
      Xvi <- X_tr[vi, , drop=F]; yvi <- y_tr[vi]
      m <- colMeans(Xti); s <- apply(Xti, 2, sd); s[s==0] <- 1
      Xti_s <- scale(Xti, m, s); Xvi_s <- scale(Xvi, m, s)
      sm <- apply_smote(Xti_s, yti, k=SMOTE_K, over_ratio=SMOTE_RATIO)
      w  <- make_weights(sm$y)
      fit <- tryCatch(suppressWarnings(
        glmnet(sm$X, sm$y, family="binomial", alpha=a,
               weights=w, standardize=FALSE)),
        error=function(e) NULL)
      if (is.null(fit) || length(unique(yvi)) < 2) return(NULL)
      preds <- predict(fit, Xvi_s, type="response")
      data.frame(
        lambda = fit$lambda,
        auc = apply(preds, 2, function(p)
          tryCatch(as.numeric(pROC::auc(
            pROC::roc(yvi, p, levels=c(0,1), direction="<", quiet=TRUE))),
            error=function(e) NA_real_)),
        alpha = a
      )
    })
    ldf <- do.call(rbind, Filter(Negate(is.null), ldf_list))
    if (is.null(ldf) || nrow(ldf) == 0) next
    avg <- ldf %>% group_by(lambda) %>%
      summarise(m = mean(auc, na.rm=TRUE), .groups="drop") %>%
      arrange(desc(m))
    if (!is.na(avg$m[1]) && avg$m[1] > best_auc) {
      best_auc <- avg$m[1]; best_a <- a; best_l <- avg$lambda[1]
    }
  }
  list(alpha = best_a, lambda = best_l)
}

# ── Tune Random Forest ────────────────────────────────────────────────────────
tune_rf <- function(X_tr, y_tr, mtry_grid, ntree = 500,
                    k = 10, seed = 42) {
  folds    <- strat_kfold(y_tr, k = k, seed = seed)
  best_auc <- -Inf; best_mtry <- mtry_grid[1]
  
  for (mt in mtry_grid) {
    mt_use <- min(mt, ncol(X_tr))
    aucs <- sapply(seq_along(folds), function(fi) {
      vi <- folds[[fi]]; ti <- setdiff(seq_len(nrow(X_tr)), vi)
      yvi_fac <- factor(ifelse(y_tr[vi]==1,"T2D","Control"),
                        levels=c("Control","T2D"))
      yti_fac <- factor(ifelse(y_tr[ti]==1,"T2D","Control"),
                        levels=c("Control","T2D"))
      if (length(unique(y_tr[vi])) < 2) return(NA)
      # Class weights for ranger
      cw <- table(yti_fac)
      fit <- tryCatch(
        ranger::ranger(x=X_tr[ti,,drop=F], y=yti_fac,
                       num.trees=ntree, mtry=mt_use,
                       probability=TRUE,
                       class.weights=c(Control=1, T2D=cw["Control"]/cw["T2D"]),
                       seed=seed+fi, verbose=FALSE),
        error=function(e) NULL)
      if (is.null(fit)) return(NA)
      p <- predict(fit, X_tr[vi,,drop=F])$predictions[,"T2D"]
      tryCatch(as.numeric(pROC::auc(
        pROC::roc(y_tr[vi], p, levels=c(0,1), direction="<", quiet=TRUE))),
        error=function(e) NA)
    })
    mean_auc <- mean(aucs, na.rm=TRUE)
    if (!is.na(mean_auc) && mean_auc > best_auc) {
      best_auc <- mean_auc; best_mtry <- mt_use
    }
  }
  list(mtry = best_mtry, ntree = ntree)
}

# ── Tune XGBoost ──────────────────────────────────────────────────────────────
tune_xgb <- function(X_tr, y_tr, xgb_grid, k = 10, seed = 42) {
  folds    <- strat_kfold(y_tr, k = k, seed = seed)
  best_auc <- -Inf
  best_params <- list(max_depth=3, eta=0.05, nrounds=100)
  scale_pos <- sum(y_tr==0) / sum(y_tr==1)
  
  for (gi in seq_len(nrow(xgb_grid))) {
    md <- xgb_grid$max_depth[gi]
    et <- xgb_grid$eta[gi]
    nr <- xgb_grid$nrounds[gi]
    
    aucs <- sapply(seq_along(folds), function(fi) {
      vi <- folds[[fi]]; ti <- setdiff(seq_len(nrow(X_tr)), vi)
      if (length(unique(y_tr[vi])) < 2) return(NA)
      dtrain <- xgboost::xgb.DMatrix(X_tr[ti,,drop=F], label=y_tr[ti])
      fit <- tryCatch(suppressWarnings(
        xgboost::xgb.train(
          params = list(objective="binary:logistic",
                        eval_metric="auc",
                        max_depth=md, eta=et,
                        scale_pos_weight=scale_pos,
                        subsample=0.8, colsample_bytree=0.8),
          data = dtrain, nrounds = nr,
          verbose = 0, nthread = 1)),
        error=function(e) NULL)
      if (is.null(fit)) return(NA)
      dtest <- xgboost::xgb.DMatrix(X_tr[vi,,drop=F])
      p <- predict(fit, dtest)
      tryCatch(as.numeric(pROC::auc(
        pROC::roc(y_tr[vi], p, levels=c(0,1), direction="<", quiet=TRUE))),
        error=function(e) NA)
    })
    mean_auc <- mean(aucs, na.rm=TRUE)
    if (!is.na(mean_auc) && mean_auc > best_auc) {
      best_auc <- mean_auc
      best_params <- list(max_depth=md, eta=et, nrounds=nr)
    }
  }
  best_params
}

# ── Tune SVM ──────────────────────────────────────────────────────────────────
tune_svm <- function(X_tr, y_tr, cost_grid, gamma_grid,
                     k = 10, seed = 42) {
  folds    <- strat_kfold(y_tr, k = k, seed = seed)
  best_auc <- -Inf
  best_cost <- 1; best_gamma <- 0.01
  y_fac_all <- factor(ifelse(y_tr==1,"T2D","Control"),
                      levels=c("Control","T2D"))
  
  for (c_val in cost_grid) {
    for (g_val in gamma_grid) {
      aucs <- sapply(seq_along(folds), function(fi) {
        vi <- folds[[fi]]; ti <- setdiff(seq_len(nrow(X_tr)), vi)
        if (length(unique(y_tr[vi])) < 2) return(NA)
        yti_fac <- y_fac_all[ti]; yvi_fac <- y_fac_all[vi]
        cw <- table(yti_fac)
        fit <- tryCatch(
          e1071::svm(x=X_tr[ti,,drop=F], y=yti_fac,
                     kernel="radial", cost=c_val, gamma=g_val,
                     probability=TRUE,
                     class.weights=c(Control=1, T2D=cw["Control"]/cw["T2D"])),
          error=function(e) NULL)
        if (is.null(fit)) return(NA)
        p <- attr(predict(fit, X_tr[vi,,drop=F],
                          probability=TRUE), "probabilities")[,"T2D"]
        tryCatch(as.numeric(pROC::auc(
          pROC::roc(y_tr[vi], p, levels=c(0,1), direction="<", quiet=TRUE))),
          error=function(e) NA)
      })
      mean_auc <- mean(aucs, na.rm=TRUE)
      if (!is.na(mean_auc) && mean_auc > best_auc) {
        best_auc <- mean_auc; best_cost <- c_val; best_gamma <- g_val
      }
    }
  }
  list(cost=best_cost, gamma=best_gamma)
}

###############################################################################
# SECTION 3: OUTER LOO LOOP
###############################################################################

cat("--- Section 3: Outer LOO loop (116 iterations x 4 models) ---\n")
cat("This will take 30-90 minutes. Progress printed every 10 folds.\n\n")

# Storage for OOF probabilities
oof_en  <- numeric(n_total)
oof_rf  <- numeric(n_total)
oof_xgb <- numeric(n_total)
oof_svm <- numeric(n_total)

# Storage for best hyperparams per fold
hp_en  <- vector("list", n_total)
hp_rf  <- vector("list", n_total)
hp_xgb <- vector("list", n_total)
hp_svm <- vector("list", n_total)

# Set up parallel backend
if (USE_PARALLEL && N_CORES > 1) {
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  cat(sprintf("Parallel: %d cores registered\n\n", N_CORES))
} else {
  cat("Running sequentially (set USE_PARALLEL=TRUE to speed up)\n\n")
}

# Run LOO loop 
loo_results <- foreach(
  i = seq_len(n_total),
  .packages = c("glmnet","ranger","xgboost","e1071","pROC","themis","tidyverse"),
  .export   = c("filter_features","apply_smote","make_weights","strat_kfold",
                "tune_en","tune_rf","tune_xgb","tune_svm",
                "ALPHA_GRID","SMOTE_K","SMOTE_RATIO","K_INNER",
                "RF_MTRY_GRID","RF_NTREE","XGB_GRID",
                "SVM_COST_GRID","SVM_GAMMA_GRID")
) %dopar% {
  
  set.seed(2025 + i)
  
  test_idx  <- i
  train_idx <- setdiff(seq_len(n_total), i)
  
  X_train <- X[train_idx, , drop = FALSE]
  y_train <- y_int[train_idx]
  X_test  <- X[test_idx,  , drop = FALSE]
  
  # ── Step 1: Feature filtering ──────────────────────────────────────────────
  filt    <- filter_features(X_train, y_train)
  feat_cols <- filt$cols
  X_tr_f  <- filt$X
  X_te_f  <- X_test[, feat_cols, drop = FALSE]
  
  # ── Step 2: Scaling (fit on train only) ────────────────────────────────────
  tr_mean <- colMeans(X_tr_f)
  tr_sd   <- apply(X_tr_f, 2, sd); tr_sd[tr_sd == 0] <- 1
  X_tr_sc <- scale(X_tr_f, tr_mean, tr_sd)
  X_te_sc <- scale(X_te_f, tr_mean, tr_sd)
  
  # ── Step 3: SMOTE on training fold ─────────────────────────────────────────
  sm      <- apply_smote(X_tr_sc, y_train, k=SMOTE_K, over_ratio=SMOTE_RATIO)
  X_sm    <- sm$X; y_sm <- sm$y
  w_sm    <- make_weights(y_sm)
  
  # ── Step 4: Tune + predict EN ──────────────────────────────────────────────
  hp_en_i <- tune_en(X_tr_sc, y_train, ALPHA_GRID, k=K_INNER, seed=i*7)
  fit_en  <- tryCatch(suppressWarnings(
    glmnet(X_sm, y_sm, family="binomial",
           alpha=hp_en_i$alpha, lambda=hp_en_i$lambda,
           weights=w_sm, standardize=FALSE)),
    error=function(e) NULL)
  p_en <- as.numeric(
    if (!is.null(fit_en))
      predict(fit_en, X_te_sc, type="response", s=hp_en_i$lambda)
    else 0.5)[1]
  
  # ── Step 5: Tune + predict RF ──────────────────────────────────────────────
  hp_rf_i <- tune_rf(X_tr_sc, y_train, RF_MTRY_GRID,
                     ntree=RF_NTREE, k=K_INNER, seed=i*7)
  y_tr_fac <- factor(ifelse(y_train==1,"T2D","Control"),
                     levels=c("Control","T2D"))
  cw_rf    <- table(y_tr_fac)
  fit_rf   <- tryCatch(
    ranger::ranger(x=X_tr_sc, y=y_tr_fac,
                   num.trees=hp_rf_i$ntree,
                   mtry=min(hp_rf_i$mtry, ncol(X_tr_sc)),
                   probability=TRUE,
                   class.weights=c(Control=1,
                                   T2D=cw_rf["Control"]/cw_rf["T2D"]),
                   seed=i, verbose=FALSE),
    error=function(e) NULL)
  p_rf <- as.numeric(
    if (!is.null(fit_rf))
      predict(fit_rf, as.data.frame(X_te_sc))$predictions[,"T2D"]
    else 0.5)[1]
  
  # ── Step 6: Tune + predict XGBoost ────────────────────────────────────────
  hp_xgb_i <- tune_xgb(X_tr_sc, y_train, XGB_GRID,
                       k=K_INNER, seed=i*7)
  scale_pos <- sum(y_train==0) / sum(y_train==1)
  dtrain_full <- xgboost::xgb.DMatrix(X_sm, label=y_sm)
  fit_xgb  <- tryCatch(suppressWarnings(
    xgboost::xgb.train(
      params = list(objective="binary:logistic", eval_metric="auc",
                    max_depth=hp_xgb_i$max_depth, eta=hp_xgb_i$eta,
                    scale_pos_weight=scale_pos,
                    subsample=0.8, colsample_bytree=0.8),
      data=dtrain_full, nrounds=hp_xgb_i$nrounds,
      verbose=0, nthread=1)),
    error=function(e) NULL)
  p_xgb <- as.numeric(
    if (!is.null(fit_xgb))
      predict(fit_xgb, xgboost::xgb.DMatrix(X_te_sc))
    else 0.5)[1]
  
  # ── Step 7: Tune + predict SVM ─────────────────────────────────────────────
  hp_svm_i <- tune_svm(X_tr_sc, y_train, SVM_COST_GRID, SVM_GAMMA_GRID,
                       k=K_INNER, seed=i*7)
  cw_svm   <- table(y_tr_fac)
  fit_svm  <- tryCatch(
    e1071::svm(x=X_tr_sc,
               y=factor(ifelse(y_train==1,"T2D","Control"),
                        levels=c("Control","T2D")),
               kernel="radial",
               cost=hp_svm_i$cost, gamma=hp_svm_i$gamma,
               probability=TRUE,
               class.weights=c(Control=1,
                               T2D=cw_svm["Control"]/cw_svm["T2D"])),
    error=function(e) NULL)
  p_svm <- as.numeric(
    if (!is.null(fit_svm))
      attr(predict(fit_svm, X_te_sc, probability=TRUE),
           "probabilities")[,"T2D"]
    else 0.5)[1]
  
  list(p_en=p_en, p_rf=p_rf, p_xgb=p_xgb, p_svm=p_svm,
       hp_en=hp_en_i, hp_rf=hp_rf_i, hp_xgb=hp_xgb_i, hp_svm=hp_svm_i,
       i=i)
}

# Stop cluster
if (USE_PARALLEL && N_CORES > 1) stopCluster(cl)

# Unpack results
for (res in loo_results) {
  i <- res$i
  oof_en[i]  <- as.numeric(res$p_en)[1]
  oof_rf[i]  <- as.numeric(res$p_rf)[1]
  oof_xgb[i] <- as.numeric(res$p_xgb)[1]
  oof_svm[i] <- as.numeric(res$p_svm)[1]
  hp_en[[i]]  <- res$hp_en
  hp_rf[[i]]  <- res$hp_rf
  hp_xgb[[i]] <- res$hp_xgb
  hp_svm[[i]] <- res$hp_svm
}

cat("LOO loop complete.\n\n")

###############################################################################
# SECTION 4: ENSEMBLE
###############################################################################

cat("--- Section 4: Computing ensemble ---\n")

# Simple average ensemble (equal weight)
oof_ensemble <- (oof_en + oof_rf + oof_xgb + oof_svm) / 4

# Weighted ensemble: weight by individual OOF AUC
compute_auc <- function(probs, y) {
  tryCatch(as.numeric(pROC::auc(
    pROC::roc(y, probs, levels=c(0,1), direction="<", quiet=TRUE))),
    error=function(e) 0.5)
}

auc_en_raw  <- compute_auc(oof_en,  y_int)
auc_rf_raw  <- compute_auc(oof_rf,  y_int)
auc_xgb_raw <- compute_auc(oof_xgb, y_int)
auc_svm_raw <- compute_auc(oof_svm, y_int)

# Weighted ensemble (weights proportional to individual AUC - 0.5)
raw_weights <- pmax(c(auc_en_raw, auc_rf_raw,
                      auc_xgb_raw, auc_svm_raw) - 0.5, 0)
if (sum(raw_weights) == 0) raw_weights <- rep(1, 4)
norm_weights <- raw_weights / sum(raw_weights)

oof_ensemble_w <- (oof_en  * norm_weights[1] +
                     oof_rf  * norm_weights[2] +
                     oof_xgb * norm_weights[3] +
                     oof_svm * norm_weights[4])

cat(sprintf("Ensemble weights: EN=%.3f RF=%.3f XGB=%.3f SVM=%.3f\n\n",
            norm_weights[1], norm_weights[2],
            norm_weights[3], norm_weights[4]))

###############################################################################
# SECTION 5: PERFORMANCE METRICS — ALL MODELS
###############################################################################

cat("--- Section 5: Performance metrics ---\n\n")

model_names <- c("Elastic Net", "Random Forest",
                 "XGBoost", "SVM RBF",
                 "Ensemble (equal)", "Ensemble (weighted)")

oof_list <- list(oof_en, oof_rf, oof_xgb, oof_svm,
                 oof_ensemble, oof_ensemble_w)

results_list <- lapply(seq_along(oof_list), function(mi) {
  probs <- oof_list[[mi]]
  name  <- model_names[mi]
  
  roc_obj <- pROC::roc(y_int, probs,
                       levels=c(0,1), direction="<", quiet=TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val  <- as.numeric(pROC::ci.auc(roc_obj, conf.level=0.95))
  
  pr_obj  <- PRROC::pr.curve(
    scores.class0 = probs[y_int==1],
    scores.class1 = probs[y_int==0],
    curve=TRUE)
  pr_auc  <- pr_obj$auc.integral
  
  coords  <- pROC::coords(roc_obj, "best",
                          ret=c("threshold","sensitivity","specificity"))
  thresh  <- if(is.finite(coords$threshold[1])) as.numeric(coords$threshold[1]) else 0.5
  pred_c  <- ifelse(probs >= thresh, 1L, 0L)
  cm      <- table(Predicted=pred_c, Observed=y_int)
  
  TP <- if("1" %in% rownames(cm) && "1" %in% colnames(cm)) cm["1","1"] else 0L
  TN <- if("0" %in% rownames(cm) && "0" %in% colnames(cm)) cm["0","0"] else 0L
  FP <- if("1" %in% rownames(cm) && "0" %in% colnames(cm)) cm["1","0"] else 0L
  FN <- if("0" %in% rownames(cm) && "1" %in% colnames(cm)) cm["0","1"] else 0L
  sens <- TP / max(TP + FN, 1)
  spec <- TN / max(TN + FP, 1)
  
  cat(sprintf("%-22s AUC=%.3f [%.3f-%.3f] PR-AUC=%.3f Sens=%.3f Spec=%.3f\n",
              name, auc_val, ci_val[1], ci_val[3],
              pr_auc, sens, spec))
  
  list(name=name, auc=auc_val, ci_lo=ci_val[1], ci_hi=ci_val[3],
       pr_auc=pr_auc, sens=sens, spec=spec, thresh=thresh,
       roc_obj=roc_obj, pr_obj=pr_obj, probs=probs)
})

names(results_list) <- model_names

# Find best model
best_idx  <- which.max(sapply(results_list, `[[`, "auc"))
best_name <- model_names[best_idx]
best_probs <- oof_list[[best_idx]]

cat(sprintf("\nBest model: %s (AUC = %.3f)\n\n", best_name,
            results_list[[best_idx]]$auc))

baseline_pr <- mean(y_int)

###############################################################################
# SECTION 6: PERMUTATION TEST ON BEST MODEL
###############################################################################

cat(sprintf("--- Section 6: Permutation test on %s (%d shuffles) ---\n",
            best_name, N_PERMS))

perm_aucs <- numeric(N_PERMS)

for (perm in seq_len(N_PERMS)) {
  set.seed(perm * 17)
  y_perm <- sample(y_int)
  perm_aucs[perm] <- tryCatch(
    as.numeric(pROC::auc(pROC::roc(
      y_perm, best_probs,
      levels=c(0,1), direction="<", quiet=TRUE))),
    error=function(e) 0.5)
  if (perm %% 25 == 0)
    cat(sprintf("  Permutation %d/%d done\n", perm, N_PERMS))
}

perm_p <- mean(perm_aucs >= results_list[[best_idx]]$auc)

cat(sprintf("\nPermutation p-value: %.4f\n", perm_p))
cat(sprintf("Null AUC: %.3f ± %.3f\n\n",
            mean(perm_aucs), sd(perm_aucs)))

###############################################################################
# SECTION 7: HYPERPARAMETER SUMMARY
###############################################################################

cat("--- Section 7: Hyperparameter summary ---\n\n")

en_alphas   <- sapply(hp_en,  `[[`, "alpha")
en_lambdas  <- sapply(hp_en,  `[[`, "lambda")
rf_mtrys    <- sapply(hp_rf,  `[[`, "mtry")
xgb_depths  <- sapply(hp_xgb, `[[`, "max_depth")
xgb_etas    <- sapply(hp_xgb, `[[`, "eta")
xgb_nrounds <- sapply(hp_xgb, `[[`, "nrounds")
svm_costs   <- sapply(hp_svm, `[[`, "cost")
svm_gammas  <- sapply(hp_svm, `[[`, "gamma")

cat("EN alpha distribution:\n");  print(table(en_alphas))
cat(sprintf("EN lambda: median=%.5f range=%.5f-%.5f\n\n",
            median(en_lambdas), min(en_lambdas), max(en_lambdas)))
cat("RF mtry distribution:\n");   print(table(rf_mtrys))
cat("\nXGB max_depth distribution:\n"); print(table(xgb_depths))
cat("XGB eta distribution:\n");   print(table(xgb_etas))
cat("XGB nrounds distribution:\n"); print(table(xgb_nrounds))
cat("\nSVM cost distribution:\n"); print(table(svm_costs))
cat("SVM gamma distribution:\n"); print(table(svm_gammas))

###############################################################################
# SECTION 8: TRAIN FINAL FROZEN MODEL (best model on all 116 samples)
###############################################################################

cat("\n--- Section 8: Training final frozen model ---\n")

# Use median hyperparameters from LOO folds
best_alpha_final  <- as.numeric(names(which.max(table(en_alphas))))
best_lambda_final <- median(en_lambdas)
best_mtry_final   <- as.numeric(names(which.max(table(rf_mtrys))))
best_depth_final  <- as.numeric(names(which.max(table(xgb_depths))))
best_eta_final    <- as.numeric(names(which.max(table(xgb_etas))))
best_nrounds_final <- as.numeric(names(which.max(table(xgb_nrounds))))
best_cost_final   <- as.numeric(names(which.max(table(svm_costs))))
best_gamma_final  <- as.numeric(names(which.max(table(svm_gammas))))

# Full data preprocessing
filt_full    <- filter_features(X, y_int)
X_full_f     <- filt_full$X
full_means   <- colMeans(X_full_f)
full_sds     <- apply(X_full_f, 2, sd); full_sds[full_sds==0] <- 1
X_full_sc    <- scale(X_full_f, full_means, full_sds)
sm_full      <- apply_smote(X_full_sc, y_int, k=SMOTE_K,
                            over_ratio=SMOTE_RATIO)
w_full       <- make_weights(sm_full$y)
y_full_fac   <- factor(ifelse(y_int==1,"T2D","Control"),
                       levels=c("Control","T2D"))
cw_full      <- table(y_full_fac)

# Train all 4 final models
final_en <- glmnet(sm_full$X, sm_full$y, family="binomial",
                   alpha=best_alpha_final, lambda=best_lambda_final,
                   weights=w_full, standardize=FALSE)

final_rf <- ranger::ranger(x=X_full_sc, y=y_full_fac,
                           num.trees=RF_NTREE,
                           mtry=min(best_mtry_final, ncol(X_full_sc)),
                           probability=TRUE,
                           class.weights=c(Control=1,
                                           T2D=cw_full["Control"]/cw_full["T2D"]),
                           importance="impurity", verbose=FALSE)

dtrain_final <- xgboost::xgb.DMatrix(sm_full$X, label=sm_full$y)
final_xgb <- xgboost::xgb.train(
  params=list(objective="binary:logistic", eval_metric="auc",
              max_depth=best_depth_final, eta=best_eta_final,
              scale_pos_weight=sum(y_int==0)/sum(y_int==1),
              subsample=0.8, colsample_bytree=0.8),
  data=dtrain_final, nrounds=best_nrounds_final, verbose=0)

final_svm <- e1071::svm(x=X_full_sc, y=y_full_fac,
                        kernel="radial",
                        cost=best_cost_final, gamma=best_gamma_final,
                        probability=TRUE,
                        class.weights=c(Control=1,
                                        "T2D"=as.numeric(cw_full["Control"])/as.numeric(cw_full["T2D"])))
cat("All 4 final models trained on full 116 samples.\n\n")
###############################################################################
# SECTION 9: SAVE ALL OUTPUTS
###############################################################################

cat("--- Section 9: Saving outputs ---\n")

# OOF predictions
oof_df <- data.frame(
  Sample_Index  = seq_len(n_total),
  True_Label    = as.character(y),
  True_Int      = y_int,
  P_T2D_EN      = round(oof_en,  4),
  P_T2D_RF      = round(oof_rf,  4),
  P_T2D_XGB     = round(oof_xgb, 4),
  P_T2D_SVM     = round(oof_svm, 4),
  P_T2D_Ensemble_Equal    = round(oof_ensemble,   4),
  P_T2D_Ensemble_Weighted = round(oof_ensemble_w, 4)
)
fwrite(oof_df, file.path(OUT_DIR, "oof_predictions.csv"))

# Performance summary
perf_df <- do.call(rbind, lapply(results_list, function(r) {
  data.frame(Model=r$name, AUC=round(r$auc,4),
             CI_lower=round(r$ci_lo,4), CI_upper=round(r$ci_hi,4),
             PR_AUC=round(r$pr_auc,4),
             Sensitivity=round(r$sens,4), Specificity=round(r$spec,4))
}))
fwrite(perf_df, file.path(OUT_DIR, "performance_summary.csv"))
cat("Performance summary:\n"); print(perf_df); cat("\n")

# Permutation results
perm_df <- data.frame(Best_Model=best_name,
                      Observed_AUC=results_list[[best_idx]]$auc,
                      Perm_p_value=perm_p,
                      Null_AUC_mean=mean(perm_aucs),
                      Null_AUC_sd=sd(perm_aucs))
fwrite(perm_df, file.path(OUT_DIR, "permutation_test.csv"))

# Hyperparameter logs
hp_df <- data.frame(
  Fold=seq_len(n_total),
  EN_alpha=en_alphas, EN_lambda=round(en_lambdas,6),
  RF_mtry=rf_mtrys,
  XGB_depth=xgb_depths, XGB_eta=xgb_etas, XGB_nrounds=xgb_nrounds,
  SVM_cost=svm_costs, SVM_gamma=svm_gammas
)
fwrite(hp_df, file.path(OUT_DIR, "hyperparameters_per_fold.csv"))

# Scaling parameters (needed for external validation)
scaling_params <- list(
  feature_names = filt_full$cols,
  means         = full_means,
  sds           = full_sds,
  smote_k       = SMOTE_K,
  smote_ratio   = SMOTE_RATIO
)
saveRDS(scaling_params, file.path(OUT_DIR, "final_scaling_params.rds"))

# Final models
saveRDS(final_en,  file.path(OUT_DIR, "final_model_EN.rds"))
saveRDS(final_rf,  file.path(OUT_DIR, "final_model_RF.rds"))
saveRDS(final_xgb, file.path(OUT_DIR, "final_model_XGB.rds"))
saveRDS(final_svm, file.path(OUT_DIR, "final_model_SVM.rds"))
saveRDS(scaling_params, file.path(OUT_DIR, "final_scaling_params.rds"))

cat(sprintf("All outputs saved to: %s\n\n", OUT_DIR))

###############################################################################
# SECTION 10: DIAGNOSTIC PLOTS (PDF)
###############################################################################

cat("--- Section 10: Generating diagnostic plots ---\n")

model_cols <- c("#E74C3C","#27AE60","#F39C12","#8E44AD","#2980B9","#1A5276")

pdf(file.path(OUT_DIR, "phase5_diagnostic_plots.pdf"),
    width=16, height=12)
par(mfrow=c(2,3), mar=c(5,4,4,2)+0.1)

# Plot 1: ROC curves all models
plot(results_list[[1]]$roc_obj, col=model_cols[1], lwd=2,
     main="ROC Curves — All Models (LOO-CV)",
     xlab="1 - Specificity", ylab="Sensitivity",
     legacy.axes=TRUE)
for (mi in 2:length(results_list)) {
  plot(results_list[[mi]]$roc_obj, col=model_cols[mi],
       lwd=2, add=TRUE)
}
abline(a=0, b=1, lty=2, col="grey60")
legend("bottomright",
       legend=sprintf("%s (%.3f)", model_names,
                      sapply(results_list, `[[`, "auc")),
       col=model_cols, lwd=2, bty="n", cex=0.75)

# Plot 2: AUC comparison bar chart
aucs_all <- sapply(results_list, `[[`, "auc")
ci_lo    <- sapply(results_list, `[[`, "ci_lo")
ci_hi    <- sapply(results_list, `[[`, "ci_hi")
bp <- barplot(aucs_all, col=model_cols, border=NA,
              ylim=c(0.4, 1.0),
              main="AUC Comparison with 95% CI",
              ylab="ROC-AUC", las=2,
              names.arg=c("EN","RF","XGB","SVM","Ens-E","Ens-W"),
              cex.names=0.85)
arrows(bp, ci_lo, bp, ci_hi, angle=90, code=3, length=0.05, lwd=1.5)
abline(h=0.5, lty=2, col="grey60")
text(bp, aucs_all+0.01, sprintf("%.3f", aucs_all),
     cex=0.75, font=2)

# Plot 3: PR curves
plot(results_list[[1]]$pr_obj, col=model_cols[1], lwd=2,
     main=sprintf("PR Curves (baseline=%.3f)", baseline_pr),
     xlab="Recall", ylab="Precision", auc.main=FALSE)
for (mi in 2:length(results_list)) {
  plot(results_list[[mi]]$pr_obj, col=model_cols[mi],
       lwd=2, add=TRUE, auc.main=FALSE)
}
abline(h=baseline_pr, lty=2, col="grey60")
legend("topright",
       legend=sprintf("%s (%.3f)",
                      c("EN","RF","XGB","SVM","Ens-E","Ens-W"),
                      sapply(results_list, `[[`, "pr_auc")),
       col=model_cols, lwd=2, bty="n", cex=0.75)

# Plot 4: Permutation test
hist(perm_aucs, breaks=20, col="#BDC3C7", border="white",
     main=sprintf("Permutation Test: %s\np = %.4f",
                  best_name, perm_p),
     xlab="AUC under null", ylab="Frequency")
abline(v=results_list[[best_idx]]$auc,
       col=model_cols[best_idx], lwd=2.5)
abline(v=mean(perm_aucs), col="#7F8C8D", lwd=1.5, lty=2)
legend("topright",
       legend=c(sprintf("Observed=%.3f",
                        results_list[[best_idx]]$auc),
                sprintf("Null=%.3f", mean(perm_aucs))),
       col=c(model_cols[best_idx],"#7F8C8D"),
       lwd=c(2.5,1.5), bty="n", cex=0.85)

# Plot 5: Predicted probability distributions (best model)
bp_probs <- split(best_probs, ifelse(y_int==1,"T2D","Control"))
boxplot(bp_probs,
        col=c("#AED6F1","#F1948A"),
        border=c("#2980B9","#E74C3C"),
        main=sprintf("Predicted P(T2D) by Class\n%s", best_name),
        ylab="Predicted P(T2D)", ylim=c(0,1), lwd=1.5)
stripchart(bp_probs, vertical=TRUE, method="jitter",
           jitter=0.08, pch=16, cex=0.6,
           col=c("#2980B950","#E74C3C50"), add=TRUE)

# Plot 6: RF feature importance (if RF is available)
if (!is.null(final_rf$variable.importance)) {
  imp <- sort(final_rf$variable.importance, decreasing=TRUE)
  imp_top <- head(imp, min(20, length(imp)))
  barplot(rev(imp_top), horiz=TRUE,
          names.arg=rev(names(imp_top)),
          las=1, cex.names=0.75, col="#85C1E9", border=NA,
          main="RF Feature Importance\n(Top 20 MEs)",
          xlab="Impurity decrease")
} else {
  plot.new()
  text(0.5, 0.5, "Feature importance\nnot available",
       cex=1.2, col="grey50")
}

dev.off()
cat(sprintf("Saved: %s/phase5_diagnostic_plots.pdf\n\n", OUT_DIR))

###############################################################################
# FINAL SUMMARY
###############################################################################

cat("=======================================================================\n")
cat("   FINAL RESULTS\n")
cat("=======================================================================\n")
cat(sprintf("  Dataset:      Blood MEs, %d samples (%d T2D, %d Control)\n",
            n_total, n_t2d, n_ctrl))
cat(sprintf("  Outer CV:     LOO (%d iterations)\n", n_total))
cat(sprintf("  Inner CV:     Stratified %d-fold\n", K_INNER))
cat(sprintf("  SMOTE:        K=%d, ratio=%.1f\n", SMOTE_K, SMOTE_RATIO))
cat(sprintf("  Models:       EN, RF, XGB, SVM + 2 ensembles\n\n"))
for (mi in seq_along(results_list)) {
  r <- results_list[[mi]]
  cat(sprintf("  %-24s AUC=%.3f [%.3f-%.3f] PR=%.3f\n",
              r$name, r$auc, r$ci_lo, r$ci_hi, r$pr_auc))
}
cat(sprintf("\n  Best model:   %s\n", best_name))
cat(sprintf("  Perm p-value: %.4f\n\n", perm_p))
cat(sprintf("  Outputs in:   %s\n", OUT_DIR))
cat("=======================================================================\n")
cat("===  STEP 5 COMPLETED SUCCESSFULLY ===\n")