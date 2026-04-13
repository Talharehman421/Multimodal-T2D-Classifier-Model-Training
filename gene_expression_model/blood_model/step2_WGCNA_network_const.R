###############################################################################
# BLOOD MODEL BUILDING — STEP 2
# WGCNA Network Construction + Module Detection
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

ROOT    <- "C:/Users/talha/OneDrive/Desktop/gene expression work"
BM_DIR  <- file.path(ROOT, "blood model building")
WGCNA_DIR <- file.path(BM_DIR, "wgcna")
dir.create(WGCNA_DIR,                          showWarnings=FALSE)
dir.create(file.path(WGCNA_DIR, "plots"),      showWarnings=FALSE)

META_FILE <- file.path(BM_DIR, "meta", "blood_samples_master.csv")
EXPR_FILE <- file.path(BM_DIR, "processed", "log_blood_expr_top5000.csv")

stopifnot(file.exists(EXPR_FILE))
stopifnot(file.exists(META_FILE))

cat("=======================================================================\n")
cat("  BLOOD MODEL — STEP 2: WGCNA NETWORK CONSTRUCTION\n")
cat("=======================================================================\n\n")

###############################################################################
# SECTION 1: LOAD DATA
###############################################################################

cat("--- Section 1: Loading data ---\n")

expr    <- fread(EXPR_FILE)
meta    <- fread(META_FILE)

# WGCNA needs samples x genes
datExpr <- t(as.matrix(expr[, -1]))
colnames(datExpr) <- expr$Gene
mode(datExpr) <- "numeric"

cat(sprintf("datExpr: %d samples x %d genes\n\n",
            nrow(datExpr), ncol(datExpr)))

###############################################################################
# SECTION 2: WGCNA QUALITY CHECK
###############################################################################

cat("--- Section 2: WGCNA goodSamplesGenes check ---\n")

gsg <- goodSamplesGenes(datExpr, verbose=3)

if (!gsg$allOK) {
  cat("Removing flagged samples/genes...\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("After QC: %d samples x %d genes\n", 
              nrow(datExpr), ncol(datExpr)))
} else {
  cat("All samples and genes passed QC.\n")
}
cat("\n")

###############################################################################
# SECTION 3: SOFT THRESHOLD SELECTION
###############################################################################

cat("--- Section 3: Soft threshold power selection ---\n")

powers <- c(1:10, seq(12, 20, by=2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector  = powers,
  networkType  = "signed hybrid",   # signed hybrid better for RNA-seq
  corFnc       = "bicor",           # bicor robust to outliers
  verbose      = 3
)

# Save results
sft_df <- data.table(
  Power            = sft$fitIndices[,1],
  ScaleFreeR2      = -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  MeanConnectivity = sft$fitIndices[,5]
)
fwrite(sft_df, file.path(WGCNA_DIR, "soft_threshold_results.csv"))


good_powers <- sft_df$Power[sft_df$ScaleFreeR2 >= 0.85]
if (length(good_powers) == 0) {

  good_powers <- sft_df$Power[sft_df$ScaleFreeR2 >= 0.80]
}
if (length(good_powers) == 0) {
  
  softPower <- sft$powerEstimate
  cat("WARNING: R2 never reached 0.80. Using pickSoftThreshold recommendation.\n")
} else {
  softPower <- min(good_powers)
}

cat(sprintf("Selected softPower: %d\n\n", softPower))

# Plot: Soft threshold selection
p_sft <- ggplot(sft_df, aes(x=Power)) +
  geom_line(aes(y=ScaleFreeR2), color="#E74C3C", linewidth=1.2) +
  geom_point(aes(y=ScaleFreeR2), color="#E74C3C", size=2.5) +
  geom_hline(yintercept=0.85, linetype="dashed", color="#E74C3C", alpha=0.5) +
  geom_hline(yintercept=0.80, linetype="dotted", color="grey50") +
  geom_vline(xintercept=softPower, linetype="dashed",
             color="#2980B9", linewidth=1) +
  annotate("text", x=softPower+0.3, y=0.2,
           label=sprintf("Selected\npower=%d", softPower),
           color="#2980B9", size=3.5, hjust=0) +
  scale_x_continuous(breaks=powers) +
  ylim(0, 1) +
  theme_minimal(base_size=12) +
  labs(title="Soft Threshold Power Selection",
       subtitle="Red dashed = R²=0.85 target | Blue = selected power",
       x="Soft Threshold Power", y="Scale-Free Topology R²")

p_mean <- ggplot(sft_df, aes(x=Power, y=MeanConnectivity)) +
  geom_line(color="#27AE60", linewidth=1.2) +
  geom_point(color="#27AE60", size=2.5) +
  geom_vline(xintercept=softPower, linetype="dashed",
             color="#2980B9", linewidth=1) +
  scale_x_continuous(breaks=powers) +
  theme_minimal(base_size=12) +
  labs(title="Mean Connectivity vs Power",
       x="Soft Threshold Power", y="Mean Connectivity")

library(gridExtra)
p_combined <- arrangeGrob(p_sft, p_mean, ncol=2)
ggsave(file.path(WGCNA_DIR, "plots", "step2_soft_threshold.png"),
       p_combined, width=12, height=5, dpi=150)

###############################################################################
# SECTION 4: NETWORK CONSTRUCTION
###############################################################################

cat("--- Section 4: Network construction (blockwiseModules) ---\n")
cat("This is the slowest step — typically 10-30 minutes.\n\n")

net <- blockwiseModules(
  datExpr,
  power           = softPower,
  networkType     = "signed hybrid",   
  corType         = "bicor",
  TOMType         = "signed",
  minModuleSize   = 20,        
  mergeCutHeight  = 0.30,      
  numericLabels   = FALSE,     
  pamRespectsDendro = TRUE,
  saveTOMs        = FALSE,
  verbose         = 5,
  nThreads        = 0          
)

# Module summary
module_colors <- net$colors
module_table  <- table(module_colors)
n_modules     <- sum(names(module_table) != "grey")
pct_grey      <- round(100 * sum(module_colors=="grey") / length(module_colors), 1)

cat(sprintf("\nModules found (excluding grey): %d\n", n_modules))
cat(sprintf("Genes in grey (unassigned): %d (%.1f%%)\n\n",
            sum(module_colors=="grey"), pct_grey))
print(sort(module_table, decreasing=TRUE))

if (pct_grey > 60) {
  cat("\nWARNING: >60% genes in grey module.\n")
  cat("Consider: lowering minModuleSize further, or checking softPower.\n\n")
}

###############################################################################
# SECTION 5: MODULE EIGENGENES
###############################################################################

cat("\n--- Section 5: Computing module eigengenes ---\n")

MEs     <- moduleEigengenes(datExpr, colors=module_colors)$eigengenes
MEs     <- orderMEs(MEs)

# Check ME variances — should NOT all be equal
me_vars <- apply(MEs, 2, var)
cat("ME variances:\n")
print(round(me_vars, 4))
cat(sprintf("\nMin ME variance: %.6f\n", min(me_vars)))

if (max(me_vars) < 0.01) {
  cat("WARNING: All ME variances near zero — check expression data quality\n")
} else {
  cat("ME variances look healthy.\n\n")
}

###############################################################################
# SECTION 6: PLOTS
###############################################################################

cat("--- Section 6: Generating plots ---\n")

# Plot: Module sizes
module_df <- as.data.frame(module_table)
colnames(module_df) <- c("Module","GeneCount")
module_df$IsGrey <- module_df$Module == "grey"
module_df <- module_df[order(-module_df$GeneCount),]

p_mod <- ggplot(module_df,
                aes(x=reorder(Module, -GeneCount),
                    y=GeneCount, fill=IsGrey)) +
  geom_col(show.legend=FALSE) +
  scale_fill_manual(values=c("FALSE"="#2980B9","TRUE"="#BDC3C7")) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title=sprintf("Module Sizes (n=%d modules + grey)",n_modules),
       subtitle=sprintf("softPower=%d, minModuleSize=20, mergeCutHeight=0.30",
                        softPower),
       x="Module", y="Number of Genes")
ggsave(file.path(WGCNA_DIR, "plots", "step2_module_sizes.png"),
       p_mod, width=10, height=5, dpi=150)

# Plot: ME correlation heatmap
png(file.path(WGCNA_DIR, "plots", "step2_ME_heatmap.png"),
    width=800, height=700)
plotEigengeneNetworks(MEs, "Module Eigengene Network",
                      marHeatmap=c(3,4,2,2),
                      plotDendrograms=TRUE,
                      xLabelsAngle=90)
dev.off()

# Plot: ME variance bar chart
me_var_df <- data.frame(
  ME       = names(me_vars),
  Variance = as.numeric(me_vars)
)
me_var_df <- me_var_df[order(-me_var_df$Variance),]

p_var <- ggplot(me_var_df, aes(x=reorder(ME, -Variance), y=Variance)) +
  geom_col(fill="#27AE60") +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Module Eigengene Variances",
       subtitle="All should be > 0 for a healthy network",
       x="Module Eigengene", y="Variance")
ggsave(file.path(WGCNA_DIR, "plots", "step2_ME_variances.png"),
       p_var, width=10, height=5, dpi=150)

###############################################################################
# SECTION 7: SAVE OUTPUTS
###############################################################################

cat("--- Section 7: Saving outputs ---\n")

# Gene-module assignments
gene_modules <- data.table(
  Gene   = colnames(datExpr),
  Module = module_colors
)
fwrite(gene_modules,
       file.path(WGCNA_DIR, "gene_module_assignments.csv"))

# Module eigengenes (samples x MEs)
MEs_out <- data.table(Sample_ID=rownames(MEs), MEs)
fwrite(MEs_out, file.path(WGCNA_DIR, "module_eigengenes.csv"))

# Full WGCNA object
save(net, datExpr, softPower, MEs, module_colors,
     file=file.path(WGCNA_DIR, "WGCNA_blood_network.RData"))

cat(sprintf("Saved: WGCNA_blood_network.RData\n"))
cat(sprintf("Saved: gene_module_assignments.csv (%d genes)\n",
            nrow(gene_modules)))
cat(sprintf("Saved: module_eigengenes.csv (%d samples x %d MEs)\n\n",
            nrow(MEs_out), ncol(MEs)-1))

###############################################################################
# FINAL SUMMARY
###############################################################################

cat("=======================================================================\n")
cat("  STEP 2 SUMMARY\n")
cat("=======================================================================\n")
cat(sprintf("  Genes input:      5000\n"))
cat(sprintf("  Soft power:       %d\n", softPower))
cat(sprintf("  Network type:     signed hybrid\n"))
cat(sprintf("  minModuleSize:    20\n"))
cat(sprintf("  mergeCutHeight:   0.30\n"))
cat(sprintf("  Modules found:    %d (excl. grey)\n", n_modules))
cat(sprintf("  Genes in grey:    %.1f%%\n", pct_grey))
cat(sprintf("  ME variance OK:   %s\n",
            ifelse(min(me_vars)>0.001,"YES","CHECK CAREFULLY")))
cat("=======================================================================\n")
cat("=== STEP 2 COMPLETED ===\n")
