###############################################################################
#STEP 3
# WGCNA with automatic batch confounder diagnostics 
###############################################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(matrixStats)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

cat("=== PHASE 1 — STEP 3 : WGCNA STARTED ===\n")

###############################################################################
# 0) Load Step 2 outputs
###############################################################################
load("processed/step2_outputs_optionC.RData")


datExpr <- t(vst_wgcna)   

cat("Expression matrix:", dim(datExpr)[1], "samples x",
    dim(datExpr)[2], "genes\n")

###############################################################################
# 1) Sample QC & outlier diagnostics 
###############################################################################
cat("Running sample clustering QC...\n")

sampleTree <- hclust(dist(datExpr), method = "average")

dir.create("wgcna/plots", recursive = TRUE, showWarnings = FALSE)

png("wgcna/plots/sample_clustering.png", width=1200, height=800)
plot(sampleTree, main="Sample clustering (QC)", sub="", xlab="")
dev.off()

###############################################################################
# PCA-based outlier check (NO auto-removal)
###############################################################################

cat("Running PCA-based outlier diagnostics...\n")


dir.create("wgcna/logs", recursive = TRUE, showWarnings = FALSE)

pca <- prcomp(datExpr, scale. = TRUE)
z_scores <- scale(pca$x[,1])

outliers <- which(abs(z_scores) > 3)

write.table(
  data.frame(
    GSM_ID = rownames(datExpr)[outliers],
    Z_PC1  = z_scores[outliers]
  ),
  "wgcna/logs/potential_outliers_PC1.txt",
  row.names = FALSE,
  quote = FALSE
)

cat("Outlier diagnostics saved (no samples removed).\n")
cat("Number of potential outliers:", length(outliers), "\n\n")


###############################################################################
# 2) Trait preparation
###############################################################################
cat("Preparing trait matrix...\n")

traits <- meta %>%
  transmute(
    Disease = as.numeric(factor(Disease_Status)),
    Age = Age,
    BMI = BMI,
    Sex = ifelse(Sex == "male", 1, 0)
  )

rownames(traits) <- meta$GSM_ID

###############################################################################
# 3) Batch / technical confounder diagnostics 
###############################################################################
cat("Checking batch confounding...\n")

batch_cols <- intersect(c("Batch", "Run", "Lane", "Flowcell"), colnames(meta))

batch_confounded <- FALSE

if (length(batch_cols) > 0) {
  for (b in batch_cols) {
    if (length(unique(meta[[b]])) > 1) {
      tab <- table(meta[[b]], meta$Disease_Status)
      p <- suppressWarnings(chisq.test(tab)$p.value)
      if (!is.na(p) && p < 0.05) {
        batch_confounded <- TRUE
        cat("Batch variable", b, "is confounded with disease (p =", p, ")\n")
      }
    }
  }
}

if (batch_confounded) {
  cat("Decision: DO NOT remove batch (risk of signal loss).\n\n")
} else {
  cat("No strong batch-disease confounding detected.\n")
  cat("Decision: retain original expression (no correction applied).\n\n")
}

###############################################################################
# 4) Soft-threshold selection 
###############################################################################
cat("Selecting soft-threshold power...\n")

powers <- c(1:20)
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

png("wgcna/plots/soft_threshold.png", width=1200, height=600)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n",
     main="Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     col="red")
abline(h=0.80, col="blue")

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n",
     main="Mean connectivity")
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels=powers,
     col="red")
dev.off()

softPower <- min(sft$fitIndices[sft$fitIndices[,2] >= 0.80, 1])
if (is.infinite(softPower)) softPower <- 8

cat("Chosen softPower =", softPower, "\n\n")

###############################################################################
# 5) Network construction & module detection
###############################################################################
cat("Running blockwiseModules...\n")

net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  deepSplit = 2,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- net$colors
MEs <- net$MEs

###############################################################################
# 6) Module–trait relationships
###############################################################################
cat("Computing module–trait correlations...\n")

MEs_ordered <- orderMEs(MEs)
moduleTraitCor <- cor(MEs_ordered, traits, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

png("wgcna/plots/module_trait_heatmap.png", width=1200, height=800)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits),
  yLabels = names(MEs_ordered),
  ySymbols = names(MEs_ordered),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = signif(moduleTraitCor, 2),
  setStdMargins = FALSE,
  main = "Module–trait relationships"
)
dev.off()

###############################################################################
# 7) Save outputs (reproducibility)
###############################################################################
dir.create("wgcna/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("wgcna/logs", recursive = TRUE, showWarnings = FALSE)

write.table(
  data.frame(Gene = colnames(datExpr),
             Module = moduleColors),
  "wgcna/tables/gene_module_membership.txt",
  row.names = FALSE,
  quote = FALSE
)

save(
  net,
  moduleColors,
  MEs,
  moduleTraitCor,
  moduleTraitP,
  traits,
  file = "wgcna/step3_wgcna_results.RData"
)

cat("=== STEP 3 COMPLETED SUCCESSFULLY ===\n")
