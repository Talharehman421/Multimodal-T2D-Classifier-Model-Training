###############################################################################
# BLOOD MODEL BUILDING — STEP 3
# Hub Gene Identification
#
# INPUT:  WGCNA_blood_network.RData  (from Step 2)
# OUTPUT: blood_hub_genes.csv
#         plots: kME distributions, hub counts per module, connectivity
#
# KEY FIXES vs old pipeline:
#   - kME threshold 0.65 (was 0.60 — slightly stricter for cleaner hubs)
#   - Reports hub count per module for transparency
#   - Saves full kME table for downstream use
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
dir.create(HUB_DIR,                        showWarnings=FALSE)
dir.create(file.path(HUB_DIR, "plots"),    showWarnings=FALSE)

stopifnot(file.exists(file.path(WGCNA_DIR, "WGCNA_blood_network.RData")))

cat("=======================================================================\n")
cat("  BLOOD MODEL — STEP 3: HUB GENE IDENTIFICATION\n")
cat("=======================================================================\n\n")

###############################################################################
# SECTION 1: LOAD WGCNA OBJECT
###############################################################################

cat("--- Section 1: Loading WGCNA object ---\n")

load(file.path(WGCNA_DIR, "WGCNA_blood_network.RData"))
# loads: net, datExpr, softPower, MEs, module_colors

cat(sprintf("Loaded: %d samples x %d genes\n",
            nrow(datExpr), ncol(datExpr)))

modules_all  <- unique(module_colors)
modules_real <- modules_all[modules_all != "grey"]
cat(sprintf("Modules (excl. grey): %d\n\n", length(modules_real)))

###############################################################################
# SECTION 2: MODULE MEMBERSHIP (kME)
###############################################################################

cat("--- Section 2: Computing module membership (kME) ---\n")

kME_mat <- as.data.frame(
  cor(datExpr, MEs, use="pairwise.complete.obs", method="pearson")
)

# Rename columns: MEturquoise → kMEturquoise
colnames(kME_mat) <- paste0("kME", substring(colnames(kME_mat), 3))

kME_dt <- as.data.table(kME_mat, keep.rownames="Gene")
kME_dt$Module <- module_colors[kME_dt$Gene]

fwrite(kME_dt, file.path(HUB_DIR, "full_kME_table.csv"))
cat(sprintf("kME computed for %d genes across %d MEs\n\n",
            nrow(kME_dt), ncol(MEs)))

###############################################################################
# SECTION 3: INTRAMODULAR CONNECTIVITY (kWithin)
###############################################################################

cat("--- Section 3: Computing intramodular connectivity ---\n")

ADJ <- adjacency(
  datExpr,
  power   = softPower,
  type    = "signed",
  corFnc  = "bicor"
)

IMConn <- intramodularConnectivity(ADJ, colors=module_colors)
IMConn <- as.data.table(IMConn, keep.rownames="Gene")

cat("Intramodular connectivity computed.\n\n")

###############################################################################
# SECTION 4: MERGE + SELECT HUBS
###############################################################################

cat("--- Section 4: Selecting hub genes ---\n")

hub_df <- merge(kME_dt, IMConn, by="Gene")
hub_df <- hub_df[Module != "grey"]

KME_THRESH   <- 0.65   # module membership threshold
KWITHIN_PCTL <- 0.75   # top 25% connectivity within module

hub_list <- list()

for (mod in modules_real) {
  
  sub      <- hub_df[Module == mod]
  kME_col  <- paste0("kME", mod)
  
  if (!kME_col %in% names(sub)) {
    cat(sprintf("  SKIP %s — kME column not found\n", mod))
    next
  }
  
  kw_thresh <- quantile(sub$kWithin, KWITHIN_PCTL, na.rm=TRUE)
  
  hubs <- sub[abs(get(kME_col)) >= KME_THRESH & kWithin >= kw_thresh]
  
  cat(sprintf("  Module %-15s: %4d genes → %3d hubs (kME>=%.2f, kWithin>=%.3f)\n",
              mod, nrow(sub), nrow(hubs), KME_THRESH, kw_thresh))
  
  if (nrow(hubs) > 0) {
    # Add the kME value for this module as its own column
    hubs$kME_own <- hubs[[kME_col]]
    hub_list[[mod]] <- hubs[, .(Gene, Module, kME_own, kWithin,
                                kTotal, kDiff)]
  }
}

blood_hubs <- rbindlist(hub_list, use.names=TRUE, fill=TRUE)
setnames(blood_hubs, "kME_own", "kME")

cat(sprintf("\nTotal hub genes: %d across %d modules\n\n",
            nrow(blood_hubs), length(hub_list)))

if (nrow(blood_hubs) == 0) {
  stop("ERROR: No hub genes found. Lower KME_THRESH or KwithinPctl in Step 3.")
}

###############################################################################
# SECTION 5: PLOTS
###############################################################################

cat("--- Section 5: Generating plots ---\n")

# Plot 1: Hub count per module
hub_counts <- blood_hubs[, .N, by=Module]
hub_counts[, Module := factor(Module, levels=Module[order(-N)])]

p1 <- ggplot(hub_counts, aes(x=Module, y=N, fill=Module)) +
  geom_col(show.legend=FALSE) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title=sprintf("Hub Genes per Module (kME≥%.2f, kWithin top 25%%)",
                     KME_THRESH),
       x="Module", y="Number of Hub Genes")
ggsave(file.path(HUB_DIR, "plots", "step3_hub_counts.png"),
       p1, width=9, height=5, dpi=150)

# Plot 2: kME distribution per module (violin)
hub_df_nongrey <- hub_df[Module != "grey"]

kme_rows <- lapply(unique(hub_df_nongrey$Module), function(mod) {
  sub     <- as.data.frame(hub_df_nongrey[Module == mod])
  kME_col <- paste0("kME", mod)
  if (!kME_col %in% names(sub)) return(NULL)
  data.frame(Gene    = sub$Gene,
             Module  = mod,
             kME_own = abs(sub[[kME_col]]))
})
hub_df_nongrey_kme <- do.call(rbind, Filter(Negate(is.null), kme_rows))

p2 <- ggplot(hub_df_nongrey_kme[!is.na(hub_df_nongrey_kme$kME_own), ],
             aes(x=Module, y=kME_own, fill=Module)) +
  geom_violin(show.legend=FALSE, alpha=0.7, trim=TRUE) +
  geom_hline(yintercept=KME_THRESH, linetype="dashed",
             color="#E74C3C", linewidth=1) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Module Membership (|kME|) Distribution per Module",
       subtitle=sprintf("Red dashed = hub threshold (kME=%.2f)", KME_THRESH),
       x="Module", y="|kME| (Module Membership)")
ggsave(file.path(HUB_DIR, "plots", "step3_kME_distributions.png"),
       p2, width=10, height=6, dpi=150)
# Plot 3: kME vs kWithin scatter for hub genes
p3 <- ggplot(blood_hubs, aes(x=kME, y=kWithin, color=Module)) +
  geom_point(alpha=0.7, size=2) +
  theme_minimal(base_size=12) +
  labs(title="Hub Genes: kME vs Intramodular Connectivity",
       x="Module Membership (kME)", y="Intramodular Connectivity (kWithin)")
ggsave(file.path(HUB_DIR, "plots", "step3_kME_vs_kWithin.png"),
       p3, width=9, height=6, dpi=150)

###############################################################################
# SECTION 6: SAVE
###############################################################################

cat("--- Section 6: Saving hub genes ---\n")

blood_hubs <- blood_hubs[order(Module, -kME)]
fwrite(blood_hubs, file.path(HUB_DIR, "blood_hub_genes.csv"))

cat(sprintf("Saved: blood_hub_genes.csv (%d genes)\n\n",
            nrow(blood_hubs)))

###############################################################################
# FINAL SUMMARY
###############################################################################

cat("=======================================================================\n")
cat("  STEP 3 SUMMARY\n")
cat("=======================================================================\n")
cat(sprintf("  kME threshold:    %.2f\n", KME_THRESH))
cat(sprintf("  kWithin pctile:   %.0f%%\n", KWITHIN_PCTL*100))
cat(sprintf("  Total hub genes:  %d\n", nrow(blood_hubs)))
cat(sprintf("  Modules with hubs: %d\n", length(hub_list)))
cat("\n  Hubs per module:\n")
for(mod in names(hub_list)) {
  cat(sprintf("    %-15s: %d\n", mod, nrow(hub_list[[mod]])))
}
cat("=======================================================================\n")
cat("=== STEP 3 COMPLETED ===\n")