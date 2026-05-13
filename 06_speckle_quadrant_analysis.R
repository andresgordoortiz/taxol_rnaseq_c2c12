#!/usr/bin/env Rscript
# ============================================================================
# 06_speckle_quadrant_analysis.R
# Quadrant-based analysis of splicing events: speckle features & taxol effect
#
# PURPOSE:
#   1. Build the Differentiation scatter (DMSO vs Taxol context)
#   2. Classify events into quadrants (Blocked by Taxol, Taxol-enabled, etc.)
#   3. For each quadrant, pull Matt exon features
#   4. Compare features to nuclear speckle-associated gene properties
#   5. Output comprehensive statistics
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat(strrep("=", 70), "\n")
cat("06_speckle_quadrant_analysis.R\n")
cat(strrep("=", 70), "\n\n")

setwd("/Users/andres.ortiz/Projects/CRG/taxol_rnaseq_c2c12")

RESULTS_DIR <- file.path(getwd(), "results")

# ============================================================================
# 1. LOAD SPLICING FDR RESULTS
# ============================================================================

cat("--- Loading splicing results ---\n")
diff_dmso  <- read.csv(file.path(RESULTS_DIR, "Differentiation_DMSO.csv"))
diff_taxol <- read.csv(file.path(RESULTS_DIR, "Differentiation_Taxol.csv"))
myo_taxol  <- read.csv(file.path(RESULTS_DIR, "Myoblast_Taxol_vs_DMSO.csv"))
tube_taxol <- read.csv(file.path(RESULTS_DIR, "Myotube_Taxol_vs_DMSO.csv"))

cat(sprintf("  Diff DMSO:  %d events\n", nrow(diff_dmso)))
cat(sprintf("  Diff Taxol: %d events\n", nrow(diff_taxol)))
cat(sprintf("  Myo Taxol:  %d events\n", nrow(myo_taxol)))
cat(sprintf("  Tube Taxol: %d events\n", nrow(tube_taxol)))

# ============================================================================
# 2. LOAD MATT FEATURE FILES
# ============================================================================

cat("\n--- Loading Matt feature files ---\n")
matt_diff_dmso  <- read_tsv("Differentiation_DMSO_exons_with_efeatures.tab", show_col_types = FALSE)
matt_diff_taxol <- read_tsv("Differentiation_Taxol_exons_with_efeatures.tab", show_col_types = FALSE)
matt_myo_taxol  <- read_tsv("Myoblast_Taxol_vs_DMSO_exons_with_efeatures.tab", show_col_types = FALSE)
matt_tube_taxol <- read_tsv("Myotube_Taxol_vs_DMSO_exons_with_efeatures.tab", show_col_types = FALSE)

# ============================================================================
# 3. LOAD NUCLEAR SPECKLE GENES
# ============================================================================

cat("\n--- Loading nuclear speckle gene list ---\n")
speckle_genes <- scan("nuclear_speckle_associated_genes.txt", what = "character", quiet = TRUE)
cat(sprintf("  Speckle genes: %d\n", length(speckle_genes)))

# ============================================================================
# 4. LOAD DESEQ2 RESULTS
# ============================================================================

cat("\n--- Loading DESeq2 results ---\n")
deseq_myo  <- read.csv(file.path(RESULTS_DIR, "deseq2_Myoblast_Taxol_vs_DMSO.csv"))
deseq_tube <- read.csv(file.path(RESULTS_DIR, "deseq2_Myotube_Taxol_vs_DMSO.csv"))
deseq_diff_dmso  <- read.csv(file.path(RESULTS_DIR, "deseq2_Differentiation_DMSO.csv"))
deseq_diff_taxol <- read.csv(file.path(RESULTS_DIR, "deseq2_Differentiation_Taxol.csv"))

# ============================================================================
# 5. DEFINE THRESHOLDS & QUADRANT CLASSIFICATION
# ============================================================================

FDR_THRESH  <- 0.05
DPSI_THRESH <- 0.1

is_sig <- function(fdr, dpsi) {
  !is.na(fdr) & fdr <= FDR_THRESH & !is.na(dpsi) & abs(dpsi) >= DPSI_THRESH
}

# ============================================================================
# 5A. DIFFERENTIATION SCATTER: DMSO vs TAXOL CONTEXT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("DIFFERENTIATION SCATTER: DMSO vs TAXOL\n")
cat(strrep("=", 70), "\n\n")

dd <- diff_dmso %>% select(EVENT, GENE, dpsi_dmso = deltapsi, FDR_dmso = FDR)
dt <- diff_taxol %>% select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)

diff_scatter <- inner_join(dd, dt, by = "EVENT") %>%
  mutate(
    sig_dmso  = is_sig(FDR_dmso,  dpsi_dmso),
    sig_taxol = is_sig(FDR_taxol, dpsi_taxol),
    quadrant = case_when(
      sig_dmso & sig_taxol & sign(dpsi_dmso) == sign(dpsi_taxol) ~ "Both_concordant",
      sig_dmso & sig_taxol & sign(dpsi_dmso) != sign(dpsi_taxol) ~ "Reversed",
      sig_dmso & !sig_taxol ~ "Blocked_by_Taxol",
      !sig_dmso & sig_taxol ~ "Taxol_enabled",
      TRUE ~ "NS"
    ),
    # For blocked events: was the event skipped or included during differentiation?
    diff_direction = case_when(
      dpsi_dmso > 0 ~ "Included_in_diff",
      dpsi_dmso < 0 ~ "Skipped_in_diff",
      TRUE ~ "No_change"
    )
  )

cat("Differentiation scatter quadrant counts:\n")
print(table(diff_scatter$quadrant))
cat("\n")

# Blocked events detail: direction of change
blocked <- diff_scatter %>% filter(quadrant == "Blocked_by_Taxol")
cat(sprintf("Blocked by Taxol: %d events\n", nrow(blocked)))
cat("  Direction in DMSO differentiation:\n")
print(table(blocked$diff_direction))

both_conc <- diff_scatter %>% filter(quadrant == "Both_concordant")
cat(sprintf("\nBoth concordant: %d events\n", nrow(both_conc)))

enabled <- diff_scatter %>% filter(quadrant == "Taxol_enabled")
cat(sprintf("Taxol-enabled: %d events\n", nrow(enabled)))

reversed <- diff_scatter %>% filter(quadrant == "Reversed")
cat(sprintf("Reversed: %d events\n", nrow(reversed)))

# ============================================================================
# 5B. TAXOL SCATTER: MYOBLAST vs MYOTUBE
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TAXOL SCATTER: MYOBLAST vs MYOTUBE\n")
cat(strrep("=", 70), "\n\n")

tm <- myo_taxol %>% select(EVENT, GENE, dpsi_myo = deltapsi, FDR_myo = FDR)
tt <- tube_taxol %>% select(EVENT, dpsi_tube = deltapsi, FDR_tube = FDR)

taxol_scatter <- inner_join(tm, tt, by = "EVENT") %>%
  mutate(
    sig_myo  = is_sig(FDR_myo, dpsi_myo),
    sig_tube = is_sig(FDR_tube, dpsi_tube),
    quadrant = case_when(
      sig_myo & sig_tube & sign(dpsi_myo) == sign(dpsi_tube) ~ "Both_concordant",
      sig_myo & sig_tube & sign(dpsi_myo) != sign(dpsi_tube) ~ "Discordant",
      sig_myo & !sig_tube ~ "Myoblast_only",
      !sig_myo & sig_tube ~ "Myotube_only",
      TRUE ~ "NS"
    )
  )

cat("Taxol scatter quadrant counts:\n")
print(table(taxol_scatter$quadrant))

# ============================================================================
# 6. GENE-LEVEL SPECKLE ENRICHMENT PER QUADRANT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SPECKLE GENE ENRICHMENT PER QUADRANT (DIFFERENTIATION)\n")
cat(strrep("=", 70), "\n\n")

# Note: speckle genes are human symbols; our data is mouse (GENE column has mouse symbols)
# Need to check if there's mapping or if gene names are conserved enough

# Check overlap of gene names
all_genes_in_data <- unique(diff_scatter$GENE)
speckle_overlap <- intersect(toupper(all_genes_in_data), toupper(speckle_genes))
cat(sprintf("Total genes in diff scatter: %d\n", length(all_genes_in_data)))
cat(sprintf("Speckle genes (human list): %d\n", length(speckle_genes)))
cat(sprintf("Overlap (case-insensitive): %d\n", length(speckle_overlap)))

# Also try direct match
speckle_direct <- intersect(all_genes_in_data, speckle_genes)
cat(sprintf("Direct overlap: %d\n", length(speckle_direct)))

# Use upper-case matching since mouse symbols may differ in case
diff_scatter$GENE_upper <- toupper(diff_scatter$GENE)
speckle_upper <- toupper(speckle_genes)

# Background: all genes in the scatter
bg_genes <- unique(diff_scatter$GENE_upper)
n_bg <- length(bg_genes)
n_speckle_in_bg <- sum(bg_genes %in% speckle_upper)

cat(sprintf("\nBackground genes: %d\n", n_bg))
cat(sprintf("Speckle genes in background: %d (%.1f%%)\n",
            n_speckle_in_bg, 100*n_speckle_in_bg/n_bg))

# Enrichment per quadrant
for (q in c("Both_concordant", "Blocked_by_Taxol", "Taxol_enabled", "Reversed")) {
  q_genes <- unique(diff_scatter$GENE_upper[diff_scatter$quadrant == q])
  n_q <- length(q_genes)
  n_speckle_q <- sum(q_genes %in% speckle_upper)

  if (n_q > 0 && n_speckle_in_bg > 0) {
    # Fisher's exact test
    n_not_q_speckle <- n_speckle_in_bg - n_speckle_q
    n_q_not_speckle <- n_q - n_speckle_q
    n_neither <- n_bg - n_q - n_speckle_in_bg + n_speckle_q

    # Clamp to avoid negative values
    n_neither <- max(0, n_neither)

    cont_mat <- matrix(c(n_speckle_q, n_q_not_speckle,
                          n_not_q_speckle, n_neither), nrow = 2)
    ft <- tryCatch(fisher.test(cont_mat), error = function(e) list(p.value = NA, estimate = NA))

    cat(sprintf("\n  %s: %d genes, %d speckle (%.1f%%, OR=%.2f, p=%.4f)\n",
                q, n_q, n_speckle_q, 100*n_speckle_q/max(n_q,1),
                ifelse(is.na(ft$estimate), 0, ft$estimate),
                ifelse(is.na(ft$p.value), 1, ft$p.value)))
  } else {
    cat(sprintf("\n  %s: %d genes, %d speckle\n", q, n_q, n_speckle_q))
  }
}

# Also for taxol scatter
cat("\n\nSPECKLE GENE ENRICHMENT PER QUADRANT (TAXOL SCATTER)\n")
taxol_scatter_g <- taxol_scatter
taxol_scatter_g$GENE_upper <- toupper(taxol_scatter_g$GENE)

for (q in c("Both_concordant", "Myoblast_only", "Myotube_only", "Discordant")) {
  q_genes <- unique(taxol_scatter_g$GENE_upper[taxol_scatter_g$quadrant == q])
  n_q <- length(q_genes)
  n_speckle_q <- sum(q_genes %in% speckle_upper)
  cat(sprintf("  %s: %d genes, %d speckle (%.1f%%)\n",
              q, n_q, n_speckle_q, 100*n_speckle_q/max(n_q,1)))
}

# ============================================================================
# 7. MATT FEATURE ANALYSIS PER QUADRANT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("MATT FEATURES PER DIFFERENTIATION QUADRANT\n")
cat(strrep("=", 70), "\n\n")

# Key features à la speckle paper: exon length, GC content, intron length,
# splice site strength, branch point scores
key_features <- c(
  "EXON_LENGTH", "EXON_GCC",
  "UPINTRON_MEDIANLENGTH", "DOINTRON_MEDIANLENGTH",
  "UPINTRON_GCC", "DOINTRON_GCC",
  "MAXENTSCR_HSAMODEL_3SS", "MAXENTSCR_HSAMODEL_5SS",
  "MAXENTSCR_MMUMODEL_3SS", "MAXENTSCR_MMUMODEL_5SS",
  "UPEXON_GCC", "DOEXON_GCC",
  "MEDIAN_EXON_NUMBER", "MEDIAN_TR_LENGTH",
  "RATIO_UPINTRON_EXON_LENGTH", "RATIO_DOINTRON_EXON_LENGTH",
  "BPSCORE_MAXBP_UPINTRON", "BPSCORE_MAXBP_DOINTRON",
  "NUM_PREDICTED_BPS_UPINTRON", "NUM_PREDICTED_BPS_DOINTRON",
  "POLYPYRITRAC_LEN_MAXBP_UPINTRON", "POLYPYRITRAC_LEN_MAXBP_DOINTRON",
  "GCC_3SS_20INT10EX", "GCC_5SS_20INT10EX"
)

# Use Matt Diff_DMSO file — it has GROUP (up/down/ndiff) and features
# Merge the quadrant info with Matt data by EVENT
matt_with_quad <- matt_diff_dmso %>%
  left_join(
    diff_scatter %>% select(EVENT, quadrant, diff_direction, sig_dmso, sig_taxol,
                            dpsi_dmso, dpsi_taxol),
    by = "EVENT"
  )

# Check available features
avail_feats <- intersect(key_features, colnames(matt_with_quad))
cat(sprintf("Available key features: %d / %d\n", length(avail_feats), length(key_features)))
cat(paste(avail_feats, collapse = ", "), "\n\n")

# For each quadrant, compute median of key features and compare to NS background
quad_stats <- list()

for (q in c("Both_concordant", "Blocked_by_Taxol", "Taxol_enabled", "Reversed")) {
  q_events <- matt_with_quad %>% filter(quadrant == q)
  ns_events <- matt_with_quad %>% filter(quadrant == "NS" | is.na(quadrant))

  cat(sprintf("\n--- %s: %d events ---\n", q, nrow(q_events)))

  for (feat in avail_feats) {
    vals_q  <- q_events[[feat]]
    vals_ns <- ns_events[[feat]]
    vals_q  <- vals_q[!is.na(vals_q)]
    vals_ns <- vals_ns[!is.na(vals_ns)]

    if (length(vals_q) >= 3 && length(vals_ns) >= 3) {
      wt <- tryCatch(
        wilcox.test(vals_q, vals_ns, exact = FALSE),
        error = function(e) list(p.value = NA, statistic = NA)
      )

      quad_stats[[paste(q, feat, sep = "__")]] <- data.frame(
        quadrant = q,
        feature  = feat,
        n_quad   = length(vals_q),
        n_ns     = length(vals_ns),
        median_quad = median(vals_q),
        median_ns   = median(vals_ns),
        diff_medians = median(vals_q) - median(vals_ns),
        p.value = wt$p.value,
        stringsAsFactors = FALSE
      )

      if (!is.na(wt$p.value) && wt$p.value < 0.05) {
        cat(sprintf("  * %s: med_q=%.3f, med_ns=%.3f, diff=%.3f, p=%.4f\n",
                    feat, median(vals_q), median(vals_ns),
                    median(vals_q) - median(vals_ns), wt$p.value))
      }
    }
  }
}

quad_stats_df <- bind_rows(quad_stats)
quad_stats_df <- quad_stats_df %>%
  group_by(quadrant) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(quadrant, p.value)

write_csv(quad_stats_df, file.path(RESULTS_DIR, "quadrant_feature_statistics.csv"))
cat(sprintf("\nSaved: quadrant_feature_statistics.csv (%d tests)\n", nrow(quad_stats_df)))

# ============================================================================
# 7B. BLOCKED vs BOTH CONCORDANT: Direct feature comparison
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("BLOCKED vs BOTH_CONCORDANT: Direct feature comparison\n")
cat(strrep("=", 70), "\n\n")

blocked_matt <- matt_with_quad %>% filter(quadrant == "Blocked_by_Taxol")
both_matt    <- matt_with_quad %>% filter(quadrant == "Both_concordant")

for (feat in avail_feats) {
  v_blocked <- blocked_matt[[feat]][!is.na(blocked_matt[[feat]])]
  v_both    <- both_matt[[feat]][!is.na(both_matt[[feat]])]

  if (length(v_blocked) >= 3 && length(v_both) >= 3) {
    wt <- tryCatch(wilcox.test(v_blocked, v_both, exact = FALSE),
                   error = function(e) list(p.value = NA))
    if (!is.na(wt$p.value) && wt$p.value < 0.1) {
      cat(sprintf("  %s: blocked=%.3f, both=%.3f, p=%.4f\n",
                  feat, median(v_blocked), median(v_both), wt$p.value))
    }
  }
}

# ============================================================================
# 8. SIGNATURE OF "SPECKLE EXONS" (from paper)
#    According to Małszycki et al., speckle-dependent exons have:
#    - Short exons
#    - Weak splice sites (low MaxEnt scores)
#    - Located in GC-rich isochores
#    - Flanked by short introns (relative to non-speckle)
#    Compare this signature to blocked-by-Taxol exons
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SPECKLE EXON SIGNATURE vs QUADRANT FEATURES\n")
cat(strrep("=", 70), "\n\n")

# Define "speckle signature" features and expected directions
# (from Małszycki et al.):
# Speckle exons tend to be:
# - Shorter exon length
# - Higher GC content (GC-rich isochores)
# - Weaker splice sites (lower MaxEnt scores)

speckle_signature <- data.frame(
  feature = c("EXON_LENGTH", "EXON_GCC",
              "MAXENTSCR_HSAMODEL_3SS", "MAXENTSCR_HSAMODEL_5SS",
              "UPINTRON_MEDIANLENGTH", "DOINTRON_MEDIANLENGTH",
              "UPINTRON_GCC", "DOINTRON_GCC"),
  expected_in_speckle = c("lower", "higher",
                           "lower", "lower",
                           "shorter", "shorter",
                           "higher", "higher"),
  stringsAsFactors = FALSE
)

cat("Speckle exon signature (from Małszycki et al.):\n")
print(speckle_signature, row.names = FALSE)

cat("\nComparing blocked-by-Taxol exons to NS (background):\n")
blocked_vs_ns <- quad_stats_df %>% filter(quadrant == "Blocked_by_Taxol")
for (i in 1:nrow(speckle_signature)) {
  feat <- speckle_signature$feature[i]
  row <- blocked_vs_ns %>% filter(feature == feat)
  if (nrow(row) > 0) {
    direction <- ifelse(row$diff_medians > 0, "higher", "lower")
    matches <- direction == speckle_signature$expected_in_speckle[i] ||
               (speckle_signature$expected_in_speckle[i] == "shorter" && direction == "lower")
    cat(sprintf("  %s: %s in blocked (diff=%.3f, p=%.4f) — expected %s — %s\n",
                feat, direction, row$diff_medians, row$p.value,
                speckle_signature$expected_in_speckle[i],
                ifelse(matches, "MATCHES", "does not match")))
  }
}

# Same for Taxol-enabled
cat("\nComparing Taxol-enabled exons to NS:\n")
enabled_vs_ns <- quad_stats_df %>% filter(quadrant == "Taxol_enabled")
for (i in 1:nrow(speckle_signature)) {
  feat <- speckle_signature$feature[i]
  row <- enabled_vs_ns %>% filter(feature == feat)
  if (nrow(row) > 0) {
    direction <- ifelse(row$diff_medians > 0, "higher", "lower")
    cat(sprintf("  %s: %s in enabled (diff=%.3f, p=%.4f)\n",
                feat, direction, row$diff_medians, row$p.value))
  }
}

# ============================================================================
# 9. EVENT COUNTS & DIRECTIONS IN EACH QUADRANT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("DETAILED QUADRANT BREAKDOWN\n")
cat(strrep("=", 70), "\n\n")

# Blocked-by-Taxol: events that change during normal differentiation but NOT
# when Taxol is present
blocked_detail <- diff_scatter %>%
  filter(quadrant == "Blocked_by_Taxol") %>%
  mutate(
    direction = case_when(
      dpsi_dmso > DPSI_THRESH  ~ "Included in differentiation (blocked by Taxol)",
      dpsi_dmso < -DPSI_THRESH ~ "Skipped in differentiation (blocked by Taxol)",
      TRUE ~ "other"
    )
  )

cat("BLOCKED BY TAXOL - direction breakdown:\n")
print(table(blocked_detail$direction))

# Speckle gene overlap per direction
cat("\nSpeckle gene overlap per direction (blocked):\n")
for (d in unique(blocked_detail$direction)) {
  d_genes <- unique(toupper(blocked_detail$GENE[blocked_detail$direction == d]))
  n_speckle <- sum(d_genes %in% speckle_upper)
  cat(sprintf("  %s: %d genes, %d speckle (%.1f%%)\n",
              d, length(d_genes), n_speckle, 100*n_speckle/max(length(d_genes),1)))
}

# Both concordant: detail
both_detail <- diff_scatter %>%
  filter(quadrant == "Both_concordant") %>%
  mutate(
    direction = case_when(
      dpsi_dmso > DPSI_THRESH  ~ "Included in differentiation (preserved)",
      dpsi_dmso < -DPSI_THRESH ~ "Skipped in differentiation (preserved)",
      TRUE ~ "other"
    )
  )

cat("\nBOTH CONCORDANT - direction breakdown:\n")
print(table(both_detail$direction))

# ============================================================================
# 10. MATT FEATURES: BLOCKED-SKIPPED vs BLOCKED-INCLUDED
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("BLOCKED EXONS: SKIPPED vs INCLUDED sub-groups\n")
cat(strrep("=", 70), "\n\n")

blocked_included_events <- diff_scatter$EVENT[diff_scatter$quadrant == "Blocked_by_Taxol" &
                                                diff_scatter$dpsi_dmso > DPSI_THRESH]
blocked_skipped_events  <- diff_scatter$EVENT[diff_scatter$quadrant == "Blocked_by_Taxol" &
                                                diff_scatter$dpsi_dmso < -DPSI_THRESH]

cat(sprintf("Blocked-Included: %d events\n", length(blocked_included_events)))
cat(sprintf("Blocked-Skipped:  %d events\n", length(blocked_skipped_events)))

matt_blocked_inc <- matt_with_quad %>% filter(EVENT %in% blocked_included_events)
matt_blocked_skp <- matt_with_quad %>% filter(EVENT %in% blocked_skipped_events)

cat("\nFeature comparison (blocked-included vs blocked-skipped):\n")
for (feat in avail_feats) {
  v_inc <- matt_blocked_inc[[feat]][!is.na(matt_blocked_inc[[feat]])]
  v_skp <- matt_blocked_skp[[feat]][!is.na(matt_blocked_skp[[feat]])]

  if (length(v_inc) >= 3 && length(v_skp) >= 3) {
    wt <- tryCatch(wilcox.test(v_inc, v_skp, exact = FALSE),
                   error = function(e) list(p.value = NA))
    if (!is.na(wt$p.value) && wt$p.value < 0.1) {
      cat(sprintf("  %s: included=%.3f, skipped=%.3f, p=%.4f\n",
                  feat, median(v_inc), median(v_skp), wt$p.value))
    }
  }
}

# ============================================================================
# 11. SUMMARY TABLE OF ALL RELEVANT SPLICING EVENTS
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("OVERALL SIGNIFICANCE STATISTICS\n")
cat(strrep("=", 70), "\n\n")

# Counts per comparison
for (comp in list(
  list(name = "Myoblast Taxol vs DMSO", df = myo_taxol),
  list(name = "Myotube Taxol vs DMSO", df = tube_taxol),
  list(name = "Differentiation DMSO", df = diff_dmso),
  list(name = "Differentiation Taxol", df = diff_taxol)
)) {
  df <- comp$df
  n_sig <- sum(is_sig(df$FDR, df$deltapsi))
  n_up  <- sum(is_sig(df$FDR, df$deltapsi) & df$deltapsi > 0)
  n_down <- sum(is_sig(df$FDR, df$deltapsi) & df$deltapsi < 0)
  cat(sprintf("  %s: %d sig total (%d included, %d skipped)\n",
              comp$name, n_sig, n_up, n_down))
}

# ============================================================================
# 12. EXPRESSION OF BLOCKED GENES
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("EXPRESSION OF BLOCKED-BY-TAXOL GENES\n")
cat(strrep("=", 70), "\n\n")

blocked_genes <- unique(diff_scatter$GENE[diff_scatter$quadrant == "Blocked_by_Taxol"])
cat(sprintf("Blocked-by-Taxol genes: %d\n", length(blocked_genes)))

# Are these genes differentially expressed in any comparison?
for (comp in list(
  list(name = "Myoblast Taxol vs DMSO", df = deseq_myo),
  list(name = "Myotube Taxol vs DMSO", df = deseq_tube),
  list(name = "Differentiation DMSO", df = deseq_diff_dmso),
  list(name = "Differentiation Taxol", df = deseq_diff_taxol)
)) {
  de_genes <- comp$df %>%
    filter(!is.na(significance) & significance != "Not Significant") %>%
    pull(gene_symbol)
  overlap <- sum(blocked_genes %in% de_genes)
  cat(sprintf("  %s: %d / %d blocked genes are DE (%.1f%%)\n",
              comp$name, overlap, length(blocked_genes),
              100*overlap/max(length(blocked_genes),1)))
}

# ============================================================================
# 13. EXISTING MATT FEATURE STATISTICS (from script 04)
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("EXISTING MATT FEATURE STATISTICS\n")
cat(strrep("=", 70), "\n\n")

matt_stats <- read_csv(file.path(RESULTS_DIR, "matt_feature_statistics.csv"),
                        show_col_types = FALSE)
cat(sprintf("Total tests: %d\n", nrow(matt_stats)))
cat(sprintf("Significant: %d\n", sum(matt_stats$is_significant, na.rm = TRUE)))

# Show top features
top_feats <- matt_stats %>%
  filter(is_significant) %>%
  arrange(p.value) %>%
  head(30)

cat("\nTop 30 significant features:\n")
print(as.data.frame(top_feats %>% select(feature, comparison, treatment,
                                          median_dir, median_bg, diff, p.value)),
      row.names = FALSE)

# ============================================================================
# 14. MATT FEATURES: Speckle-relevant features by treatment
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SPECKLE-RELEVANT FEATURES BY TREATMENT\n")
cat(strrep("=", 70), "\n\n")

# Compare speckle-relevant features between Differentiation_DMSO (up/down)
# and how they differ when looking at Taxol context

speckle_feats <- c("EXON_LENGTH", "EXON_GCC",
                    "UPINTRON_MEDIANLENGTH", "DOINTRON_MEDIANLENGTH",
                    "MAXENTSCR_HSAMODEL_3SS", "MAXENTSCR_HSAMODEL_5SS",
                    "UPINTRON_GCC", "DOINTRON_GCC")

relevant_matt <- matt_stats %>%
  filter(feature %in% speckle_feats) %>%
  select(feature, comparison, treatment, median_dir, median_bg, diff, p.value, is_significant) %>%
  arrange(feature, treatment, comparison)

cat("Speckle-relevant features across all comparisons:\n")
print(as.data.frame(relevant_matt), row.names = FALSE)

# ============================================================================
# 15. FINAL COMPREHENSIVE SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("COMPREHENSIVE SUMMARY\n")
cat(strrep("=", 70), "\n\n")

# Overall numbers
n_diff_sig_dmso  <- sum(is_sig(diff_dmso$FDR, diff_dmso$deltapsi))
n_diff_sig_taxol <- sum(is_sig(diff_taxol$FDR, diff_taxol$deltapsi))
n_myo_sig  <- sum(is_sig(myo_taxol$FDR, myo_taxol$deltapsi))
n_tube_sig <- sum(is_sig(tube_taxol$FDR, tube_taxol$deltapsi))

cat(sprintf("Differentiation (DMSO): %d significant events\n", n_diff_sig_dmso))
cat(sprintf("Differentiation (Taxol): %d significant events\n", n_diff_sig_taxol))
cat(sprintf("  Taxol changes the differentiation programme: %d -> %d events\n",
            n_diff_sig_dmso, n_diff_sig_taxol))
cat(sprintf("Myoblast Taxol effect: %d significant events\n", n_myo_sig))
cat(sprintf("Myotube Taxol effect: %d significant events\n", n_tube_sig))

cat("\n--- Differentiation quadrants ---\n")
cat(sprintf("  Both concordant: %d events (differentiation splicing preserved)\n",
            sum(diff_scatter$quadrant == "Both_concordant")))
cat(sprintf("  Blocked by Taxol: %d events (differentiation splicing SUPPRESSED)\n",
            sum(diff_scatter$quadrant == "Blocked_by_Taxol")))
cat(sprintf("  Taxol-enabled: %d events (new splicing events under Taxol)\n",
            sum(diff_scatter$quadrant == "Taxol_enabled")))
cat(sprintf("  Reversed: %d events (opposite direction)\n",
            sum(diff_scatter$quadrant == "Reversed")))

cat("\nDone.\n")
