#!/usr/bin/env Rscript
# ============================================================================
# 07_intron_retention_quadrant_analysis.R
# Same quadrant/speckle analysis as 06, but specifically for intron retention
# events (MmuINT). Extracts intron features from EVENT_INFO and sequences.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(Biostrings)
})

cat(strrep("=", 70), "\n")
cat("07_intron_retention_quadrant_analysis.R\n")
cat(strrep("=", 70), "\n\n")

setwd("/Users/andres.ortiz/Projects/CRG/taxol_rnaseq_c2c12")
RESULTS_DIR <- file.path(getwd(), "results")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("--- Loading splicing results (RDS for complete data) ---\n")
diff_dmso  <- readRDS(file.path(RESULTS_DIR, "Differentiation_DMSO.rds"))
diff_taxol <- readRDS(file.path(RESULTS_DIR, "Differentiation_Taxol.rds"))
myo_taxol  <- readRDS(file.path(RESULTS_DIR, "Myoblast_Taxol_vs_DMSO.rds"))
tube_taxol <- readRDS(file.path(RESULTS_DIR, "Myotube_Taxol_vs_DMSO.rds"))

# Filter to IR events only
filter_ir <- function(df) df %>% filter(grepl("MmuINT", EVENT))

diff_dmso_ir  <- filter_ir(diff_dmso)
diff_taxol_ir <- filter_ir(diff_taxol)
myo_taxol_ir  <- filter_ir(myo_taxol)
tube_taxol_ir <- filter_ir(tube_taxol)

# Also keep EX-only for comparison
filter_ex <- function(df) df %>% filter(grepl("MmuEX", EVENT))
diff_dmso_ex  <- filter_ex(diff_dmso)
diff_taxol_ex <- filter_ex(diff_taxol)

cat(sprintf("  IR events: Diff_DMSO=%d, Diff_Taxol=%d, Myo=%d, Tube=%d\n",
            nrow(diff_dmso_ir), nrow(diff_taxol_ir), nrow(myo_taxol_ir), nrow(tube_taxol_ir)))
cat(sprintf("  EX events: Diff_DMSO=%d, Diff_Taxol=%d\n",
            nrow(diff_dmso_ex), nrow(diff_taxol_ex)))

# Load EVENT_INFO for intron features
cat("\n--- Loading EVENT_INFO ---\n")
event_info <- read.table("EVENT_INFO-mm10.tab", header = TRUE, sep = "\t",
                          quote = "", stringsAsFactors = FALSE)
ir_info <- event_info %>% filter(grepl("MmuINT", EVENT))
cat(sprintf("  IR events in EVENT_INFO: %d\n", nrow(ir_info)))

# Load speckle genes
speckle_genes <- scan("nuclear_speckle_associated_genes.txt", what = "character", quiet = TRUE)
speckle_upper <- toupper(speckle_genes)
cat(sprintf("  Speckle genes: %d\n", length(speckle_genes)))

# Load DESeq2
deseq_myo  <- read.csv(file.path(RESULTS_DIR, "deseq2_Myoblast_Taxol_vs_DMSO.csv"))
deseq_tube <- read.csv(file.path(RESULTS_DIR, "deseq2_Myotube_Taxol_vs_DMSO.csv"))
deseq_diff_d <- read.csv(file.path(RESULTS_DIR, "deseq2_Differentiation_DMSO.csv"))
deseq_diff_t <- read.csv(file.path(RESULTS_DIR, "deseq2_Differentiation_Taxol.csv"))

# ============================================================================
# 2. COMPUTE INTRON FEATURES FROM EVENT_INFO
# ============================================================================

cat("\n--- Computing intron features from sequences ---\n")

gc_content <- function(seq_str) {
  if (is.na(seq_str) || seq_str == "") return(NA_real_)
  seq_str <- toupper(seq_str)
  chars <- strsplit(seq_str, "")[[1]]
  gc <- sum(chars %in% c("G", "C"))
  total <- length(chars)
  if (total == 0) return(NA_real_)
  return(gc / total)
}

# Extract features for IR events
ir_features <- ir_info %>%
  mutate(
    INTRON_LENGTH = nchar(Seq_A),
    INTRON_GCC = sapply(Seq_A, gc_content),
    UPSTREAM_EXON_LENGTH = nchar(Seq_C1),
    DOWNSTREAM_EXON_LENGTH = nchar(Seq_C2),
    UPSTREAM_EXON_GCC = sapply(Seq_C1, gc_content),
    DOWNSTREAM_EXON_GCC = sapply(Seq_C2, gc_content),
    # 5' splice site (first 9nt of intron: 3 exonic + 6 intronic for GT-AG)
    SS5_SEQ = ifelse(nchar(Seq_A) >= 6, substr(Seq_A, 1, 6), NA_character_),
    # 3' splice site (last 20nt of intron)
    SS3_SEQ = ifelse(nchar(Seq_A) >= 20,
                     substr(Seq_A, nchar(Seq_A) - 19, nchar(Seq_A)),
                     NA_character_),
    # Polypyrimidine tract: count Y (C/T) in last 40nt of intron (before AG)
    PPT_REGION = ifelse(nchar(Seq_A) >= 40,
                        substr(Seq_A, nchar(Seq_A) - 39, nchar(Seq_A) - 2),
                        ifelse(nchar(Seq_A) >= 3,
                               substr(Seq_A, 1, nchar(Seq_A) - 2), NA_character_)),
    PPT_PYRIMIDINE_FRAC = sapply(PPT_REGION, function(s) {
      if (is.na(s) || s == "") return(NA_real_)
      chars <- strsplit(toupper(s), "")[[1]]
      sum(chars %in% c("C", "T")) / length(chars)
    }),
    # GC content of splice site regions
    SS5_GCC = sapply(SS5_SEQ, gc_content),
    SS3_GCC = sapply(SS3_SEQ, gc_content),
    # Intron type (GT-AG canonical?)
    DINUCLEOTIDE_5 = ifelse(nchar(Seq_A) >= 2, toupper(substr(Seq_A, 1, 2)), NA_character_),
    DINUCLEOTIDE_3 = ifelse(nchar(Seq_A) >= 2,
                            toupper(substr(Seq_A, nchar(Seq_A) - 1, nchar(Seq_A))),
                            NA_character_),
    IS_CANONICAL = (DINUCLEOTIDE_5 == "GT" & DINUCLEOTIDE_3 == "AG")
  ) %>%
  select(GENE, EVENT, INTRON_LENGTH, INTRON_GCC,
         UPSTREAM_EXON_LENGTH, DOWNSTREAM_EXON_LENGTH,
         UPSTREAM_EXON_GCC, DOWNSTREAM_EXON_GCC,
         PPT_PYRIMIDINE_FRAC, SS5_GCC, SS3_GCC,
         IS_CANONICAL, DINUCLEOTIDE_5, DINUCLEOTIDE_3)

cat(sprintf("  Computed features for %d IR events\n", nrow(ir_features)))
cat(sprintf("  Canonical GT-AG: %d (%.1f%%)\n",
            sum(ir_features$IS_CANONICAL, na.rm = TRUE),
            100 * mean(ir_features$IS_CANONICAL, na.rm = TRUE)))

# Save precomputed features for downstream scripts (09, 10)
write.csv(ir_features, file.path(RESULTS_DIR, "ir_features_precomputed.csv"), row.names = FALSE)
cat(sprintf("  Saved ir_features_precomputed.csv (%d rows)\n", nrow(ir_features)))

# ============================================================================
# 3. THRESHOLDS & SIGNIFICANCE
# ============================================================================

FDR_THRESH  <- 0.05
DPSI_THRESH <- 0.1  # Note: for IR, deltapsi is deltaPIR

is_sig <- function(fdr, dpsi) {
  !is.na(fdr) & fdr <= FDR_THRESH & !is.na(dpsi) & abs(dpsi) >= DPSI_THRESH
}

# ============================================================================
# 4. IR SIGNIFICANCE COUNTS
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("IR EVENT SIGNIFICANCE COUNTS\n")
cat(strrep("=", 70), "\n\n")

for (x in list(
  list("Myoblast Taxol vs DMSO (IR)", myo_taxol_ir),
  list("Myotube Taxol vs DMSO (IR)", tube_taxol_ir),
  list("Differentiation DMSO (IR)", diff_dmso_ir),
  list("Differentiation Taxol (IR)", diff_taxol_ir)
)) {
  df <- x[[2]]
  n_sig <- sum(is_sig(df$FDR, df$deltapsi))
  n_up  <- sum(is_sig(df$FDR, df$deltapsi) & df$deltapsi > 0)
  n_dn  <- sum(is_sig(df$FDR, df$deltapsi) & df$deltapsi < 0)
  cat(sprintf("  %s: %d sig (%d retention_up, %d retention_down)\n",
              x[[1]], n_sig, n_up, n_dn))
}

# Also for EX events
cat("\n  For comparison (EX events):\n")
for (x in list(
  list("Myoblast Taxol vs DMSO (EX)", filter_ex(myo_taxol)),
  list("Myotube Taxol vs DMSO (EX)", filter_ex(tube_taxol)),
  list("Differentiation DMSO (EX)", diff_dmso_ex),
  list("Differentiation Taxol (EX)", diff_taxol_ex)
)) {
  df <- x[[2]]
  n_sig <- sum(is_sig(df$FDR, df$deltapsi))
  n_up  <- sum(is_sig(df$FDR, df$deltapsi) & df$deltapsi > 0)
  n_dn  <- sum(is_sig(df$FDR, df$deltapsi) & df$deltapsi < 0)
  cat(sprintf("  %s: %d sig (%d included, %d skipped)\n",
              x[[1]], n_sig, n_up, n_dn))
}

# ============================================================================
# 5. DIFFERENTIATION SCATTER — IR EVENTS ONLY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("DIFFERENTIATION SCATTER — IR EVENTS ONLY\n")
cat(strrep("=", 70), "\n\n")

dd_ir <- diff_dmso_ir %>% select(EVENT, GENE, dpir_dmso = deltapsi, FDR_dmso = FDR)
dt_ir <- diff_taxol_ir %>% select(EVENT, dpir_taxol = deltapsi, FDR_taxol = FDR)

ir_scatter <- inner_join(dd_ir, dt_ir, by = "EVENT") %>%
  mutate(
    sig_dmso  = is_sig(FDR_dmso, dpir_dmso),
    sig_taxol = is_sig(FDR_taxol, dpir_taxol),
    quadrant = case_when(
      sig_dmso & sig_taxol & sign(dpir_dmso) == sign(dpir_taxol) ~ "Both_concordant",
      sig_dmso & sig_taxol & sign(dpir_dmso) != sign(dpir_taxol) ~ "Reversed",
      sig_dmso & !sig_taxol ~ "Blocked_by_Taxol",
      !sig_dmso & sig_taxol ~ "Taxol_enabled",
      TRUE ~ "NS"
    ),
    direction = case_when(
      dpir_dmso > 0 ~ "Retention_increased",
      dpir_dmso < 0 ~ "Retention_decreased",
      TRUE ~ "No_change"
    )
  )

cat("IR Differentiation scatter quadrant counts:\n")
print(table(ir_scatter$quadrant))

blocked_ir <- ir_scatter %>% filter(quadrant == "Blocked_by_Taxol")
cat(sprintf("\nBlocked by Taxol (IR): %d events\n", nrow(blocked_ir)))
cat("  Direction in DMSO differentiation:\n")
print(table(blocked_ir$direction))

both_ir <- ir_scatter %>% filter(quadrant == "Both_concordant")
cat(sprintf("\nBoth concordant (IR): %d events\n", nrow(both_ir)))
cat("  Direction:\n")
print(table(both_ir$direction))

enabled_ir <- ir_scatter %>% filter(quadrant == "Taxol_enabled")
cat(sprintf("\nTaxol-enabled (IR): %d events\n", nrow(enabled_ir)))

# ============================================================================
# 5B. ALSO DO EX-ONLY SCATTER FOR COMPARISON
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("DIFFERENTIATION SCATTER — EX EVENTS ONLY\n")
cat(strrep("=", 70), "\n\n")

dd_ex <- diff_dmso_ex %>% select(EVENT, GENE, dpsi_dmso = deltapsi, FDR_dmso = FDR)
dt_ex <- diff_taxol_ex %>% select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)

ex_scatter <- inner_join(dd_ex, dt_ex, by = "EVENT") %>%
  mutate(
    sig_dmso  = is_sig(FDR_dmso, dpsi_dmso),
    sig_taxol = is_sig(FDR_taxol, dpsi_taxol),
    quadrant = case_when(
      sig_dmso & sig_taxol & sign(dpsi_dmso) == sign(dpsi_taxol) ~ "Both_concordant",
      sig_dmso & sig_taxol & sign(dpsi_dmso) != sign(dpsi_taxol) ~ "Reversed",
      sig_dmso & !sig_taxol ~ "Blocked_by_Taxol",
      !sig_dmso & sig_taxol ~ "Taxol_enabled",
      TRUE ~ "NS"
    )
  )

cat("EX Differentiation scatter quadrant counts:\n")
print(table(ex_scatter$quadrant))

# ============================================================================
# 6. SPECKLE ENRICHMENT — IR QUADRANTS
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SPECKLE GENE ENRICHMENT — IR QUADRANTS\n")
cat(strrep("=", 70), "\n\n")

ir_scatter$GENE_upper <- toupper(ir_scatter$GENE)

bg_genes_ir <- unique(ir_scatter$GENE_upper)
n_bg_ir <- length(bg_genes_ir)
n_speckle_bg_ir <- sum(bg_genes_ir %in% speckle_upper)

cat(sprintf("IR Background: %d genes, %d speckle (%.1f%%)\n",
            n_bg_ir, n_speckle_bg_ir, 100 * n_speckle_bg_ir / n_bg_ir))

for (q in c("Both_concordant", "Blocked_by_Taxol", "Taxol_enabled", "Reversed")) {
  q_genes <- unique(ir_scatter$GENE_upper[ir_scatter$quadrant == q])
  n_q <- length(q_genes)
  n_sp <- sum(q_genes %in% speckle_upper)

  if (n_q > 2 && n_speckle_bg_ir > 0) {
    n_neither <- max(0, n_bg_ir - n_q - n_speckle_bg_ir + n_sp)
    cont <- matrix(c(n_sp, n_q - n_sp, n_speckle_bg_ir - n_sp, n_neither), nrow = 2)
    ft <- fisher.test(cont)
    cat(sprintf("  %s: %d genes, %d speckle (%.1f%%), OR=%.2f, p=%.2e\n",
                q, n_q, n_sp, 100 * n_sp / max(n_q, 1), ft$estimate, ft$p.value))
  } else {
    cat(sprintf("  %s: %d genes, %d speckle (%.1f%%)\n",
                q, n_q, n_sp, 100 * n_sp / max(n_q, 1)))
  }
}

# Sub-groups within blocked
blocked_ret_up <- unique(toupper(blocked_ir$GENE[blocked_ir$dpir_dmso > DPSI_THRESH]))
blocked_ret_dn <- unique(toupper(blocked_ir$GENE[blocked_ir$dpir_dmso < -DPSI_THRESH]))
cat(sprintf("\nBlocked-Retention_up: %d genes, %d speckle (%.1f%%)\n",
            length(blocked_ret_up), sum(blocked_ret_up %in% speckle_upper),
            100 * sum(blocked_ret_up %in% speckle_upper) / max(length(blocked_ret_up), 1)))
cat(sprintf("Blocked-Retention_down: %d genes, %d speckle (%.1f%%)\n",
            length(blocked_ret_dn), sum(blocked_ret_dn %in% speckle_upper),
            100 * sum(blocked_ret_dn %in% speckle_upper) / max(length(blocked_ret_dn), 1)))

# ============================================================================
# 7. INTRON FEATURES PER QUADRANT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("INTRON FEATURES PER DIFFERENTIATION QUADRANT (IR)\n")
cat(strrep("=", 70), "\n\n")

# Merge features with scatter
ir_with_feat <- ir_scatter %>%
  left_join(ir_features, by = "EVENT")

feature_cols <- c("INTRON_LENGTH", "INTRON_GCC",
                   "UPSTREAM_EXON_LENGTH", "DOWNSTREAM_EXON_LENGTH",
                   "UPSTREAM_EXON_GCC", "DOWNSTREAM_EXON_GCC",
                   "PPT_PYRIMIDINE_FRAC", "SS5_GCC", "SS3_GCC")

avail_feats <- intersect(feature_cols, colnames(ir_with_feat))
cat(sprintf("Available features: %d\n\n", length(avail_feats)))

ir_quad_stats <- list()

for (q in c("Both_concordant", "Blocked_by_Taxol", "Taxol_enabled")) {
  q_events <- ir_with_feat %>% filter(quadrant == q)
  ns_events <- ir_with_feat %>% filter(quadrant == "NS")

  cat(sprintf("--- %s: %d events (vs %d NS) ---\n", q, nrow(q_events), nrow(ns_events)))

  for (feat in avail_feats) {
    vals_q  <- q_events[[feat]][!is.na(q_events[[feat]])]
    vals_ns <- ns_events[[feat]][!is.na(ns_events[[feat]])]

    if (length(vals_q) >= 3 && length(vals_ns) >= 3) {
      wt <- tryCatch(
        wilcox.test(vals_q, vals_ns, exact = FALSE),
        error = function(e) list(p.value = NA)
      )

      ir_quad_stats[[paste(q, feat, sep = "__")]] <- data.frame(
        quadrant = q, feature = feat,
        n_quad = length(vals_q), n_ns = length(vals_ns),
        median_quad = median(vals_q), median_ns = median(vals_ns),
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
  cat("\n")
}

ir_quad_stats_df <- bind_rows(ir_quad_stats)
ir_quad_stats_df <- ir_quad_stats_df %>%
  group_by(quadrant) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(quadrant, p.value)

write_csv(ir_quad_stats_df, file.path(RESULTS_DIR, "ir_quadrant_feature_statistics.csv"))
cat(sprintf("Saved: ir_quadrant_feature_statistics.csv (%d tests)\n", nrow(ir_quad_stats_df)))

# ============================================================================
# 8. BLOCKED vs CONCORDANT — DIRECT IR FEATURE COMPARISON
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("BLOCKED vs CONCORDANT — IR FEATURES\n")
cat(strrep("=", 70), "\n\n")

blocked_ir_feat <- ir_with_feat %>% filter(quadrant == "Blocked_by_Taxol")
both_ir_feat    <- ir_with_feat %>% filter(quadrant == "Both_concordant")

cat(sprintf("Blocked IR: %d events, Concordant IR: %d events\n",
            nrow(blocked_ir_feat), nrow(both_ir_feat)))

for (feat in avail_feats) {
  v_blocked <- blocked_ir_feat[[feat]][!is.na(blocked_ir_feat[[feat]])]
  v_both    <- both_ir_feat[[feat]][!is.na(both_ir_feat[[feat]])]

  if (length(v_blocked) >= 3 && length(v_both) >= 3) {
    wt <- tryCatch(wilcox.test(v_blocked, v_both, exact = FALSE),
                   error = function(e) list(p.value = NA))
    cat(sprintf("  %s: blocked=%.4f, concordant=%.4f, p=%.4f %s\n",
                feat, median(v_blocked), median(v_both), wt$p.value,
                ifelse(!is.na(wt$p.value) && wt$p.value < 0.05, "*", "")))
  }
}

# ============================================================================
# 9. SPECKLE INTRON SIGNATURE
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SPECKLE INTRON SIGNATURE ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

# According to Małszycki et al., speckle-associated introns in GC-rich isochores:
# - Higher GC content
# - Shorter introns (relative to gene size)
# - Flanked by short exons with weak splice sites
# For retained introns specifically: higher GC means harder to distinguish from exons

speckle_ir_signature <- data.frame(
  feature = c("INTRON_LENGTH", "INTRON_GCC",
              "UPSTREAM_EXON_GCC", "DOWNSTREAM_EXON_GCC",
              "PPT_PYRIMIDINE_FRAC", "SS5_GCC", "SS3_GCC"),
  expected_in_speckle = c("shorter", "higher",
                          "higher", "higher",
                          "lower", "higher", "higher"),
  stringsAsFactors = FALSE
)

cat("Speckle IR signature expectations:\n")
print(speckle_ir_signature, row.names = FALSE)

cat("\nComparing BLOCKED IR to NS:\n")
blocked_ir_vs_ns <- ir_quad_stats_df %>% filter(quadrant == "Blocked_by_Taxol")
for (i in seq_len(nrow(speckle_ir_signature))) {
  feat <- speckle_ir_signature$feature[i]
  row <- blocked_ir_vs_ns %>% filter(feature == feat)
  if (nrow(row) > 0) {
    direction <- ifelse(row$diff_medians > 0, "higher", ifelse(row$diff_medians < 0, "lower/shorter", "equal"))
    expected <- speckle_ir_signature$expected_in_speckle[i]
    matches <- (expected == "higher" && row$diff_medians > 0) ||
               (expected %in% c("lower", "shorter") && row$diff_medians < 0)
    cat(sprintf("  %s: %s (diff=%.4f, p=%.4f) — expected %s — %s\n",
                feat, direction, row$diff_medians, row$p.value,
                expected, ifelse(matches, "MATCHES", "does NOT match")))
  }
}

cat("\nComparing CONCORDANT IR to NS:\n")
conc_ir_vs_ns <- ir_quad_stats_df %>% filter(quadrant == "Both_concordant")
for (i in seq_len(nrow(speckle_ir_signature))) {
  feat <- speckle_ir_signature$feature[i]
  row <- conc_ir_vs_ns %>% filter(feature == feat)
  if (nrow(row) > 0) {
    direction <- ifelse(row$diff_medians > 0, "higher", ifelse(row$diff_medians < 0, "lower/shorter", "equal"))
    expected <- speckle_ir_signature$expected_in_speckle[i]
    matches <- (expected == "higher" && row$diff_medians > 0) ||
               (expected %in% c("lower", "shorter") && row$diff_medians < 0)
    cat(sprintf("  %s: %s (diff=%.4f, p=%.4f) — expected %s — %s\n",
                feat, direction, row$diff_medians, row$p.value,
                expected, ifelse(matches, "MATCHES", "does NOT match")))
  }
}

# ============================================================================
# 10. RETENTION UP vs DOWN WITHIN BLOCKED
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("BLOCKED IR: RETENTION_UP vs RETENTION_DOWN\n")
cat(strrep("=", 70), "\n\n")

ret_up   <- ir_with_feat %>% filter(quadrant == "Blocked_by_Taxol", dpir_dmso > DPSI_THRESH)
ret_down <- ir_with_feat %>% filter(quadrant == "Blocked_by_Taxol", dpir_dmso < -DPSI_THRESH)

cat(sprintf("Retention_up (blocked): %d events\n", nrow(ret_up)))
cat(sprintf("Retention_down (blocked): %d events\n", nrow(ret_down)))

if (nrow(ret_up) >= 3 && nrow(ret_down) >= 3) {
  for (feat in avail_feats) {
    v_up <- ret_up[[feat]][!is.na(ret_up[[feat]])]
    v_dn <- ret_down[[feat]][!is.na(ret_down[[feat]])]

    if (length(v_up) >= 3 && length(v_dn) >= 3) {
      wt <- tryCatch(wilcox.test(v_up, v_dn, exact = FALSE),
                     error = function(e) list(p.value = NA))
      if (!is.na(wt$p.value) && wt$p.value < 0.1) {
        cat(sprintf("  %s: up=%.4f, down=%.4f, p=%.4f\n",
                    feat, median(v_up), median(v_dn), wt$p.value))
      }
    }
  }
} else {
  cat("  Not enough events in sub-groups for comparison\n")
}

# ============================================================================
# 11. CANONICAL vs NON-CANONICAL SPLICE SITES PER QUADRANT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("CANONICAL SPLICE SITE USAGE PER IR QUADRANT\n")
cat(strrep("=", 70), "\n\n")

for (q in c("NS", "Both_concordant", "Blocked_by_Taxol", "Taxol_enabled")) {
  q_events <- ir_with_feat %>% filter(quadrant == q)
  q_feats <- q_events %>% left_join(ir_features %>% select(EVENT, IS_CANONICAL), by = "EVENT")
  n_canon <- sum(q_feats$IS_CANONICAL, na.rm = TRUE)
  n_total <- sum(!is.na(q_feats$IS_CANONICAL))
  cat(sprintf("  %s: %d / %d canonical (%.1f%%)\n",
              q, n_canon, n_total, 100 * n_canon / max(n_total, 1)))
}

# ============================================================================
# 12. EXPRESSION OF IR-BLOCKED GENES
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("EXPRESSION OF IR-BLOCKED GENES\n")
cat(strrep("=", 70), "\n\n")

ir_blocked_genes <- unique(ir_scatter$GENE[ir_scatter$quadrant == "Blocked_by_Taxol"])
cat(sprintf("IR-Blocked genes: %d\n", length(ir_blocked_genes)))

for (x in list(
  list("Myo Taxol DE", deseq_myo),
  list("Tube Taxol DE", deseq_tube),
  list("Diff DMSO DE", deseq_diff_d),
  list("Diff Taxol DE", deseq_diff_t)
)) {
  de <- x[[2]] %>%
    filter(!is.na(significance) & significance != "Not Significant") %>%
    pull(gene_symbol)
  ov <- sum(ir_blocked_genes %in% de)
  cat(sprintf("  %s: %d / %d ir-blocked genes DE (%.1f%%)\n",
              x[[1]], ov, length(ir_blocked_genes),
              100 * ov / max(length(ir_blocked_genes), 1)))
}

# ============================================================================
# 13. TAXOL SCATTER — IR EVENTS (Myoblast vs Myotube)
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TAXOL SCATTER — IR: MYOBLAST vs MYOTUBE\n")
cat(strrep("=", 70), "\n\n")

tm_ir <- myo_taxol_ir %>% select(EVENT, GENE, dpir_myo = deltapsi, FDR_myo = FDR)
tt_ir <- tube_taxol_ir %>% select(EVENT, dpir_tube = deltapsi, FDR_tube = FDR)

taxol_ir_scatter <- inner_join(tm_ir, tt_ir, by = "EVENT") %>%
  mutate(
    sig_myo  = is_sig(FDR_myo, dpir_myo),
    sig_tube = is_sig(FDR_tube, dpir_tube),
    quadrant = case_when(
      sig_myo & sig_tube & sign(dpir_myo) == sign(dpir_tube) ~ "Both_concordant",
      sig_myo & sig_tube & sign(dpir_myo) != sign(dpir_tube) ~ "Discordant",
      sig_myo & !sig_tube ~ "Myoblast_only",
      !sig_myo & sig_tube ~ "Myotube_only",
      TRUE ~ "NS"
    )
  )

cat("IR Taxol scatter quadrant counts:\n")
print(table(taxol_ir_scatter$quadrant))

# Speckle enrichment in taxol IR scatter
taxol_ir_scatter$GENE_upper <- toupper(taxol_ir_scatter$GENE)
for (q in c("Both_concordant", "Myoblast_only", "Myotube_only", "Discordant")) {
  q_genes <- unique(taxol_ir_scatter$GENE_upper[taxol_ir_scatter$quadrant == q])
  n_q <- length(q_genes)
  n_sp <- sum(q_genes %in% speckle_upper)
  cat(sprintf("  %s: %d genes, %d speckle (%.1f%%)\n",
              q, n_q, n_sp, 100 * n_sp / max(n_q, 1)))
}

# ============================================================================
# 14. COMPREHENSIVE SUMMARY TABLE
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("IR vs EX COMPREHENSIVE COMPARISON\n")
cat(strrep("=", 70), "\n\n")

# Side-by-side: EX quadrants vs IR quadrants
cat("Event Type | Quadrant          | Count | Speckle%\n")
cat(strrep("-", 55), "\n")

for (etype in c("EX", "IR")) {
  scatter <- if (etype == "EX") ex_scatter else ir_scatter
  scatter$GENE_upper <- toupper(scatter$GENE)
  for (q in c("Both_concordant", "Blocked_by_Taxol", "Taxol_enabled")) {
    n_events <- sum(scatter$quadrant == q)
    q_genes <- unique(scatter$GENE_upper[scatter$quadrant == q])
    n_sp <- sum(q_genes %in% speckle_upper)
    pct <- 100 * n_sp / max(length(q_genes), 1)
    cat(sprintf("%-9s | %-17s | %5d | %.1f%%\n", etype, q, n_events, pct))
  }
  cat(strrep("-", 55), "\n")
}

# Percentage blocked
ir_sig_dmso <- sum(is_sig(diff_dmso_ir$FDR, diff_dmso_ir$deltapsi))
ex_sig_dmso <- sum(is_sig(diff_dmso_ex$FDR, diff_dmso_ex$deltapsi))
ir_blocked_n <- sum(ir_scatter$quadrant == "Blocked_by_Taxol")
ex_blocked_n <- sum(ex_scatter$quadrant == "Blocked_by_Taxol")

cat(sprintf("\nEX: %d sig in DMSO diff, %d blocked by Taxol (%.1f%%)\n",
            ex_sig_dmso, ex_blocked_n, 100 * ex_blocked_n / max(ex_sig_dmso, 1)))
cat(sprintf("IR: %d sig in DMSO diff, %d blocked by Taxol (%.1f%%)\n",
            ir_sig_dmso, ir_blocked_n, 100 * ir_blocked_n / max(ir_sig_dmso, 1)))

cat("\nDone.\n")
