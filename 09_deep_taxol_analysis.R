#!/usr/bin/env Rscript
# ============================================================================
# 09_deep_taxol_analysis.R
# Deep analysis of taxol effects: myotube vs myoblast, IR vs EX,
# forces/speckle model. Text-only output (no plots).
# Publication figures based on these analyses are in 10_final_figures.R.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

setwd("/Users/andres.ortiz/Projects/CRG/taxol_rnaseq_c2c12")
RESULTS_DIR <- "results"

# ============================================================================
# 1. LOAD ALL DATA
# ============================================================================

cat("=== Loading data ===\n")
diff_dmso  <- readRDS(file.path(RESULTS_DIR, "Differentiation_DMSO.rds"))
diff_taxol <- readRDS(file.path(RESULTS_DIR, "Differentiation_Taxol.rds"))
myo_taxol  <- readRDS(file.path(RESULTS_DIR, "Myoblast_Taxol_vs_DMSO.rds"))
tube_taxol <- readRDS(file.path(RESULTS_DIR, "Myotube_Taxol_vs_DMSO.rds"))

speckle_genes <- scan("nuclear_speckle_associated_genes.txt", what = "character", quiet = TRUE)
speckle_upper <- toupper(speckle_genes)

ir_features <- read.csv("results/ir_features_precomputed.csv", stringsAsFactors = FALSE)

FDR_THRESH <- 0.05; DPSI_THRESH <- 0.1
is_sig <- function(fdr, dpsi) !is.na(fdr) & fdr <= FDR_THRESH & !is.na(dpsi) & abs(dpsi) >= DPSI_THRESH
filter_ir <- function(df) df %>% filter(grepl("MmuINT", EVENT))
filter_ex <- function(df) df %>% filter(grepl("MmuEX", EVENT))

# ============================================================================
# 2. TAXOL SCATTER: MYOBLAST vs MYOTUBE — ALL, EX, IR
# ============================================================================

cat("\n=== TAXOL SCATTER: MYOBLAST vs MYOTUBE ===\n\n")

# Build unified scatter with event types
build_taxol_scatter <- function(myo_df, tube_df, type_label) {
  m <- myo_df %>% select(EVENT, GENE, dpsi_myo = deltapsi, FDR_myo = FDR)
  t <- tube_df %>% select(EVENT, dpsi_tube = deltapsi, FDR_tube = FDR)
  inner_join(m, t, by = "EVENT") %>%
    mutate(
      sig_myo = is_sig(FDR_myo, dpsi_myo),
      sig_tube = is_sig(FDR_tube, dpsi_tube),
      quadrant = case_when(
        sig_myo & sig_tube & sign(dpsi_myo) == sign(dpsi_tube) ~ "Both",
        sig_myo & sig_tube & sign(dpsi_myo) != sign(dpsi_tube) ~ "Discordant",
        sig_myo & !sig_tube ~ "Myoblast_only",
        !sig_myo & sig_tube ~ "Myotube_only",
        TRUE ~ "NS"),
      event_type = type_label,
      is_speckle = toupper(GENE) %in% speckle_upper
    )
}

tx_all <- build_taxol_scatter(myo_taxol, tube_taxol, "All")
tx_ir  <- build_taxol_scatter(filter_ir(myo_taxol), filter_ir(tube_taxol), "IR")
tx_ex  <- build_taxol_scatter(filter_ex(myo_taxol), filter_ex(tube_taxol), "EX")

for (label in c("All", "EX", "IR")) {
  df <- switch(label, "All" = tx_all, "EX" = tx_ex, "IR" = tx_ir)
  cat(sprintf("--- %s events: Taxol scatter ---\n", label))
  cat("Quadrant counts:\n")
  print(table(df$quadrant))

  # Direction within each quadrant
  for (q in c("Myoblast_only", "Myotube_only", "Both", "Discordant")) {
    qd <- df %>% filter(quadrant == q)
    if (nrow(qd) == 0) next

    if (label == "IR") {
      n_ret_up_myo <- sum(qd$dpsi_myo > 0)
      n_ret_dn_myo <- sum(qd$dpsi_myo < 0)
      n_ret_up_tube <- sum(qd$dpsi_tube > 0)
      n_ret_dn_tube <- sum(qd$dpsi_tube < 0)
      cat(sprintf("  %s (%d): myo_ret↑=%d myo_ret↓=%d | tube_ret↑=%d tube_ret↓=%d\n",
                  q, nrow(qd), n_ret_up_myo, n_ret_dn_myo, n_ret_up_tube, n_ret_dn_tube))
    } else if (label == "EX") {
      n_inc_myo <- sum(qd$dpsi_myo > 0)
      n_skp_myo <- sum(qd$dpsi_myo < 0)
      n_inc_tube <- sum(qd$dpsi_tube > 0)
      n_skp_tube <- sum(qd$dpsi_tube < 0)
      cat(sprintf("  %s (%d): myo_inc=%d myo_skp=%d | tube_inc=%d tube_skp=%d\n",
                  q, nrow(qd), n_inc_myo, n_skp_myo, n_inc_tube, n_skp_tube))
    }

    # Speckle enrichment
    qg <- unique(toupper(qd$GENE))
    n_sp <- sum(qg %in% speckle_upper)
    cat(sprintf("    Speckle: %d/%d genes (%.1f%%)\n", n_sp, length(qg), 100*n_sp/max(length(qg),1)))
  }
  cat("\n")
}

# ============================================================================
# 3. KEY QUESTION: Myotube-specific taxol IR events — are these speckle/GC-rich?
# ============================================================================

cat("\n=== MYOTUBE-SPECIFIC TAXOL IR EVENTS ===\n")
cat("(These are the events where forces should matter most)\n\n")

tube_only_ir <- tx_ir %>% filter(quadrant == "Myotube_only")
cat(sprintf("Myotube-only IR events: %d\n", nrow(tube_only_ir)))
cat(sprintf("  Retention increased (introns STAY IN under taxol): %d\n",
            sum(tube_only_ir$dpsi_tube > 0)))
cat(sprintf("  Retention decreased (introns get EXCISED under taxol): %d\n",
            sum(tube_only_ir$dpsi_tube < 0)))

# Speckle enrichment of myotube-only IR events
tube_only_ir_genes <- unique(toupper(tube_only_ir$GENE))
bg_ir_genes <- unique(toupper(tx_ir$GENE))
n_bg <- length(bg_ir_genes)
n_sp_bg <- sum(bg_ir_genes %in% speckle_upper)

cat(sprintf("\nSpeckle enrichment:\n"))
cat(sprintf("  IR Background: %d genes, %d speckle (%.1f%%)\n",
            n_bg, n_sp_bg, 100*n_sp_bg/n_bg))

# By direction
tube_ret_up <- tube_only_ir %>% filter(dpsi_tube > 0)
tube_ret_dn <- tube_only_ir %>% filter(dpsi_tube < 0)

for (x in list(
  list("Myotube_only (all IR)", tube_only_ir),
  list("Myotube_only retention↑ (kept in)", tube_ret_up),
  list("Myotube_only retention↓ (excised)", tube_ret_dn)
)) {
  g <- unique(toupper(x[[2]]$GENE))
  n_g <- length(g); n_sp <- sum(g %in% speckle_upper)
  if (n_g >= 2 & n_sp_bg > 0) {
    cont <- matrix(c(n_sp, n_g - n_sp, n_sp_bg - n_sp, max(0, n_bg - n_g - n_sp_bg + n_sp)), nrow=2)
    ft <- tryCatch(fisher.test(cont), error = function(e) list(estimate=NA, p.value=1))
    cat(sprintf("  %s: %d genes, %d speckle (%.1f%%), OR=%.2f, p=%.3g\n",
                x[[1]], n_g, n_sp, 100*n_sp/max(n_g,1), unname(ft$estimate), ft$p.value))
  } else {
    cat(sprintf("  %s: %d genes, %d speckle (%.1f%%)\n",
                x[[1]], n_g, n_sp, 100*n_sp/max(n_g,1)))
  }
}

# Features of myotube-only retention↑ IR events
cat("\n--- Features of Myotube-only retention↑ IR events ---\n")
tube_ret_up_feat <- tube_ret_up %>% left_join(ir_features, by = "EVENT")
ns_ir_feat <- tx_ir %>% filter(quadrant == "NS") %>% left_join(ir_features, by = "EVENT")

for (feat in c("INTRON_LENGTH", "INTRON_GCC", "UPSTREAM_EXON_GCC", "DOWNSTREAM_EXON_GCC")) {
  v_q <- tube_ret_up_feat[[feat]][!is.na(tube_ret_up_feat[[feat]])]
  v_ns <- ns_ir_feat[[feat]][!is.na(ns_ir_feat[[feat]])]
  if (length(v_q) >= 3) {
    wt <- wilcox.test(v_q, v_ns, exact = FALSE)
    cat(sprintf("  %s: tube_ret↑=%.3f, NS=%.3f, p=%.4g\n",
                feat, median(v_q), median(v_ns), wt$p.value))
  }
}

# ============================================================================
# 4. MYOBLAST-SPECIFIC TAXOL IR EVENTS
# ============================================================================

cat("\n=== MYOBLAST-SPECIFIC TAXOL IR EVENTS ===\n\n")

myo_only_ir <- tx_ir %>% filter(quadrant == "Myoblast_only")
cat(sprintf("Myoblast-only IR events: %d\n", nrow(myo_only_ir)))
cat(sprintf("  Retention increased: %d\n", sum(myo_only_ir$dpsi_myo > 0)))
cat(sprintf("  Retention decreased: %d\n", sum(myo_only_ir$dpsi_myo < 0)))

myo_ret_up <- myo_only_ir %>% filter(dpsi_myo > 0)
myo_ret_dn <- myo_only_ir %>% filter(dpsi_myo < 0)

for (x in list(
  list("Myoblast_only (all IR)", myo_only_ir),
  list("Myoblast_only retention↑", myo_ret_up),
  list("Myoblast_only retention↓", myo_ret_dn)
)) {
  g <- unique(toupper(x[[2]]$GENE))
  n_g <- length(g); n_sp <- sum(g %in% speckle_upper)
  if (n_g >= 2 & n_sp_bg > 0) {
    cont <- matrix(c(n_sp, n_g - n_sp, n_sp_bg - n_sp, max(0, n_bg - n_g - n_sp_bg + n_sp)), nrow=2)
    ft <- tryCatch(fisher.test(cont), error = function(e) list(estimate=NA, p.value=1))
    cat(sprintf("  %s: %d genes, %d speckle (%.1f%%), OR=%.2f, p=%.3g\n",
                x[[1]], n_g, n_sp, 100*n_sp/max(n_g,1), unname(ft$estimate), ft$p.value))
  }
}

cat("\n--- Features of Myoblast-only retention↓ IR events ---\n")
myo_ret_dn_feat <- myo_ret_dn %>% left_join(ir_features, by = "EVENT")

for (feat in c("INTRON_LENGTH", "INTRON_GCC", "UPSTREAM_EXON_GCC", "DOWNSTREAM_EXON_GCC")) {
  v_q <- myo_ret_dn_feat[[feat]][!is.na(myo_ret_dn_feat[[feat]])]
  v_ns <- ns_ir_feat[[feat]][!is.na(ns_ir_feat[[feat]])]
  if (length(v_q) >= 3) {
    wt <- wilcox.test(v_q, v_ns, exact = FALSE)
    cat(sprintf("  %s: myo_ret↓=%.3f, NS=%.3f, p=%.4g\n",
                feat, median(v_q), median(v_ns), wt$p.value))
  }
}

# ============================================================================
# 5. DISCORDANT EVENTS (opposite direction myo vs tube)
# ============================================================================

cat("\n=== DISCORDANT TAXOL IR EVENTS ===\n\n")
disc_ir <- tx_ir %>% filter(quadrant == "Discordant")
cat(sprintf("Discordant IR events: %d\n", nrow(disc_ir)))
if (nrow(disc_ir) > 0) {
  cat("  Myo_ret↑ + Tube_ret↓:", sum(disc_ir$dpsi_myo > 0 & disc_ir$dpsi_tube < 0), "\n")
  cat("  Myo_ret↓ + Tube_ret↑:", sum(disc_ir$dpsi_myo < 0 & disc_ir$dpsi_tube > 0), "\n")
}

# ============================================================================
# 6. COMPREHENSIVE DIRECTION TABLE
# ============================================================================

cat("\n=== COMPREHENSIVE DIRECTION TABLE ===\n\n")

cat("Table: Taxol effect on intron retention by cell type\n")
cat("(Retention↑ = intron STAYS IN = less efficient splicing)\n")
cat("(Retention↓ = intron EXCISED = more efficient splicing)\n\n")

cat(sprintf("%-20s | %6s %6s | %6s %6s | %8s\n",
            "Comparison", "Ret↑", "Ret↓", "%Ret↑", "%Ret↓", "Speckle%"))
cat(strrep("-", 75), "\n")

for (x in list(
  list("Myo Taxol (all)", filter_ir(myo_taxol)),
  list("Tube Taxol (all)", filter_ir(tube_taxol)),
  list("Myo Taxol (sig)", filter_ir(myo_taxol) %>% filter(is_sig(FDR, deltapsi))),
  list("Tube Taxol (sig)", filter_ir(tube_taxol) %>% filter(is_sig(FDR, deltapsi)))
)) {
  df <- x[[2]]
  sig <- df %>% filter(is_sig(FDR, deltapsi))
  n_up <- sum(sig$deltapsi > 0)
  n_dn <- sum(sig$deltapsi < 0)
  n_tot <- n_up + n_dn
  g <- unique(toupper(sig$GENE))
  sp_pct <- 100 * sum(g %in% speckle_upper) / max(length(g), 1)
  cat(sprintf("%-20s | %6d %6d | %5.1f%% %5.1f%% | %7.1f%%\n",
              x[[1]], n_up, n_dn, 100*n_up/max(n_tot,1), 100*n_dn/max(n_tot,1), sp_pct))
}

# ============================================================================
# 7. DIFFERENTIATION QUADRANTS — now with correct terminology
# ============================================================================

cat("\n=== DIFFERENTIATION QUADRANTS — CORRECTED TERMINOLOGY ===\n\n")

# Build IR differentiation scatter
dd_ir <- filter_ir(diff_dmso) %>% select(EVENT, GENE, dpsi_dmso = deltapsi, FDR_dmso = FDR)
dt_ir <- filter_ir(diff_taxol) %>% select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)
ir_diff <- inner_join(dd_ir, dt_ir, by = "EVENT") %>%
  mutate(
    sig_dmso = is_sig(FDR_dmso, dpsi_dmso),
    sig_taxol = is_sig(FDR_taxol, dpsi_taxol),
    quadrant = case_when(
      sig_dmso & sig_taxol & sign(dpsi_dmso) == sign(dpsi_taxol) ~ "Concordant",
      sig_dmso & sig_taxol & sign(dpsi_dmso) != sign(dpsi_taxol) ~ "Reversed",
      sig_dmso & !sig_taxol ~ "Blocked",
      !sig_dmso & sig_taxol ~ "Enabled",
      TRUE ~ "NS"),
    is_speckle = toupper(GENE) %in% speckle_upper
  )

blocked_ir <- ir_diff %>% filter(quadrant == "Blocked")
cat("Blocked IR events (introns that should change during differentiation but don't under taxol):\n")
cat(sprintf("  Total: %d\n", nrow(blocked_ir)))
cat(sprintf("  Retention↓ blocked (introns that should be EXCISED during diff, but STAY IN): %d (%.0f%%)\n",
            sum(blocked_ir$dpsi_dmso < 0), 100*sum(blocked_ir$dpsi_dmso < 0)/nrow(blocked_ir)))
cat(sprintf("  Retention↑ blocked (introns that should STAY IN during diff, but get excised): %d (%.0f%%)\n",
            sum(blocked_ir$dpsi_dmso > 0), 100*sum(blocked_ir$dpsi_dmso > 0)/nrow(blocked_ir)))

cat("\nMeaning: During normal differentiation, 297 introns become MORE EFFICIENTLY SPLICED OUT.")
cat("\nTaxol PREVENTS this efficient excision → introns remain retained in myotubes.\n")
cat("This is the hallmark of impaired speckle function: the spliceosome cannot\n")
cat("excise GC-rich introns without the concentrated splicing factors from enlarged speckles.\n")

# ============================================================================
# 8. OVERLAP: blocked differentiation IR ∩ myotube taxol IR
# ============================================================================

cat("\n=== OVERLAP: BLOCKED DIFF IR ∩ MYOTUBE TAXOL IR ===\n\n")

blocked_diff_events <- blocked_ir$EVENT
tube_sig_ir <- filter_ir(tube_taxol) %>% filter(is_sig(FDR, deltapsi))
myo_sig_ir <- filter_ir(myo_taxol) %>% filter(is_sig(FDR, deltapsi))

ov_tube <- sum(blocked_diff_events %in% tube_sig_ir$EVENT)
ov_myo  <- sum(blocked_diff_events %in% myo_sig_ir$EVENT)

cat(sprintf("Blocked diff IR events also sig by Taxol in myotubes: %d/%d (%.1f%%)\n",
            ov_tube, length(blocked_diff_events), 100*ov_tube/length(blocked_diff_events)))
cat(sprintf("Blocked diff IR events also sig by Taxol in myoblasts: %d/%d (%.1f%%)\n",
            ov_myo, length(blocked_diff_events), 100*ov_myo/length(blocked_diff_events)))

# What about at the gene level?
blocked_genes <- unique(toupper(blocked_ir$GENE))
tube_sig_genes <- unique(toupper(tube_sig_ir$GENE))
myo_sig_genes <- unique(toupper(myo_sig_ir$GENE))

cat(sprintf("\nAt gene level:\n"))
cat(sprintf("  Blocked diff IR genes: %d\n", length(blocked_genes)))
cat(sprintf("  Myotube taxol IR genes: %d, overlap with blocked: %d (%.1f%%)\n",
            length(tube_sig_genes), sum(blocked_genes %in% tube_sig_genes),
            100*sum(blocked_genes %in% tube_sig_genes)/length(blocked_genes)))
cat(sprintf("  Myoblast taxol IR genes: %d, overlap with blocked: %d (%.1f%%)\n",
            length(myo_sig_genes), sum(blocked_genes %in% myo_sig_genes),
            100*sum(blocked_genes %in% myo_sig_genes)/length(blocked_genes)))

# ============================================================================
# 9. THE FORCE STORY
# ============================================================================

cat("\n=== THE FORCE NARRATIVE ===\n\n")

# In myotubes: forces are maximal, speckles should be large
# Taxol in myotubes → blocks forces → speckles shrink → speckle-dependent IR fails
# Prediction: Taxol retention↑ events in myotubes should be speckle-enriched and GC-rich

tube_taxol_ir_sig <- filter_ir(tube_taxol) %>% filter(is_sig(FDR, deltapsi))
tube_ret_up_all <- tube_taxol_ir_sig %>% filter(deltapsi > 0)  # retention increases = introns stay
tube_ret_dn_all <- tube_taxol_ir_sig %>% filter(deltapsi < 0)  # retention decreases = excised

cat("Taxol effect on IR in MYOTUBES (where forces are maximal):\n")
cat(sprintf("  Retention↑ (introns STAY IN = speckle failure): %d events\n", nrow(tube_ret_up_all)))
cat(sprintf("  Retention↓ (introns excised = enhanced splicing): %d events\n", nrow(tube_ret_dn_all)))

# Features
tube_up_feat <- tube_ret_up_all %>% left_join(ir_features, by = "EVENT")
tube_dn_feat <- tube_ret_dn_all %>% left_join(ir_features, by = "EVENT")

cat("\nFeatures of retention↑ (speckle failure) in myotubes:\n")
for (feat in c("INTRON_LENGTH", "INTRON_GCC", "UPSTREAM_EXON_GCC", "DOWNSTREAM_EXON_GCC")) {
  v_up <- tube_up_feat[[feat]][!is.na(tube_up_feat[[feat]])]
  v_ns <- ns_ir_feat[[feat]][!is.na(ns_ir_feat[[feat]])]
  if (length(v_up) >= 3) {
    wt <- wilcox.test(v_up, v_ns, exact = FALSE)
    cat(sprintf("  %s: ret↑=%.3f, NS=%.3f, diff=%.3f, p=%.4g\n",
                feat, median(v_up), median(v_ns), median(v_up) - median(v_ns), wt$p.value))
  }
}

cat("\nTaxol effect on IR in MYOBLASTS (where forces are minimal):\n")
myo_taxol_ir_sig <- filter_ir(myo_taxol) %>% filter(is_sig(FDR, deltapsi))
myo_ret_up_all <- myo_taxol_ir_sig %>% filter(deltapsi > 0)
myo_ret_dn_all <- myo_taxol_ir_sig %>% filter(deltapsi < 0)
cat(sprintf("  Retention↑ (introns stay): %d events\n", nrow(myo_ret_up_all)))
cat(sprintf("  Retention↓ (introns excised): %d events\n", nrow(myo_ret_dn_all)))

cat("\nFeatures of retention↓ (paradoxical enhanced splicing) in myoblasts:\n")
myo_dn_feat <- myo_ret_dn_all %>% left_join(ir_features, by = "EVENT")
for (feat in c("INTRON_LENGTH", "INTRON_GCC", "UPSTREAM_EXON_GCC", "DOWNSTREAM_EXON_GCC")) {
  v_dn <- myo_dn_feat[[feat]][!is.na(myo_dn_feat[[feat]])]
  v_ns <- ns_ir_feat[[feat]][!is.na(ns_ir_feat[[feat]])]
  if (length(v_dn) >= 3) {
    wt <- wilcox.test(v_dn, v_ns, exact = FALSE)
    cat(sprintf("  %s: ret↓=%.3f, NS=%.3f, diff=%.3f, p=%.4g\n",
                feat, median(v_dn), median(v_ns), median(v_dn) - median(v_ns), wt$p.value))
  }
}

# ============================================================================
# 10. COMPARISON: EX direction reversal between cell types
# ============================================================================

cat("\n=== EX DIRECTION REVERSAL: MYOBLAST vs MYOTUBE ===\n\n")

myo_ex_sig <- filter_ex(myo_taxol) %>% filter(is_sig(FDR, deltapsi))
tube_ex_sig <- filter_ex(tube_taxol) %>% filter(is_sig(FDR, deltapsi))

cat(sprintf("Myoblast taxol EX: %d sig (%d inc, %d skp)\n",
            nrow(myo_ex_sig), sum(myo_ex_sig$deltapsi > 0), sum(myo_ex_sig$deltapsi < 0)))
cat(sprintf("Myotube taxol EX: %d sig (%d inc, %d skp)\n",
            nrow(tube_ex_sig), sum(tube_ex_sig$deltapsi > 0), sum(tube_ex_sig$deltapsi < 0)))

cat("\nIn myoblasts: taxol → skipping (48) ≈ inclusion (46) — balanced\n")
cat("In myotubes: taxol → inclusion (113) >> skipping (37) — 3:1 inclusion\n")
cat("This reversal parallels the IR reversal.\n")

# The combined story
cat("\n=== COMBINED EX + IR TAXOL DIRECTION BY CELL TYPE ===\n\n")
cat(sprintf("%-20s | %10s | %10s | %10s | %10s\n",
            "", "EX inc", "EX skp", "IR ret↑", "IR ret↓"))
cat(strrep("-", 75), "\n")
cat(sprintf("%-20s | %10d | %10d | %10d | %10d\n",
            "Myoblast Taxol",
            sum(myo_ex_sig$deltapsi > 0), sum(myo_ex_sig$deltapsi < 0),
            nrow(myo_ret_up_all), nrow(myo_ret_dn_all)))
cat(sprintf("%-20s | %10d | %10d | %10d | %10d\n",
            "Myotube Taxol",
            sum(tube_ex_sig$deltapsi > 0), sum(tube_ex_sig$deltapsi < 0),
            nrow(tube_ret_up_all), nrow(tube_ret_dn_all)))

cat("\n=== KEY INTERPRETATION ===\n")
cat("In MYOTUBES (max forces, max speckles):\n")
cat("  Taxol → 49 introns become MORE retained (speckle failure)\n")
cat("  Taxol → 113 exons become MORE included (loss of exclusion capacity)\n")
cat("  = Taxol impairs the splicing EFFICIENCY that large speckles provide\n\n")
cat("In MYOBLASTS (min forces, smaller speckles):\n")
cat("  Taxol → 77 introns become LESS retained (paradoxical improved excision)\n")
cat("  Taxol → 46 exons included, 48 skipped (balanced)\n")
cat("  = Different mechanism — possibly microtubule-dependent transport effects\n")

cat("\nDone with deep analysis.\n")
