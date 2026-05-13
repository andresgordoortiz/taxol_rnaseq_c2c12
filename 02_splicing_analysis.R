# ============================================================================
# 02_splicing_analysis.R
# Taxol Effect on mRNA Splicing in C2C12 Myoblasts and Myotubes
# Comprehensive Splicing Analysis & Visualization
#
# PURPOSE: Load FDR results from 01_fdr_calculation.R and perform a thorough
#          splicing analysis including QC, visualization, enrichment, and
#          feature-level characterization.
#
# EXPERIMENTAL DESIGN:
#   12 samples — 4 conditions x 3 replicates:
#     Myoblast_DMSO, Myoblast_Taxol, Myotube_DMSO, Myotube_Taxol
#
# REQUIRES: Results from 01_fdr_calculation.R in results/
#
# Author: Andrés Gordo Ortiz
# ============================================================================

cat(strrep("=", 70), "\n")
cat("02_splicing_analysis.R — Taxol Splicing Analysis\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(betAS)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(patchwork)
  library(showtext)
  library(sysfonts)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(colorRamp2)
  library(dendextend)
  library(UpSetR)
  library(grid)
  library(enrichR)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(biomaRt)
  library(ggpubr)
  library(viridis)
  library(scales)
})

# ============================================================================
# FONT SETUP (Courier Prime — consistent with lab styling)
# ============================================================================

# Test if Cairo actually works (capabilities() can lie on macOS without XQuartz)
HAS_CAIRO <- tryCatch({
  tmp <- tempfile(fileext = ".pdf")
  cairo_pdf(tmp, width = 1, height = 1)
  dev.off()
  unlink(tmp)
  TRUE
}, warning = function(w) FALSE, error = function(e) FALSE)

if (HAS_CAIRO) {
  font_add_google(name = "Courier Prime", family = "Courier Prime")
  showtext_auto()
  showtext_opts(dpi = 300)
  FONT_FAMILY <- "Courier Prime"
  SAVE_DEVICE <- cairo_pdf
  cat("  Font backend: showtext + cairo_pdf (Courier Prime)\n")
} else {
  FONT_FAMILY <- "Courier"
  SAVE_DEVICE <- grDevices::pdf
  cat("  Font backend: base pdf (Courier) — install XQuartz for Courier Prime\n")
}

open_pdf <- function(file, ...) {
  if (HAS_CAIRO) cairo_pdf(file, ...) else pdf(file, ...)
}

# ============================================================================
# SHARED THEME (adapted from q3_sb50_comparison.R styling)
# ============================================================================

theme_taxol <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = FONT_FAMILY) %+replace%
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold",
                                       size = rel(1.2), color = "black"),
      plot.subtitle    = element_text(hjust = 0.5, size = rel(0.9),
                                       color = "grey40"),
      axis.title       = element_text(face = "bold", size = rel(1.0)),
      axis.text        = element_text(size = rel(0.9), color = "black"),
      legend.title     = element_text(face = "bold"),
      legend.text      = element_text(size = rel(0.85)),
      legend.position  = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(fill = NA,
                                       color = alpha("black", 0.3),
                                       linewidth = 0.5),
      strip.text       = element_text(face = "bold", size = rel(1.0))
    )
}

# ============================================================================
# COLOR PALETTES
# ============================================================================

# Condition colors
condition_colors <- c(
  "Myoblast_DMSO"  = "#4BA3C3",
  "Myoblast_Taxol" = "#D62839",
  "Myotube_DMSO"   = "#2166AC",
  "Myotube_Taxol"  = "#B2182B"
)

# Significance colors
sig_colors <- c(
  "Included"        = "#D62839",
  "Skipped"         = "#4BA3C3",
  "Not Significant" = "grey80"
)

# Intron retention colors
ir_colors <- c(
  "Retained"        = "#D62839",
  "Spliced Out"     = "#4BA3C3",
  "Not Significant" = "grey80"
)

# ============================================================================
# PATHS
# ============================================================================

RESULTS_DIR  <- file.path(getwd(), "results")
PLOTS_DIR    <- file.path(getwd(), "plots")
METADATA_FILE <- file.path(getwd(), "metadata", "metadata.csv")
INCLUSION_TABLE <- file.path(getwd(),
                             "INCLUSION_LEVELS_FULL-mm10-12.tab")

if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("--- Loading inclusion table ---\n")
data <- getDataset(pathTables = INCLUSION_TABLE, tool = "vast-tools")
all_events <- filterEvents(getEvents(data, tool = "vast-tools"), N = 10)
exons   <- filterEvents(all_events, types = c("C1", "C2", "C3", "S", "MIC"), N = 10)
introns <- filterEvents(all_events, types = c("IR"), N = 10)
alt_ss  <- filterEvents(all_events, types = c("ANN", "ALTD", "ALTA"), N = 10)

metadata <- read.csv(METADATA_FILE, stringsAsFactors = FALSE) %>% arrange(condition)

cat("--- Loading FDR results ---\n")
result_files <- list.files(RESULTS_DIR, pattern = "\\.csv$", full.names = TRUE)
result_files <- result_files[!grepl("summary|enrichr_|deseq2_|integrated_|matt_|event_class|force_gene|interaction_|ir_features|ir_quadrant|nuclear_speckle|quadrant_feature|sf_corr|taxol_responsive|vst_matrix", result_files)]

fdr_results <- list()
for (f in result_files) {
  name <- tools::file_path_sans_ext(basename(f))
  fdr_results[[name]] <- read.csv(f)
  cat(sprintf("  Loaded: %s (%d events)\n", name, nrow(fdr_results[[name]])))
}

# Build groupList for betAS functions
groups <- unique(metadata$condition)
groupList <- lapply(seq_along(groups), function(i) {
  list(
    name    = groups[i],
    samples = metadata$sample_id[metadata$condition == groups[i]],
    color   = condition_colors[groups[i]]
  )
})
names(groupList) <- groups

# ============================================================================
# ============================================================================
#  SECTION 1: QUALITY CONTROL
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 1: QUALITY CONTROL\n")
cat(strrep("=", 70), "\n\n")

# ----------------------------------------------------------------------------
# 1.1 PCA on PSI values
# ----------------------------------------------------------------------------

cat("--- 1.1 PCA on PSI values ---\n")

pca_mat <- na.omit(data[, c("EVENT", all_events$Samples)])
rownames(pca_mat) <- pca_mat$EVENT
pca_mat <- t(pca_mat[, -1])

# Map sample IDs to conditions
condition_pca <- metadata$condition[match(rownames(pca_mat), metadata$sample_id)]

pca_result <- prcomp(pca_mat)
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample    <- rownames(pca_data)
pca_data$Condition <- condition_pca

var_pc1 <- round(100 * summary(pca_result)$importance[2, 1], 1)
var_pc2 <- round(100 * summary(pca_result)$importance[2, 2], 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA of Splicing Inclusion Levels (PSI)",
    x     = sprintf("PC1 (%s%%)", var_pc1),
    y     = sprintf("PC2 (%s%%)", var_pc2),
    color = NULL
  ) +
  coord_fixed() +
  theme_taxol(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(PLOTS_DIR, "01_pca_psi.pdf"), p_pca, width = 7, height = 5, device = SAVE_DEVICE)
cat("  Saved: 01_pca_psi.pdf\n")

# ----------------------------------------------------------------------------
# 1.2 PSI Distribution (U-shape check for exons)
# ----------------------------------------------------------------------------

cat("--- 1.2 PSI distributions ---\n")

# Flatten PSI values for each condition
psi_long <- exons$PSI %>%
  dplyr::select(EVENT, all_of(metadata$sample_id)) %>%
  pivot_longer(-EVENT, names_to = "sample_id", values_to = "PSI") %>%
  left_join(metadata[, c("sample_id", "condition")], by = "sample_id") %>%
  filter(!is.na(PSI))

p_dist <- ggplot(psi_long, aes(x = PSI, fill = condition)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity", color = "grey30",
                 linewidth = 0.2) +
  facet_wrap(~condition, ncol = 2) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Exon PSI Distribution per Condition",
    subtitle = "U-shape indicates proper exon quantification",
    x = "PSI", y = "Event Count"
  ) +
  theme_taxol() +
  theme(legend.position = "none")

ggsave(file.path(PLOTS_DIR, "02_psi_distribution.pdf"), p_dist,
       width = 8, height = 6, device = SAVE_DEVICE)
cat("  Saved: 02_psi_distribution.pdf\n")

# ============================================================================
# ============================================================================
#  SECTION 2: VOLCANO PLOTS
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 2: VOLCANO PLOTS\n")
cat(strrep("=", 70), "\n\n")

# Function to build a volcano plot (consistent styling)
make_volcano <- function(df, title_str, is_intron = FALSE) {
  df <- na.omit(df)
  max_y <- 2

  if (is_intron) {
    df$Significant <- case_when(
      df$FDR <= 0.05 & df$deltapsi >=  0.1 ~ "Retained",
      df$FDR <= 0.05 & df$deltapsi <= -0.1 ~ "Spliced Out",
      TRUE ~ "Not Significant"
    )
    cols <- ir_colors
  } else {
    df$Significant <- case_when(
      df$FDR <= 0.05 & df$deltapsi >=  0.1 ~ "Included",
      df$FDR <= 0.05 & df$deltapsi <= -0.1 ~ "Skipped",
      TRUE ~ "Not Significant"
    )
    cols <- sig_colors
  }

  df$negLogFDR <- -log10(df$FDR)

  set.seed(123)
  df <- df %>%
    mutate(
      plot_y = ifelse(is.infinite(negLogFDR) | negLogFDR > max_y, max_y, negLogFDR),
      is_capped = (is.infinite(negLogFDR) | negLogFDR > max_y),
      plot_y_jitter = ifelse(is_capped, max_y - runif(n(), 0, 0.1), plot_y)
    )

  n_sig <- sum(df$Significant != "Not Significant")

  ggplot(df, aes(x = deltapsi)) +
    geom_point(data = filter(df, !is_capped),
               aes(y = plot_y, color = Significant, alpha = Significant, size = Significant)) +
    geom_point(data = filter(df, is_capped),
               aes(y = plot_y_jitter, color = Significant, alpha = Significant, size = Significant),
               stroke = 0.5, show.legend = FALSE) +
    scale_color_manual(values = cols) +
    scale_alpha_manual(values = c(setNames(rep(0.6, 2), names(cols)[1:2]),
                                  "Not Significant" = 0.2)) +
    scale_size_manual(values = c(setNames(rep(2.5, 2), names(cols)[1:2]),
                                 "Not Significant" = 1.5)) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    scale_y_continuous(limits = c(0, max_y), expand = c(0, 0),
                       labels = function(x) ifelse(x == max_y, paste0("\u2265", max_y), x)) +
    labs(
      title = title_str,
      subtitle = sprintf("%d significant events (FDR \u2264 0.05, |\u0394PSI| \u2265 0.1)", n_sig),
      x = "\u0394PSI",
      y = expression(-log[10](FDR))
    ) +
    theme_taxol(base_size = 11) +
    theme(
      legend.title   = element_blank(),
      legend.position = "top",
      axis.line      = element_line(color = "black"),
      panel.border   = element_blank()
    )
}

# Generate volcano plots for each comparison
comparison_labels <- c(
  "Myoblast_Taxol_vs_DMSO" = "Myoblast: Taxol vs DMSO",
  "Myotube_Taxol_vs_DMSO"  = "Myotube: Taxol vs DMSO",
  "Differentiation_DMSO"   = "Differentiation (DMSO)",
  "Differentiation_Taxol"  = "Differentiation (Taxol)"
)

for (comp_name in names(comparison_labels)) {
  if (comp_name %in% names(fdr_results)) {
    title <- comparison_labels[comp_name]
    p <- make_volcano(fdr_results[[comp_name]], title)
    fname <- sprintf("03_volcano_%s.pdf", comp_name)
    ggsave(file.path(PLOTS_DIR, fname), p, width = 5.5, height = 4, device = SAVE_DEVICE)
    cat(sprintf("  Saved: %s\n", fname))
  }
}

# ============================================================================
# ============================================================================
#  SECTION 3: DIFFERENTIAL SPLICING OVERVIEW
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 3: DIFFERENTIAL SPLICING OVERVIEW\n")
cat(strrep("=", 70), "\n\n")

# ----------------------------------------------------------------------------
# 3.1 Summary bar chart of DE events per comparison
# ----------------------------------------------------------------------------

extract_sig <- function(df, name) {
  sig <- df %>%
    filter(!is.na(FDR), FDR <= 0.05, abs(deltapsi) >= 0.1) %>%
    mutate(
      direction = ifelse(deltapsi > 0, "Included/Retained", "Skipped/Spliced Out"),
      event_type = case_when(
        grepl("^MmuEX",  EVENT) ~ "Cassette Exon (EX)",
        grepl("^MmuINT", EVENT) ~ "Intron Retention (IR)",
        TRUE                    ~ "Other"
      )
    )
  data.frame(
    comparison = name,
    direction  = sig$direction,
    event_type = sig$event_type,
    stringsAsFactors = FALSE
  )
}

sig_summary_list <- list()
for (name in names(fdr_results)) {
  sig_summary_list[[name]] <- extract_sig(fdr_results[[name]], name)
}
sig_summary <- bind_rows(sig_summary_list)

# Parse comparison name into parts
sig_summary <- sig_summary %>%
  mutate(
    comparison_label = comparison_labels[comparison]
  )

p_summary <- ggplot(sig_summary, aes(x = comparison_label, fill = direction)) +
  geom_bar(position = "stack", width = 0.7, color = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = c("Included/Retained" = "#D62839",
                                "Skipped/Spliced Out" = "#4BA3C3")) +
  facet_wrap(~ event_type, scales = "free_y") +
  labs(
    title = "Differential Splicing Events per Comparison",
    subtitle = "FDR <= 0.05, |dPSI| >= 0.1 -- separated by event type",
    x = NULL, y = "Number of Events", fill = NULL
  ) +
  theme_taxol() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = rel(0.8)),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = rel(1.1))
  )

ggsave(file.path(PLOTS_DIR, "04_de_events_summary.pdf"), p_summary,
       width = 12, height = 5, device = SAVE_DEVICE)
cat("  Saved: 04_de_events_summary.pdf\n")

# ----------------------------------------------------------------------------
# 3.2 Taxol responsiveness: myoblasts vs myotubes
#     (q3-style combined scatter + stacked bar figure)
# ----------------------------------------------------------------------------
#
# HOW TO READ:
#   Scatter X-axis: Taxol effect in myoblasts (\u0394PSI)
#   Scatter Y-axis: Taxol effect in myotubes (\u0394PSI)
#
#   Along X-axis → event responds to Taxol in myoblasts ONLY
#   Along Y-axis → event responds to Taxol in myotubes ONLY
#   On diagonal  → event responds similarly in BOTH cell types
#   Off-diagonal → discordant response (opposite directions)
# ----------------------------------------------------------------------------

cat("--- 3.2 Taxol responsiveness: myoblasts vs myotubes (q3 style) ---\n")

DPSI_THRESH <- 0.1
FDR_THRESH  <- 0.05
is_sig_spl <- function(fdr, dpsi) {
  !is.na(fdr) & fdr <= FDR_THRESH & !is.na(dpsi) & abs(dpsi) >= DPSI_THRESH
}

if (all(c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO") %in%
        names(fdr_results))) {

  taxol_myo <- fdr_results[["Myoblast_Taxol_vs_DMSO"]] %>%
    dplyr::select(EVENT, dpsi_myoblast = deltapsi, FDR_myoblast = FDR)

  taxol_tube <- fdr_results[["Myotube_Taxol_vs_DMSO"]] %>%
    dplyr::select(EVENT, dpsi_myotube = deltapsi, FDR_myotube = FDR)

  scatter_data <- inner_join(taxol_myo, taxol_tube, by = "EVENT") %>%
    mutate(
      sig_myoblast = is_sig_spl(FDR_myoblast, dpsi_myoblast),
      sig_myotube  = is_sig_spl(FDR_myotube,  dpsi_myotube),
      category = case_when(
        sig_myoblast & sig_myotube &
          sign(dpsi_myoblast) == sign(dpsi_myotube)  ~ "Both (concordant)",
        sig_myoblast & sig_myotube &
          sign(dpsi_myoblast) != sign(dpsi_myotube)  ~ "Discordant",
        sig_myoblast & !sig_myotube                  ~ "Myoblast only",
        !sig_myoblast & sig_myotube                  ~ "Myotube only",
        TRUE                                         ~ "NS"
      )
    )

  # Category colours (q3-inspired palette)
  cat_colors_scatter <- c(
    "Both (concordant)" = "#762A83",
    "Discordant"        = "#1B7837",
    "Myoblast only"     = "#D62839",
    "Myotube only"      = "#2166AC",
    "NS"                = "grey85"
  )

  # ---- SCATTER PLOT ----
  n_both  <- sum(scatter_data$category == "Both (concordant)")
  n_disc  <- sum(scatter_data$category == "Discordant")
  n_myo   <- sum(scatter_data$category == "Myoblast only")
  n_tube  <- sum(scatter_data$category == "Myotube only")

  p_taxol_scatter <- ggplot(scatter_data,
                            aes(x = dpsi_myoblast, y = dpsi_myotube,
                                color = category)) +
    geom_point(data = filter(scatter_data, category == "NS"),
               alpha = 0.1, size = 0.5) +
    geom_point(data = filter(scatter_data, category != "NS"),
               alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = c(-DPSI_THRESH, DPSI_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-DPSI_THRESH, DPSI_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red",
                linetype = "dotted", linewidth = 0.6) +
    scale_color_manual(
      values = cat_colors_scatter,
      breaks = c("Both (concordant)", "Discordant", "Myoblast only", "Myotube only")
    ) +
    labs(
      title    = "Taxol Splicing Effect",
      subtitle = sprintf("Both: %d | Myoblast: %d | Myotube: %d | Discordant: %d",
                          n_both, n_myo, n_tube, n_disc),
      x = expression(Delta*PSI~"(Myoblasts: Taxol vs DMSO)"),
      y = expression(Delta*PSI~"(Myotubes: Taxol vs DMSO)"),
      color = NULL
    ) +
    coord_fixed(ratio = 1, xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_taxol(base_size = 10) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50"),
      legend.position = "none"
    )

  ggsave(file.path(PLOTS_DIR, "05_taxol_scatter.pdf"), p_taxol_scatter,
         width = 6, height = 6, device = SAVE_DEVICE)
  cat("  Saved: 05_taxol_scatter.pdf\n")

  # ---- 100% STACKED BARPLOT ----
  scatter_prop <- scatter_data %>%
    filter(category != "NS") %>%
    mutate(category = factor(category,
                              levels = c("Myotube only", "Discordant",
                                         "Both (concordant)", "Myoblast only"))) %>%
    dplyr::count(category, .drop = FALSE) %>%
    mutate(total = sum(n), pct = 100 * n / total)

  p_taxol_bar <- ggplot(scatter_prop, aes(x = "Taxol", y = pct, fill = category)) +
    geom_col(position = "stack", width = 0.55, color = "grey30", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(pct > 5, sprintf("%d\n(%.0f%%)", n, pct), "")),
      position = position_stack(vjust = 0.5),
      size = 2.8, color = "white", fontface = "bold", family = FONT_FAMILY
    ) +
    scale_fill_manual(
      values = cat_colors_scatter,
      breaks = c("Both (concordant)", "Discordant", "Myoblast only", "Myotube only")
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(
      title    = "Category proportions",
      x = NULL, y = "% of significant events", fill = NULL
    ) +
    theme_taxol(base_size = 9) +
    theme(
      plot.title     = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      legend.text     = element_text(size = 7),
      panel.grid.major.x = element_blank(),
      axis.text.x     = element_text(size = 8),
      aspect.ratio    = 1
    ) +
    guides(fill = guide_legend(nrow = 2))

  ggsave(file.path(PLOTS_DIR, "06_taxol_category_bar.pdf"), p_taxol_bar,
         width = 5, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 06_taxol_category_bar.pdf\n")
}

# ----------------------------------------------------------------------------
# 3.3 Differentiation responsiveness: DMSO vs Taxol context
#     Does Taxol change the differentiation splicing programme?
# ----------------------------------------------------------------------------
#
# HOW TO READ:
#   X-axis: differentiation effect WITHOUT Taxol (\u0394PSI DMSO)
#   Y-axis: differentiation effect WITH Taxol (\u0394PSI Taxol)
#
#   On diagonal → Taxol does NOT alter differentiation-associated splicing
#   Along X-axis → differentiation event is BLOCKED by Taxol
#   Along Y-axis → differentiation event REQUIRES Taxol
# ----------------------------------------------------------------------------

cat("--- 3.3 Differentiation responsiveness: DMSO vs Taxol context ---\n")

if (all(c("Differentiation_DMSO", "Differentiation_Taxol") %in%
        names(fdr_results))) {

  diff_dmso <- fdr_results[["Differentiation_DMSO"]] %>%
    dplyr::select(EVENT, dpsi_dmso = deltapsi, FDR_dmso = FDR)

  diff_taxol <- fdr_results[["Differentiation_Taxol"]] %>%
    dplyr::select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)

  diff_scatter <- inner_join(diff_dmso, diff_taxol, by = "EVENT") %>%
    mutate(
      sig_dmso  = is_sig_spl(FDR_dmso,  dpsi_dmso),
      sig_taxol = is_sig_spl(FDR_taxol, dpsi_taxol),
      category = case_when(
        sig_dmso & sig_taxol &
          sign(dpsi_dmso) == sign(dpsi_taxol) ~ "Both (concordant)",
        sig_dmso & sig_taxol &
          sign(dpsi_dmso) != sign(dpsi_taxol) ~ "Reversed",
        sig_dmso & !sig_taxol                 ~ "Blocked by Taxol",
        !sig_dmso & sig_taxol                 ~ "Taxol-enabled",
        TRUE                                  ~ "NS"
      )
    )

  # Differentiation colours
  diff_cat_colors <- c(
    "Both (concordant)" = "#762A83",
    "Reversed"          = "#1B7837",
    "Blocked by Taxol"  = "#2166AC",
    "Taxol-enabled"     = "#D62839",
    "NS"                = "grey85"
  )

  n_both_d  <- sum(diff_scatter$category == "Both (concordant)")
  n_rev_d   <- sum(diff_scatter$category == "Reversed")
  n_block_d <- sum(diff_scatter$category == "Blocked by Taxol")
  n_enab_d  <- sum(diff_scatter$category == "Taxol-enabled")

  p_diff_scatter <- ggplot(diff_scatter,
                            aes(x = dpsi_dmso, y = dpsi_taxol,
                                color = category)) +
    geom_point(data = filter(diff_scatter, category == "NS"),
               alpha = 0.1, size = 0.5) +
    geom_point(data = filter(diff_scatter, category != "NS"),
               alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = c(-DPSI_THRESH, DPSI_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-DPSI_THRESH, DPSI_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red",
                linetype = "dotted", linewidth = 0.6) +
    scale_color_manual(
      values = diff_cat_colors,
      breaks = c("Both (concordant)", "Reversed", "Blocked by Taxol", "Taxol-enabled")
    ) +
    labs(
      title    = "Differentiation Splicing Programme",
      subtitle = sprintf("Conserved: %d | Blocked: %d | Enabled: %d | Reversed: %d",
                          n_both_d, n_block_d, n_enab_d, n_rev_d),
      x = expression(Delta*PSI~"(Differentiation in DMSO)"),
      y = expression(Delta*PSI~"(Differentiation in Taxol)"),
      color = NULL
    ) +
    coord_fixed(ratio = 1, xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_taxol(base_size = 10) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50"),
      legend.position = "none"
    )

  ggsave(file.path(PLOTS_DIR, "05b_diff_scatter.pdf"), p_diff_scatter,
         width = 6, height = 6, device = SAVE_DEVICE)
  cat("  Saved: 05b_diff_scatter.pdf\n")

  # ---- STACKED BAR ----
  diff_prop <- diff_scatter %>%
    filter(category != "NS") %>%
    mutate(category = factor(category,
                              levels = c("Taxol-enabled", "Reversed",
                                         "Both (concordant)", "Blocked by Taxol"))) %>%
    dplyr::count(category, .drop = FALSE) %>%
    mutate(total = sum(n), pct = 100 * n / total)

  p_diff_bar <- ggplot(diff_prop, aes(x = "Differentiation", y = pct, fill = category)) +
    geom_col(position = "stack", width = 0.55, color = "grey30", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(pct > 5, sprintf("%d\n(%.0f%%)", n, pct), "")),
      position = position_stack(vjust = 0.5),
      size = 2.8, color = "white", fontface = "bold", family = FONT_FAMILY
    ) +
    scale_fill_manual(
      values = diff_cat_colors,
      breaks = c("Both (concordant)", "Reversed", "Blocked by Taxol", "Taxol-enabled")
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(
      title    = "Category proportions",
      x = NULL, y = "% of significant events", fill = NULL
    ) +
    theme_taxol(base_size = 9) +
    theme(
      plot.title     = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      legend.text     = element_text(size = 7),
      panel.grid.major.x = element_blank(),
      axis.text.x     = element_text(size = 8),
      aspect.ratio    = 1
    ) +
    guides(fill = guide_legend(nrow = 2))

  ggsave(file.path(PLOTS_DIR, "06b_diff_category_bar.pdf"), p_diff_bar,
         width = 5, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 06b_diff_category_bar.pdf\n")
}

# ----------------------------------------------------------------------------
# 3.4 Combined 2×2 figure (q3-style)
# ----------------------------------------------------------------------------

cat("--- 3.4 Combined taxol × differentiation figure ---\n")

if (exists("p_taxol_scatter") && exists("p_diff_scatter") &&
    exists("p_taxol_bar")     && exists("p_diff_bar")) {

  combined_q3 <- (p_taxol_scatter + p_diff_scatter) /
                 (p_taxol_bar     + p_diff_bar) +
    plot_layout(widths = c(1, 1), heights = c(1, 1)) +
    plot_annotation(
      title = "Taxol × Differentiation Splicing Analysis",
      subtitle = paste0("Left: Cell-type specificity of Taxol effect  |  ",
                        "Right: Taxol impact on differentiation programme"),
      tag_levels = "A",
      theme = theme(
        plot.title    = element_text(size = 16, face = "bold",
                                      hjust = 0.5, family = FONT_FAMILY),
        plot.subtitle = element_text(size = 9, hjust = 0.5,
                                      color = "grey40", family = FONT_FAMILY)
      )
    )

  ggsave(file.path(PLOTS_DIR, "05c_combined_taxol_diff.pdf"), combined_q3,
         width = 10, height = 10, device = SAVE_DEVICE)
  cat("  Saved: 05c_combined_taxol_diff.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 4: UPSET PLOTS — EVENT & GENE LEVEL OVERLAP
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 4: UPSET PLOTS\n")
cat(strrep("=", 70), "\n\n")

get_sig_events <- function(name) {
  if (name %in% names(fdr_results)) {
    df <- fdr_results[[name]]
    na.omit(df[df$FDR <= 0.05 & abs(df$deltapsi) >= 0.1, ])$EVENT
  } else character(0)
}

get_sig_genes <- function(name) {
  if (name %in% names(fdr_results)) {
    df <- fdr_results[[name]]
    unique(na.omit(df[df$FDR <= 0.05 & abs(df$deltapsi) >= 0.1, ])$GENE)
  } else character(0)
}

# Event-level UpSet
event_list <- list(
  "Myoblast: Taxol vs DMSO"   = get_sig_events("Myoblast_Taxol_vs_DMSO"),
  "Myotube: Taxol vs DMSO"    = get_sig_events("Myotube_Taxol_vs_DMSO"),
  "Differentiation (DMSO)"    = get_sig_events("Differentiation_DMSO"),
  "Differentiation (Taxol)"   = get_sig_events("Differentiation_Taxol")
)

# Only plot if we have data
nonempty <- sapply(event_list, length) > 0
if (sum(nonempty) >= 2) {
  open_pdf(file.path(PLOTS_DIR, "07_upset_events.pdf"), width = 8, height = 5)
  binary_mat <- fromList(event_list[nonempty])
  UpSetR::upset(
    binary_mat,
    sets = names(event_list[nonempty]),
    sets.bar.color = "#2166AC",
    order.by = "freq",
    mainbar.y.label = "Shared Splicing Events",
    sets.x.label = "Total Significant Events",
    keep.order = TRUE,
    text.scale = c(1.3, 1.2, 1.0, 1.0, 1.3, 1.2),
    mb.ratio = c(0.65, 0.35),
    main.bar.color = "#4BA3C3",
    point.size = 3, line.size = 0.8
  )
  grid.text("Intersection of Significant Splicing Events",
            x = 0.5, y = 0.97,
            gp = gpar(fontsize = 12, fontface = "bold", fontfamily = FONT_FAMILY))
  dev.off()
  cat("  Saved: 07_upset_events.pdf\n")
}

# Gene-level UpSet
gene_list <- list(
  "Myoblast: Taxol vs DMSO"   = get_sig_genes("Myoblast_Taxol_vs_DMSO"),
  "Myotube: Taxol vs DMSO"    = get_sig_genes("Myotube_Taxol_vs_DMSO"),
  "Differentiation (DMSO)"    = get_sig_genes("Differentiation_DMSO"),
  "Differentiation (Taxol)"   = get_sig_genes("Differentiation_Taxol")
)

nonempty_g <- sapply(gene_list, length) > 0
if (sum(nonempty_g) >= 2) {
  open_pdf(file.path(PLOTS_DIR, "08_upset_genes.pdf"), width = 8, height = 5)
  binary_mat_g <- fromList(gene_list[nonempty_g])
  UpSetR::upset(
    binary_mat_g,
    sets = names(gene_list[nonempty_g]),
    sets.bar.color = "#B2182B",
    order.by = "freq",
    mainbar.y.label = "Shared Genes",
    sets.x.label = "Total Significant Genes",
    keep.order = TRUE,
    text.scale = c(1.3, 1.2, 1.0, 1.0, 1.3, 1.2),
    mb.ratio = c(0.65, 0.35),
    main.bar.color = "#D62839",
    point.size = 3, line.size = 0.8
  )
  grid.text("Intersection of Significantly Spliced Genes",
            x = 0.5, y = 0.97,
            gp = gpar(fontsize = 12, fontface = "bold", fontfamily = FONT_FAMILY))
  dev.off()
  cat("  Saved: 08_upset_genes.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 5: PSI DISTRIBUTION & SKEWNESS BY GENE COMPLEXITY
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 5: EXON NUMBER vs SPLICING SKEWNESS\n")
cat(strrep("=", 70), "\n\n")

# Pick one primary comparison for this analysis (Taxol effect in myoblasts)
primary_comp <- "Myoblast_Taxol_vs_DMSO"

if (primary_comp %in% names(fdr_results)) {
  df_primary <- fdr_results[[primary_comp]] %>% filter(!is.na(deltapsi))

  # 5.1 dPSI histograms — all 4 comparisons
  cat("--- 5.1 dPSI distribution (all comparisons) ---\n")

  # Single primary comparison histogram
  pct_below <- round(mean(df_primary$deltapsi < 0) * 100, 1)
  pct_above <- round(mean(df_primary$deltapsi > 0) * 100, 1)

  p_dpsi <- ggplot(df_primary, aes(x = deltapsi)) +
    geom_histogram(aes(y = after_stat(density), fill = deltapsi > 0),
                   bins = 80, alpha = 0.85, color = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = c("#4BA3C3", "#D62839"),
                      labels = c("\u0394PSI < 0", "\u0394PSI > 0")) +
    annotate("label", x = -0.7, y = Inf, vjust = 1.5,
             label = sprintf("\u0394PSI < 0: %.1f%%", pct_below),
             fill = alpha("white", 0.8), color = "#4BA3C3",
             size = 4, fontface = "bold", label.size = 0.4,
             family = FONT_FAMILY) +
    annotate("label", x = 0.7, y = Inf, vjust = 1.5,
             label = sprintf("\u0394PSI > 0: %.1f%%", pct_above),
             fill = alpha("white", 0.8), color = "#D62839",
             size = 4, fontface = "bold", label.size = 0.4,
             family = FONT_FAMILY) +
    labs(
      title = "Exon \u0394PSI Distribution: Taxol vs DMSO (Myoblasts)",
      x = "\u0394PSI", y = "Density", fill = NULL
    ) +
    theme_taxol() +
    theme(legend.position = "right")

  ggsave(file.path(PLOTS_DIR, "10_dpsi_distribution.pdf"), p_dpsi,
         width = 8, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 10_dpsi_distribution.pdf\n")

  # Faceted histogram for all 4 comparisons
  all_dpsi_list <- list()
  for (cn in names(comparison_labels)) {
    if (cn %in% names(fdr_results)) {
      df_cn <- fdr_results[[cn]] %>%
        filter(!is.na(deltapsi)) %>%
        mutate(comparison = comparison_labels[cn])
      all_dpsi_list[[cn]] <- df_cn
    }
  }

  if (length(all_dpsi_list) > 0) {
    all_dpsi <- bind_rows(all_dpsi_list) %>%
      mutate(comparison = factor(comparison, levels = unname(comparison_labels)))

    # Compute per-facet percentages
    pct_annot <- all_dpsi %>%
      group_by(comparison) %>%
      summarise(
        pct_neg = round(mean(deltapsi < 0) * 100, 1),
        pct_pos = round(mean(deltapsi > 0) * 100, 1),
        .groups = "drop"
      )

    p_dpsi_all <- ggplot(all_dpsi, aes(x = deltapsi)) +
      geom_histogram(aes(y = after_stat(density), fill = deltapsi > 0),
                     bins = 80, alpha = 0.85, color = "grey30", linewidth = 0.2) +
      geom_label(data = pct_annot, aes(x = -0.7, y = Inf,
                                        label = sprintf("\u0394PSI < 0: %.1f%%", pct_neg)),
                 vjust = 1.5, fill = alpha("white", 0.8), color = "#4BA3C3",
                 size = 3, fontface = "bold", label.size = 0.4,
                 family = FONT_FAMILY, inherit.aes = FALSE) +
      geom_label(data = pct_annot, aes(x = 0.7, y = Inf,
                                        label = sprintf("\u0394PSI > 0: %.1f%%", pct_pos)),
                 vjust = 1.5, fill = alpha("white", 0.8), color = "#D62839",
                 size = 3, fontface = "bold", label.size = 0.4,
                 family = FONT_FAMILY, inherit.aes = FALSE) +
      scale_fill_manual(values = c("#4BA3C3", "#D62839"),
                        labels = c("\u0394PSI < 0", "\u0394PSI > 0")) +
      facet_wrap(~comparison, ncol = 2, scales = "free_y") +
      labs(
        title = "\u0394PSI Distributions Across All Comparisons",
        x = "\u0394PSI", y = "Density", fill = NULL
      ) +
      theme_taxol() +
      theme(legend.position = "right",
            strip.text = element_text(size = 10, face = "bold"))

    ggsave(file.path(PLOTS_DIR, "10b_dpsi_distribution_all.pdf"), p_dpsi_all,
           width = 10, height = 8, device = SAVE_DEVICE)
    cat("  Saved: 10b_dpsi_distribution_all.pdf\n")
  }

  # 5.2 Skewness by exon count per gene (bootstrapped)
  cat("--- 5.2 Bootstrapped skewness by gene complexity ---\n")

  # Use coding_unique_exon_counts_per_gene.txt (Ensembl ID → exon count)
  exon_count_file <- file.path(getwd(), "coding_unique_exon_counts_per_gene.txt")
  if (file.exists(exon_count_file)) {
    exon_counts_raw <- read.table(exon_count_file, header = FALSE, sep = "\t",
                                  col.names = c("ensembl_id", "exon_count"),
                                  stringsAsFactors = FALSE)
    cat(sprintf("  Loaded exon counts for %d Ensembl genes\n", nrow(exon_counts_raw)))

    # Map Ensembl IDs to gene symbols using biomaRt
    cat("  Mapping Ensembl IDs to gene symbols via biomaRt...\n")
    ensembl <- tryCatch({
      useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    }, error = function(e) {
      cat("  biomaRt primary mirror failed, trying uswest...\n")
      tryCatch(
        useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                host = "https://uswest.ensembl.org"),
        error = function(e2) NULL
      )
    })

    if (!is.null(ensembl)) {
      id_map <- tryCatch({
        getBM(
          attributes = c("ensembl_gene_id", "mgi_symbol"),
          filters    = "ensembl_gene_id",
          values     = exon_counts_raw$ensembl_id,
          mart       = ensembl
        )
      }, error = function(e) {
        cat("  getBM failed:", conditionMessage(e), "\n")
        NULL
      })
    } else {
      id_map <- NULL
    }

    if (!is.null(id_map) && nrow(id_map) > 0) {
      exon_counts_mapped <- merge(exon_counts_raw, id_map,
                                  by.x = "ensembl_id", by.y = "ensembl_gene_id") %>%
        filter(mgi_symbol != "") %>%
        group_by(mgi_symbol) %>%
        summarise(exon_count = max(exon_count), .groups = "drop")  # keep max if duplicates
      gene_exon_counts <- setNames(exon_counts_mapped$exon_count,
                                   exon_counts_mapped$mgi_symbol)
      cat(sprintf("  Mapped %d genes with exon counts\n", length(gene_exon_counts)))
    } else {
      # Fallback: use AnnotationDbi
      cat("  biomaRt unavailable, falling back to org.Mm.eg.db...\n")
      ens_ids <- exon_counts_raw$ensembl_id
      symbol_map <- tryCatch({
        AnnotationDbi::select(org.Mm.eg.db,
                              keys = ens_ids,
                              keytype = "ENSEMBL",
                              columns = "SYMBOL")
      }, error = function(e) data.frame(ENSEMBL = character(), SYMBOL = character()))
      exon_counts_mapped <- merge(exon_counts_raw, symbol_map,
                                  by.x = "ensembl_id", by.y = "ENSEMBL") %>%
        filter(!is.na(SYMBOL), SYMBOL != "") %>%
        group_by(SYMBOL) %>%
        summarise(exon_count = max(exon_count), .groups = "drop")
      gene_exon_counts <- setNames(exon_counts_mapped$exon_count,
                                   exon_counts_mapped$SYMBOL)
      cat(sprintf("  Mapped %d genes with exon counts (AnnotationDbi)\n",
                  length(gene_exon_counts)))
    }

    df_with_exons <- df_primary %>%
      mutate(numb_exons = as.numeric(gene_exon_counts[GENE])) %>%
      filter(!is.na(numb_exons), abs(deltapsi) > 0)

    skewness_manual <- function(x) {
      x <- na.omit(x)
      n <- length(x)
      if (n < 3) return(NA_real_)
      m <- mean(x); s <- sd(x)
      if (s == 0) return(NA_real_)
      sum((x - m)^3) / n / (s^3)
    }

    binned <- df_with_exons %>%
      filter(abs(deltapsi) > 0) %>%
      mutate(
        exon_bin = cut(numb_exons,
                       breaks = c(0, 10, 50, Inf),
                       labels = c("1\u201310", "11\u201350", "51+"),
                       right = TRUE)
      ) %>%
      filter(!is.na(exon_bin))

    # Run bootstrap with deterministic per-bin seeding
    bin_levels <- c("1\u201310", "11\u201350", "51+")
    boot_results_list <- list()
    for (i in seq_along(bin_levels)) {
      bl <- bin_levels[i]
      bin_data <- binned %>% filter(exon_bin == bl)
      if (nrow(bin_data) < 3) next
      set.seed(1000 + i)  # deterministic per bin
      bs <- replicate(1000, skewness_manual(sample(bin_data$deltapsi, replace = TRUE)))
      bs <- na.omit(bs)
      p_val <- if (length(bs) == 0) 1 else {
        obs <- mean(bs)
        2 * min(mean(bs <= 0), mean(bs >= 0))  # two-sided bootstrap p-value
      }
      boot_results_list[[bl]] <- tibble(
        exon_bin      = bl,
        n_events      = nrow(bin_data),
        skew_estimate = mean(bs),
        ci_lower      = quantile(bs, 0.025),
        ci_upper      = quantile(bs, 0.975),
        p_boot        = p_val,
        boot_skews    = list(bs)
      )
    }

    # Helper: tiered significance stars based on bootstrap p-value
    sig_stars <- function(p) {
      ifelse(p < 0.001, "***",
      ifelse(p < 0.01,  "**",
      ifelse(p < 0.05,  "*", "")))
    }

    boot_results <- bind_rows(boot_results_list) %>%
      mutate(
        exon_bin = factor(exon_bin, levels = bin_levels),
        stars    = sig_stars(p_boot)
      )

    # Unnest all 1000 bootstrap replicates for display
    boot_long <- boot_results %>%
      unnest(boot_skews)

    # Auto y-limits from actual data range
    y_pad <- 0.3
    y_lim <- c(
      min(boot_long$boot_skews, na.rm = TRUE) - y_pad,
      max(boot_long$boot_skews, na.rm = TRUE) + y_pad
    )

    p_skew <- ggplot() +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4, linetype = "dashed") +
      # All 1000 bootstrap dots — tiny and translucent for density cloud
      geom_jitter(data = boot_long,
                  aes(x = exon_bin, y = boot_skews, fill = exon_bin),
                  width = 0.22, shape = 21, color = NA,
                  size = 0.4, alpha = 0.04) +
      # 95% CI error bar
      geom_errorbar(data = boot_results,
                    aes(x = exon_bin, ymin = ci_lower, ymax = ci_upper),
                    width = 0.12, linewidth = 0.6, color = "grey20") +
      # Mean point
      geom_point(data = boot_results,
                 aes(x = exon_bin, y = skew_estimate, fill = exon_bin),
                 shape = 21, size = 3.5, stroke = 0.7, color = "black") +
      # Significance stars (tiered)
      geom_text(data = boot_results %>% filter(stars != ""),
                aes(x = exon_bin, y = ci_upper + 0.12, label = stars),
                size = 5, color = "#D62839", fontface = "bold",
                family = FONT_FAMILY) +
      scale_fill_manual(values = c(
        "1\u201310"  = "#d1cbe5",
        "11\u201350" = "#9671bd",
        "51+"        = "#6a408d"
      )) +
      coord_cartesian(ylim = y_lim) +
      labs(
        title = "Bootstrapped Skewness of Taxol Splicing",
        subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
        x = "Number of coding exons per gene",
        y = "Skewness of \u0394PSI"
      ) +
      theme_taxol(base_size = 11) +
      theme(aspect.ratio = 1, legend.position = "none")

    ggsave(file.path(PLOTS_DIR, "11_skewness_by_exon_count.pdf"), p_skew, device = SAVE_DEVICE,
           width = 5, height = 5)
    cat("  Saved: 11_skewness_by_exon_count.pdf\n")

    # 5.2b Stacked panel: dPSI violin + bootstrap skewness (q3 style)
    cat("--- 5.2b \u0394PSI violin + bootstrap skewness (stacked) ---\n")

    # Violin uses |dPSI| >= 0.1 so the shape is readable
    binned_violin <- binned %>% filter(abs(deltapsi) >= 0.1)

    # ---- Top: violin of raw dPSI (|dPSI| >= 0.1) ----
    p_violin_top <- ggplot(binned_violin, aes(x = exon_bin, y = deltapsi, fill = exon_bin)) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
      geom_violin(alpha = 0.50, color = "grey30", linewidth = 0.25,
                  scale = "width", trim = FALSE, width = 0.65) +
      geom_boxplot(width = 0.10, outlier.shape = NA,
                   color = "grey20", fill = "white", alpha = 0.7) +
      stat_summary(fun = median, geom = "crossbar", width = 0.10,
                   color = "white", linewidth = 0.5) +
      # n labels
      geom_text(data = binned_violin %>%
                  group_by(exon_bin) %>%
                  summarise(n_events = n(), .groups = "drop"),
                aes(x = exon_bin,
                    y = min(binned_violin$deltapsi, na.rm = TRUE) - 0.04,
                    label = paste0("n = ", n_events)),
                size = 3, color = "grey40", family = FONT_FAMILY) +
      scale_fill_manual(values = c(
        "1\u201310"  = "#d1cbe5",
        "11\u201350" = "#9671bd",
        "51+"        = "#6a408d"
      )) +
      labs(y = expression(Delta*"PSI  (|"*Delta*"PSI| \u2265 0.1)"),
           x = NULL) +
      theme_taxol(base_size = 10) +
      theme(legend.position = "none",
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin  = margin(t = 5, r = 10, b = 0, l = 10))

    # ---- Bottom: bootstrap skewness (reuse p_skew, strip title) ----
    p_skew_bot <- p_skew +
      labs(title = NULL, subtitle = NULL) +
      theme(plot.margin = margin(t = 0, r = 10, b = 5, l = 10))

    # ---- Stack ----
    p_combined_11c <- p_violin_top / p_skew_bot +
      plot_layout(heights = c(1, 1)) +
      plot_annotation(
        title    = "Taxol Splicing Effect",
        subtitle = "Top: \u0394PSI distribution  |  Bottom: bootstrapped skewness (1000\u00d7)\n* p < 0.05   ** p < 0.01   *** p < 0.001",
        theme = theme(
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 12,
                                       family = FONT_FAMILY),
          plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey50",
                                       family = FONT_FAMILY)
        )
      )

    ggsave(file.path(PLOTS_DIR, "11c_dpsi_and_skewness.pdf"), p_combined_11c,
           device = SAVE_DEVICE, width = 5, height = 8)
    cat("  Saved: 11c_dpsi_and_skewness.pdf\n")

    cat("  Primary comparison skewness:\n")
    print(as.data.frame(boot_results %>%
      dplyr::select(exon_bin, n_events, skew_estimate, ci_lower, ci_upper, p_boot, stars)),
      row.names = FALSE)
    cat("\n")

    # 5.3 Bootstrapped skewness: ALL 4 comparisons with bootstrap dots
    cat("--- 5.3 Bootstrapped skewness: all 4 comparisons ---\n")

    all_comps <- names(comparison_labels)
    all_comp_binned <- list()

    for (tc in all_comps) {
      if (tc %in% names(fdr_results)) {
        df_tc <- fdr_results[[tc]] %>%
          filter(!is.na(deltapsi)) %>%
          mutate(
            numb_exons = as.numeric(gene_exon_counts[GENE]),
            comparison = comparison_labels[tc]
          ) %>%
          filter(!is.na(numb_exons), abs(deltapsi) > 0) %>%
          mutate(
            exon_bin = cut(numb_exons,
                           breaks = c(0, 10, 50, Inf),
                           labels = c("1\u201310", "11\u201350", "51+"),
                           right = TRUE)
          ) %>%
          filter(!is.na(exon_bin))
        all_comp_binned[[tc]] <- df_tc
      }
    }

    if (length(all_comp_binned) > 0) {
      combined_binned <- bind_rows(all_comp_binned) %>%
        mutate(comparison = factor(comparison, levels = unname(comparison_labels)))

      # Deterministic bootstrap: loop per comparison × bin with fixed seed
      comp_levels <- levels(combined_binned$comparison)
      boot_list <- list()
      seed_counter <- 0
      for (cl in comp_levels) {
        for (bl in bin_levels) {
          seed_counter <- seed_counter + 1
          bin_data <- combined_binned %>%
            filter(comparison == cl, exon_bin == bl)
          if (nrow(bin_data) < 3) next
          set.seed(2000 + seed_counter)  # unique, deterministic seed
          bs <- replicate(1000, skewness_manual(sample(bin_data$deltapsi, replace = TRUE)))
          bs <- na.omit(bs)
          if (length(bs) == 0) next
          p_val <- 2 * min(mean(bs <= 0), mean(bs >= 0))  # two-sided
          boot_list[[paste0(cl, "_", bl)]] <- tibble(
            comparison    = cl,
            exon_bin      = bl,
            n_events      = nrow(bin_data),
            skew_estimate = mean(bs),
            ci_lower      = quantile(bs, 0.025),
            ci_upper      = quantile(bs, 0.975),
            p_boot        = p_val,
            boot_skews    = list(bs)
          )
        }
      }
      boot_combined <- bind_rows(boot_list) %>%
        mutate(
          comparison = factor(comparison, levels = comp_levels),
          exon_bin   = factor(exon_bin, levels = bin_levels),
          stars      = sig_stars(p_boot)
        )

      # Unnest all 1000 bootstrap replicates for display
      boot_long_all <- boot_combined %>%
        unnest(boot_skews)

      # Auto y-limits from data range
      y_pad_all <- 0.3
      y_lim_all <- c(
        min(boot_long_all$boot_skews, na.rm = TRUE) - y_pad_all,
        max(boot_long_all$boot_skews, na.rm = TRUE) + y_pad_all
      )

      p_skew_all <- ggplot() +
        geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4, linetype = "dashed") +
        # All 1000 bootstrap dots — tiny and translucent
        geom_jitter(data = boot_long_all,
                    aes(x = exon_bin, y = boot_skews, fill = exon_bin),
                    width = 0.22, shape = 21, color = NA,
                    size = 0.3, alpha = 0.04) +
        # 95% CI error bar
        geom_errorbar(data = boot_combined,
                      aes(x = exon_bin, ymin = ci_lower, ymax = ci_upper),
                      width = 0.1, linewidth = 0.5, color = "grey20") +
        # Mean point
        geom_point(data = boot_combined,
                   aes(x = exon_bin, y = skew_estimate, fill = exon_bin),
                   shape = 21, size = 2.5, stroke = 0.6, color = "black") +
        # Significance stars (tiered)
        geom_text(data = boot_combined %>% filter(stars != ""),
                  aes(x = exon_bin, y = ci_upper + 0.1, label = stars),
                  size = 4.5, color = "#D62839", fontface = "bold",
                  family = FONT_FAMILY) +
        scale_fill_manual(values = c(
          "1\u201310"  = "#d1cbe5",
          "11\u201350" = "#9671bd",
          "51+"        = "#6a408d"
        )) +
        coord_cartesian(ylim = y_lim_all) +
        facet_wrap(~comparison, ncol = 2) +
        labs(
          title = "Bootstrapped \u0394PSI Skewness by Gene Complexity",
          subtitle = "1000 bootstraps | * p<0.05  ** p<0.01  *** p<0.001",
          x = "Number of coding exons per gene",
          y = "Skewness of \u0394PSI"
        ) +
        theme_taxol(base_size = 11) +
        theme(aspect.ratio = 1,
              legend.position = "none",
              strip.text = element_text(size = 10, face = "bold"))

      ggsave(file.path(PLOTS_DIR, "11b_skewness_all_comparisons.pdf"), p_skew_all,
             device = SAVE_DEVICE, width = 9, height = 9)
      cat("  Saved: 11b_skewness_all_comparisons.pdf\n")

      # 5.3b Stacked panel: dPSI violin + bootstrap skewness, all comparisons
      cat("--- 5.3b \u0394PSI violin + bootstrap skewness (all comps, stacked) ---\n")

      # Violin uses |dPSI| >= 0.1 so the shape is readable
      combined_binned_violin <- combined_binned %>% filter(abs(deltapsi) >= 0.1)

      # ---- Top: faceted violin of dPSI (|dPSI| >= 0.1) ----
      p_violin_top_all <- ggplot(combined_binned_violin,
                                 aes(x = exon_bin, y = deltapsi, fill = exon_bin)) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
        geom_violin(alpha = 0.50, color = "grey30", linewidth = 0.25,
                    scale = "width", trim = FALSE, width = 0.65) +
        geom_boxplot(width = 0.10, outlier.shape = NA,
                     color = "grey20", fill = "white", alpha = 0.7) +
        stat_summary(fun = median, geom = "crossbar", width = 0.10,
                     color = "white", linewidth = 0.5) +
        geom_text(data = combined_binned_violin %>%
                    group_by(comparison, exon_bin) %>%
                    summarise(n_events = n(), .groups = "drop"),
                  aes(x = exon_bin,
                      y = min(combined_binned_violin$deltapsi, na.rm = TRUE) - 0.04,
                      label = paste0("n = ", n_events)),
                  size = 2.5, color = "grey40", family = FONT_FAMILY) +
        scale_fill_manual(values = c(
          "1\u201310"  = "#d1cbe5",
          "11\u201350" = "#9671bd",
          "51+"        = "#6a408d"
        )) +
        facet_wrap(~comparison, ncol = 2) +
        labs(y = expression(Delta*"PSI  (|"*Delta*"PSI| \u2265 0.1)"),
             x = NULL) +
        theme_taxol(base_size = 10) +
        theme(legend.position = "none",
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text   = element_text(face = "bold", size = 10),
              plot.margin  = margin(t = 5, r = 10, b = 0, l = 10))

      # ---- Bottom: faceted bootstrap skewness (reuse p_skew_all, drop title/strip) ----
      p_skew_all_bot <- p_skew_all +
        labs(title = NULL, subtitle = NULL) +
        theme(strip.text  = element_blank(),
              plot.margin = margin(t = 0, r = 10, b = 5, l = 10))

      # ---- Stack ----
      p_combined_11d <- p_violin_top_all / p_skew_all_bot +
        plot_layout(heights = c(1, 1)) +
        plot_annotation(
          title    = "Splicing Skewness & Effect Size by Gene Complexity",
          subtitle = "Top: \u0394PSI distribution  |  Bottom: bootstrapped skewness (1000\u00d7)\n* p < 0.05   ** p < 0.01   *** p < 0.001",
          theme = theme(
            plot.title    = element_text(hjust = 0.5, face = "bold", size = 13,
                                         family = FONT_FAMILY),
            plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey50",
                                         family = FONT_FAMILY)
          )
        )

      ggsave(file.path(PLOTS_DIR, "11d_dpsi_and_skewness_all.pdf"), p_combined_11d,
             device = SAVE_DEVICE, width = 9, height = 14)
      cat("  Saved: 11d_dpsi_and_skewness_all.pdf\n")

      # Print summary table
      boot_summary <- boot_combined %>%
        dplyr::select(comparison, exon_bin, n_events, skew_estimate, ci_lower, ci_upper, p_boot, stars)
      cat("\n  Skewness summary (all comparisons):\n")
      print(as.data.frame(boot_summary), row.names = FALSE)
      cat("\n")
    }

  } else {
    cat("  Skipping skewness analysis: coding_unique_exon_counts_per_gene.txt not found\n")
  }
}

# ============================================================================
# ============================================================================
#  SECTION 6: EXON LENGTH ANALYSIS
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 6: EXON LENGTH & GC CONTENT ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

event_info_file <- file.path(getwd(), "EVENT_INFO-mm10.tab")
if (file.exists(event_info_file)) {

  events_info <- read.table(event_info_file, header = TRUE, sep = "\t")

  gc_content <- function(seqs) {
    sapply(seqs, function(seq) {
      if (is.na(seq) || seq == "") return(NA_real_)
      seq <- toupper(seq)
      gc <- sum(strsplit(seq, "")[[1]] %in% c("G", "C"))
      total <- nchar(seq)
      if (total == 0) return(NA_real_)
      return(100 * gc / total)
    })
  }

  # Build combined diff_exons for all 4 comparisons
  all_diff_exons <- list()
  for (comp_name in names(comparison_labels)) {
    if (!(comp_name %in% names(fdr_results))) next
    df_comp <- fdr_results[[comp_name]] %>%
      filter(!is.na(FDR), FDR <= 0.05, abs(deltapsi) >= 0.1) %>%
      left_join(events_info %>% dplyr::select(EVENT, Seq_A, LE_o), by = "EVENT") %>%
      mutate(
        LENGTH = nchar(Seq_A),
        direction = ifelse(deltapsi > 0, "Included", "Skipped"),
        comparison = comparison_labels[comp_name]
      )
    all_diff_exons[[comp_name]] <- df_comp
  }

  if (length(all_diff_exons) > 0) {
    combined_diff_exons <- bind_rows(all_diff_exons) %>%
      filter(!is.na(Seq_A), Seq_A != "") %>%
      mutate(
        comparison = factor(comparison, levels = unname(comparison_labels)),
        length_bin = cut(LENGTH, breaks = c(-Inf, 27, 50, Inf),
                         labels = c("\u226427", "28\u201350", ">50")),
        GC = gc_content(Seq_A)
      )

    # 6.1 Exon length by direction — all 4 comparisons
    cat("--- 6.1 Exon length by inclusion direction (all comparisons) ---\n")

    seq_cols_len <- c("\u226427" = "#c6dbef", "28\u201350" = "#6baed6", ">50" = "#2171b5")

    p_len <- ggplot(combined_diff_exons,
                    aes(x = length_bin, y = deltapsi, fill = length_bin)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.7) +
      stat_compare_means(
        comparisons = list(
          c("\u226427", "28\u201350"),
          c("\u226427", ">50"),
          c("28\u201350", ">50")
        ),
        method = "t.test", label = "p.signif",
        hide.ns = TRUE, tip.length = 0.02, size = 3
      ) +
      scale_fill_manual(values = seq_cols_len, name = "Length bin (nt)") +
      facet_wrap(~comparison, ncol = 2) +
      labs(
        title = "\u0394PSI by Exon Length Category",
        x = "Length bin", y = "\u0394PSI"
      ) +
      theme_taxol() +
      theme(legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"))

    ggsave(file.path(PLOTS_DIR, "12_exon_length_dpsi.pdf"), p_len,
           width = 9, height = 9, device = SAVE_DEVICE)
    cat("  Saved: 12_exon_length_dpsi.pdf\n")

    # 6.2 GC content by direction — all 4 comparisons
    cat("--- 6.2 GC content by splicing direction (all comparisons) ---\n")

    p_gc <- ggplot(combined_diff_exons,
                   aes(x = direction, y = GC, fill = direction)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
      stat_compare_means(method = "wilcox.test", size = 3.5,
                         family = FONT_FAMILY) +
      scale_fill_manual(values = c("Included" = "#D62839", "Skipped" = "#4BA3C3")) +
      facet_wrap(~comparison, ncol = 2) +
      labs(
        title = "GC Content of Differentially Spliced Exons",
        x = "Splicing Direction", y = "GC Content (%)"
      ) +
      theme_taxol() +
      theme(legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"))

    ggsave(file.path(PLOTS_DIR, "13_gc_content_by_direction.pdf"), p_gc,
           width = 9, height = 9, device = SAVE_DEVICE)
    cat("  Saved: 13_gc_content_by_direction.pdf\n")
  }
}

# ============================================================================
# ============================================================================
#  SECTION 7: HEATMAP OF DIFFERENTIALLY SPLICED EVENTS
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 7: SPLICING DYNAMICS HEATMAP\n")
cat(strrep("=", 70), "\n\n")

# Collect all significant events across comparisons
all_sig_events <- unique(unlist(lapply(names(fdr_results), get_sig_events)))
cat(sprintf("  Total unique significant events across all comparisons: %d\n",
            length(all_sig_events)))

if (length(all_sig_events) > 5) {
  # Build PSI matrix from the raw data (same structure used in PCA/distribution)
  all_filtered <- filterEvents(all_events, types = c("C1","C2","C3","S","MIC","IR"), N = 10)
  heatmap_mat <- all_filtered$PSI %>%
    dplyr::select(EVENT, all_of(metadata$sample_id))
  heatmap_mat <- heatmap_mat[heatmap_mat$EVENT %in% all_sig_events, ]

  # Convert to numeric matrix
  event_ids <- heatmap_mat$EVENT
  heatmap_mat_num <- heatmap_mat %>%
    dplyr::select(-EVENT) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  rownames(heatmap_mat_num) <- event_ids

  # Remove rows with all NA
  heatmap_mat_num <- heatmap_mat_num[complete.cases(heatmap_mat_num), ]

  if (nrow(heatmap_mat_num) > 5) {
    cat(sprintf("  Matrix for heatmap: %d events x %d samples\n",
                nrow(heatmap_mat_num), ncol(heatmap_mat_num)))

    # Column dendrogram
    col_dend <- hclust(dist(t(heatmap_mat_num))) %>%
      as.dendrogram() %>%
      color_branches(k = 4)

    # Color function
    col_fun <- colorRamp2(c(0, 50, 100), viridis(3))

    # Column annotation
    ha_col <- HeatmapAnnotation(
      Condition = metadata$condition[match(colnames(heatmap_mat_num), metadata$sample_id)],
      col = list(Condition = condition_colors),
      annotation_name_gp = gpar(fontsize = 10)
    )

    # Determine n clusters (min 2, max 6)
    n_clust <- min(6, max(2, round(nrow(heatmap_mat_num) / 30)))

    ht <- Heatmap(
      heatmap_mat_num,
      name = "PSI",
      col  = col_fun,
      clustering_distance_rows = "pearson",
      cluster_columns  = col_dend,
      top_annotation   = ha_col,
      column_title     = "C2C12 Taxol Splicing Analysis",
      column_title_gp  = gpar(fontsize = 16, fontface = "bold"),
      row_title        = "Splicing Events",
      row_title_gp     = gpar(fontsize = 12, fontface = "italic"),
      column_labels    = metadata$condition[match(colnames(heatmap_mat_num), metadata$sample_id)],
      column_names_rot = 45,
      column_names_gp  = gpar(fontsize = 10),
      row_km           = n_clust,
      show_row_names   = FALSE,
      use_raster       = FALSE,
      rect_gp          = gpar(col = "grey90", lwd = 0.3),
      row_gap          = unit(3, "mm"),
      column_gap       = unit(2, "mm"),
      heatmap_legend_param = list(
        title      = "PSI",
        title_gp   = gpar(fontsize = 12, fontface = "bold"),
        labels_gp  = gpar(fontsize = 9)
      )
    )

    heatmap_file <- file.path(PLOTS_DIR, "14_heatmap_splicing_dynamics.pdf")
    pdf(heatmap_file, width = 10, height = 12)
    draw(ht, heatmap_legend_side = "right")
    dev.off()
    if (file.exists(heatmap_file) && file.size(heatmap_file) > 0) {
      cat("  Saved: 14_heatmap_splicing_dynamics.pdf\n")
    } else {
      cat("  WARNING: heatmap PDF was not written correctly\n")
    }
    cat("  Saved: 14_heatmap_splicing_dynamics.pdf\n")
  }
}

# ============================================================================
# ============================================================================
#  SECTION 8: GENE ONTOLOGY ENRICHMENT
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 8: GO ENRICHMENT ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

# Enrichment for each comparison
run_enrichment <- function(sig_genes, all_genes, comp_label) {
  if (length(sig_genes) < 5) {
    cat(sprintf("  Skipping %s: only %d significant genes\n", comp_label, length(sig_genes)))
    return(NULL)
  }

  cat(sprintf("  %s: %d significant genes\n", comp_label, length(sig_genes)))

  # Convert to ENTREZ IDs
  entrez_sig <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
  }, error = function(e) character(0))

  entrez_bg <- tryCatch({
    bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
  }, error = function(e) character(0))

  if (length(entrez_sig) < 3) return(NULL)

  ego <- tryCatch({
    enrichGO(
      gene          = entrez_sig,
      universe      = entrez_bg,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
  }, error = function(e) NULL)

  return(ego)
}

plot_enrichment <- function(ego, comp_label) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)

  df <- as.data.frame(ego) %>%
    mutate(
      log_pval = -log10(pvalue),
      Term = str_remove(Description, "\\s*\\(GO:\\d+\\)"),
      Count = as.numeric(Count)
    ) %>%
    arrange(pvalue) %>%
    slice_head(n = 15) %>%
    mutate(Term = fct_reorder(Term, log_pval))

  ggplot(df, aes(x = log_pval, y = Term, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "#2166AC", high = "#B2182B", trans = "reverse") +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      x     = expression(-log[10](p-value)),
      y     = NULL,
      color = "Adj. P-value",
      size  = "Gene Count",
      title = sprintf("GO BP: %s", comp_label)
    ) +
    theme_taxol() +
    theme(
      axis.text.y  = element_text(size = rel(0.8)),
      legend.position = "right",
      aspect.ratio = NULL
    )
}

# Background genes from all events
all_bg_genes <- unique(all_events$PSI$GENE)

# Run enrichment for each comparison
for (comp_name in names(comparison_labels)) {
  if (comp_name %in% names(fdr_results)) {
    sig_genes <- get_sig_genes(comp_name)
    ego <- run_enrichment(sig_genes, all_bg_genes, comparison_labels[comp_name])

    if (!is.null(ego)) {
      p <- plot_enrichment(ego, comparison_labels[comp_name])
      if (!is.null(p)) {
        fname <- sprintf("15_go_enrichment_%s.pdf", comp_name)
        ggsave(file.path(PLOTS_DIR, fname), p, width = 9, height = 6, device = SAVE_DEVICE)
        cat(sprintf("    Saved: %s\n", fname))
      }
    }
  }
}

# EnrichR analysis (broader databases)
cat("\n--- EnrichR analysis ---\n")
for (comp_name in names(comparison_labels)) {
  sig_genes_all <- get_sig_genes(comp_name)

  if (length(sig_genes_all) >= 5) {
    enrichr_res <- tryCatch({
      enrichr(sig_genes_all, databases = c(
        "GO_Biological_Process_2023",
        "GO_Cellular_Component_2023",
        "GO_Molecular_Function_2023",
        "KEGG_2021_Human",
        "Reactome_2022"
      ))
    }, error = function(e) NULL)

    if (!is.null(enrichr_res)) {
      for (db_name in names(enrichr_res)) {
        outfile <- file.path(RESULTS_DIR,
                             sprintf("enrichr_%s_%s.csv", comp_name,
                                     gsub(" ", "_", db_name)))
        write.csv(enrichr_res[[db_name]], outfile, row.names = FALSE)
      }
      cat(sprintf("  EnrichR saved for %s (%d genes)\n",
                  comparison_labels[comp_name], length(sig_genes_all)))
    }
  }
}

# ============================================================================
# ============================================================================
#  SECTION 9: TAXOL-SPECIFIC INTERACTION ANALYSIS
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 9: TAXOL-DIFFERENTIATION INTERACTION\n")
cat(strrep("=", 70), "\n\n")

# Compare differentiation with and without Taxol
# Which differentiation-associated splicing events are affected by Taxol?

if (all(c("Differentiation_DMSO", "Differentiation_Taxol") %in%
        names(fdr_results))) {

  diff_dmso  <- fdr_results[["Differentiation_DMSO"]] %>%
    dplyr::select(EVENT, GENE, dpsi_dmso = deltapsi, FDR_dmso = FDR)
  diff_taxol <- fdr_results[["Differentiation_Taxol"]] %>%
    dplyr::select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)

  interaction <- inner_join(diff_dmso, diff_taxol, by = "EVENT") %>%
    mutate(
      sig_dmso  = !is.na(FDR_dmso)  & FDR_dmso  <= 0.05 & abs(dpsi_dmso)  >= 0.1,
      sig_taxol = !is.na(FDR_taxol) & FDR_taxol <= 0.05 & abs(dpsi_taxol) >= 0.1,
      category = case_when(
        sig_dmso & sig_taxol & sign(dpsi_dmso) == sign(dpsi_taxol) ~ "Preserved",
        sig_dmso & !sig_taxol ~ "Disrupted by Taxol",
        sig_dmso & sig_taxol & sign(dpsi_dmso) != sign(dpsi_taxol) ~ "Reversed by Taxol",
        !sig_dmso & sig_taxol ~ "Taxol-induced",
        TRUE ~ "NS"
      )
    )

  interaction_colors <- c(
    "Preserved"         = "#2166AC",
    "Disrupted by Taxol"= "#B2182B",
    "Reversed by Taxol" = "#762A83",
    "Taxol-induced"     = "#1B7837",
    "NS"                = "grey85"
  )

  p_interaction <- ggplot(interaction,
                          aes(x = dpsi_dmso, y = dpsi_taxol, color = category)) +
    geom_point(data = filter(interaction, category == "NS"), alpha = 0.1, size = 0.5) +
    geom_point(data = filter(interaction, category != "NS"), alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted",
                linewidth = 0.6) +
    scale_color_manual(values = interaction_colors,
                       breaks = names(interaction_colors)[1:4]) +
    labs(
      title    = "Taxol Impact on Differentiation-Associated Splicing",
      subtitle = "Does Taxol disrupt normal differentiation splicing programs?",
      x = "dPSI Differentiation (DMSO)",
      y = "dPSI Differentiation (Taxol)",
      color = NULL
    ) +
    coord_fixed(ratio = 1, xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_taxol() +
    theme(legend.position = "right")

  ggsave(file.path(PLOTS_DIR, "16_interaction_scatter.pdf"), p_interaction,
         width = 7, height = 6, device = SAVE_DEVICE)
  cat("  Saved: 16_interaction_scatter.pdf\n")

  # Stacked bar of categories
  interaction_summary <- interaction %>%
    filter(category != "NS") %>%
    dplyr::count(category) %>%
    mutate(
      pct = 100 * n / sum(n),
      category = factor(category,
                        levels = c("Taxol-induced", "Reversed by Taxol",
                                   "Disrupted by Taxol", "Preserved"))
    )

  p_int_bar <- ggplot(interaction_summary,
                      aes(x = "All Events", y = pct, fill = category)) +
    geom_col(width = 0.5, color = "grey30", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold",
              family = FONT_FAMILY) +
    scale_fill_manual(values = interaction_colors) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(
      title = "How Taxol Affects Differentiation Splicing",
      x = NULL, y = "% of Significant Events", fill = NULL
    ) +
    theme_taxol() +
    theme(legend.position = "right", axis.text.x = element_blank())

  ggsave(file.path(PLOTS_DIR, "17_interaction_barplot.pdf"), p_int_bar,
         width = 5, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 17_interaction_barplot.pdf\n")

  # Save the interaction table
  write.csv(interaction %>% filter(category != "NS"),
            file.path(RESULTS_DIR, "interaction_differentiation_taxol.csv"),
            row.names = FALSE)
}

# ============================================================================
# ============================================================================
#  SECTION 10: EVENT TYPE BREAKDOWN
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 10: EVENT TYPE BREAKDOWN\n")
cat(strrep("=", 70), "\n\n")

# For each comparison, break down significant events by type (C1, C2, C3, S, MIC, IR)
type_breakdown_list <- list()

for (comp_name in names(comparison_labels)) {
  if (comp_name %in% names(fdr_results)) {
    sig <- fdr_results[[comp_name]] %>%
      filter(!is.na(FDR), FDR <= 0.05, abs(deltapsi) >= 0.1) %>%
      mutate(
        event_subtype = case_when(
          grepl("^MmuEX", EVENT) ~ "EX",
          grepl("^MmuINT", EVENT) ~ "IR",
          grepl("^MmuALT3", EVENT) ~ "Alt3",
          grepl("^MmuALT5", EVENT) ~ "Alt5",
          TRUE ~ "Other"
        ),
        comparison = comparison_labels[comp_name]
      )
    type_breakdown_list[[comp_name]] <- sig
  }
}

if (length(type_breakdown_list) > 0) {
  type_df <- bind_rows(type_breakdown_list)

  type_summary <- type_df %>%
    dplyr::count(comparison, event_subtype) %>%
    group_by(comparison) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  type_colors <- c(
    "EX"    = "#4BA3C3",
    "IR"    = "#D62839",
    "Alt3"  = "#762A83",
    "Alt5"  = "#1B7837",
    "Other" = "grey60"
  )

  p_types <- ggplot(type_summary,
                    aes(x = comparison, y = n, fill = event_subtype)) +
    geom_col(position = "dodge", width = 0.7, color = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = type_colors) +
    labs(
      title = "Breakdown of Significant Events by Type",
      x = NULL, y = "Number of Events", fill = "Event Type"
    ) +
    theme_taxol() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = rel(0.8)))

  ggsave(file.path(PLOTS_DIR, "18_event_type_breakdown.pdf"), p_types,
         width = 9, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 18_event_type_breakdown.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 11: SPLICING FACTOR CORRELATION
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 11: SPLICING FACTOR CORRELATION ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

# Known splicing factors list (common ones in muscle differentiation)
known_sfs <- c("Rbfox2", "Rbfox1", "Mbnl1", "Mbnl2", "Ptbp1", "Ptbp2",
               "Srsf1", "Srsf2", "Srsf3", "Srsf5", "Srsf6", "Srsf7",
               "Hnrnpa1", "Hnrnpa2b1", "Hnrnph1", "Hnrnpk",
               "Qki", "Esrp1", "Esrp2", "Celf1", "Celf2",
               "Tra2a", "Tra2b", "Nova1", "Nova2", "Khdrbs1",
               "Rbm3", "Rbm4", "Rbm10", "Rbm25")

# Check which SFs have significant events
all_sig_genes_combined <- unique(unlist(lapply(names(fdr_results), get_sig_genes)))
sf_with_events <- known_sfs[tolower(known_sfs) %in% tolower(all_sig_genes_combined)]

if (length(sf_with_events) > 0) {
  cat(sprintf("  Splicing factors with differential events: %s\n",
              paste(sf_with_events, collapse = ", ")))

  # Build correlation matrix for SF events
  all_filtered <- filterEvents(all_events, types = c("C1","C2","C3","S","MIC","IR"), N = 10)

  # Get events for SFs and non-SFs
  sf_events <- all_filtered$PSI %>%
    filter(tolower(GENE) %in% tolower(sf_with_events)) %>%
    filter(EVENT %in% all_sig_events)

  non_sf_events <- all_filtered$PSI %>%
    filter(!tolower(GENE) %in% tolower(known_sfs)) %>%
    filter(EVENT %in% all_sig_events)

  if (nrow(sf_events) > 1 && nrow(non_sf_events) > 1) {
    sf_mat  <- as.matrix(sf_events[, metadata$sample_id])
    rownames(sf_mat) <- sf_events$EVENT
    nsf_mat <- as.matrix(non_sf_events[, metadata$sample_id])
    rownames(nsf_mat) <- non_sf_events$EVENT

    # Pearson correlation between SFs and all other events
    results_corr <- data.frame(
      SF_Gene = character(), SF_Event = character(),
      Target_Gene = character(), Target_Event = character(),
      Correlation = numeric(), P_value = numeric(),
      stringsAsFactors = FALSE
    )

    for (i in 1:nrow(sf_mat)) {
      for (j in 1:nrow(nsf_mat)) {
        if (all(!is.na(sf_mat[i, ])) && all(!is.na(nsf_mat[j, ]))) {
          ct <- cor.test(sf_mat[i, ], nsf_mat[j, ], method = "pearson")
          if (!is.na(ct$p.value)) {
            results_corr <- rbind(results_corr, data.frame(
              SF_Gene = sf_events$GENE[i], SF_Event = rownames(sf_mat)[i],
              Target_Gene = non_sf_events$GENE[j], Target_Event = rownames(nsf_mat)[j],
              Correlation = ct$estimate, P_value = ct$p.value
            ))
          }
        }
      }
    }

    if (nrow(results_corr) > 0) {
      results_corr$Adjusted_P <- p.adjust(results_corr$P_value, method = "BH")
      sig_corr <- results_corr %>%
        filter(Adjusted_P < 0.01, abs(Correlation) >= 0.85)

      write.csv(sig_corr, file.path(RESULTS_DIR, "sf_correlations.csv"),
                row.names = FALSE)
      cat(sprintf("  Found %d significant SF-target correlations (|r| >= 0.85, FDR < 0.01)\n",
                  nrow(sig_corr)))

      if (nrow(sig_corr) > 0) {
        # Cumulative interaction plot
        cum_data <- sig_corr %>%
          dplyr::count(SF_Gene, name = "n_targets") %>%
          arrange(desc(n_targets))

        p_sf <- ggplot(cum_data, aes(x = reorder(SF_Gene, -n_targets), y = n_targets,
                                     fill = SF_Gene)) +
          geom_col(width = 0.7, color = "grey30", linewidth = 0.3, show.legend = FALSE) +
          geom_text(aes(label = n_targets), vjust = -0.3, family = FONT_FAMILY, size = 3.5) +
          scale_fill_viridis_d(option = "plasma") +
          labs(
            title = "Splicing Factor Target Interactions",
            subtitle = "Pearson |r| \u2265 0.85, FDR < 0.01",
            x = "Splicing Factor", y = "Number of Correlated Events"
          ) +
          theme_taxol() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

        ggsave(file.path(PLOTS_DIR, "19_sf_interactions.pdf"), p_sf, device = SAVE_DEVICE,
               width = 7, height = 5)
        cat("  Saved: 19_sf_interactions.pdf\n")
      }
    }
  }
} else {
  cat("  No known splicing factors found in significant events\n")
}

# ============================================================================
# ============================================================================
#  SECTION 12: COMBINED TAXOL EFFECT FIGURE
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 12: COMBINED FIGURE\n")
cat(strrep("=", 70), "\n\n")

# Build a 2x2 combined figure similar to q3 blocking analysis

if (all(c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO") %in%
        names(fdr_results))) {

  p1 <- make_volcano(fdr_results[["Myoblast_Taxol_vs_DMSO"]],
                     "Myoblast: Taxol vs DMSO")
  p2 <- make_volcano(fdr_results[["Myotube_Taxol_vs_DMSO"]],
                     "Myotube: Taxol vs DMSO")

  # Additional panels if differentiation data exists
  extra_panels <- list()
  if ("Differentiation_DMSO" %in% names(fdr_results)) {
    extra_panels[["p3"]] <- make_volcano(
      fdr_results[["Differentiation_DMSO"]],
      "Differentiation (DMSO)")
  }
  if ("Differentiation_Taxol" %in% names(fdr_results)) {
    extra_panels[["p4"]] <- make_volcano(
      fdr_results[["Differentiation_Taxol"]],
      "Differentiation (Taxol)")
  }

  if (length(extra_panels) == 2) {
    combined <- (p1 + p2) / (extra_panels[["p3"]] + extra_panels[["p4"]]) +
      plot_layout(widths = c(1, 1), heights = c(1, 1)) +
      plot_annotation(
        title    = "Taxol Effect on Alternative Splicing in C2C12 Cells",
        subtitle = "All comparisons: FDR \u2264 0.05, |\u0394PSI| \u2265 0.1",
        theme    = theme(
          plot.title    = element_text(size = 14, face = "bold", hjust = 0.5,
                                       family = FONT_FAMILY),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40",
                                       family = FONT_FAMILY)
        )
      )
  } else {
    combined <- p1 + p2 +
      plot_annotation(
        title    = "Taxol Effect on Exon Splicing in C2C12",
        theme    = theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                                    family = FONT_FAMILY)
        )
      )
  }

  ggsave(file.path(PLOTS_DIR, "20_combined_volcano.pdf"), combined,
         width = 10, height = 10, device = SAVE_DEVICE)
  cat("  Saved: 20_combined_volcano.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 13: NUCLEAR SPECKLE ENRICHMENT ANALYSIS
# ============================================================================
# ============================================================================
#
# HYPOTHESIS: Taxol stabilises microtubules → alters cytoplasmic forces →
#   changes nuclear envelope mechanics → perturbs nucleoplasmic agitation →
#   differentially affects splicing at nuclear speckles / speckle-associated
#   genes. If true, taxol-responsive splicing events should be enriched in
#   genes whose mRNAs are processed at or near nuclear speckles.
#
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 13: NUCLEAR SPECKLE GENE ENRICHMENT\n")
cat(strrep("=", 70), "\n\n")

speckle_file <- file.path(getwd(), "nuclear_speckle_associated_genes.txt")

if (file.exists(speckle_file)) {
  speckle_genes_raw <- readLines(speckle_file)
  speckle_genes_raw <- speckle_genes_raw[speckle_genes_raw != ""]
  cat(sprintf("  Loaded %d nuclear speckle-associated genes\n",
              length(speckle_genes_raw)))

  # The speckle gene list is likely human symbols — convert to mouse equivalents
  # Mouse gene symbols: first letter upper, rest lower
  speckle_genes_mouse <- str_to_title(tolower(speckle_genes_raw))
  # Keep original too for matching
  speckle_genes_all <- unique(c(speckle_genes_raw, speckle_genes_mouse))

  # Background: all genes in the dataset
  all_dataset_genes <- unique(all_events$PSI$GENE)

  # 13.1 Test enrichment for EACH taxol comparison
  taxol_comparisons <- c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO")
  diff_comparisons  <- c("Differentiation_DMSO", "Differentiation_Taxol")
  all_comps <- c(taxol_comparisons, diff_comparisons)

  speckle_enrichment_results <- list()

  for (comp_name in all_comps) {
    if (!(comp_name %in% names(fdr_results))) next
    sig_genes <- get_sig_genes(comp_name)

    # 2x2 contingency table
    # Sig & Speckle | Sig & NOT Speckle
    # NS  & Speckle | NS  & NOT Speckle
    in_speckle_sig  <- sum(sig_genes %in% speckle_genes_all)
    in_speckle_bg   <- sum(all_dataset_genes %in% speckle_genes_all)
    not_speckle_sig <- length(sig_genes) - in_speckle_sig
    not_speckle_bg  <- length(all_dataset_genes) - in_speckle_bg

    mat_2x2 <- matrix(c(
      in_speckle_sig,
      not_speckle_sig,
      in_speckle_bg  - in_speckle_sig,
      not_speckle_bg - not_speckle_sig
    ), nrow = 2, byrow = TRUE,
    dimnames = list(c("Significant", "Not Significant"),
                    c("Speckle", "Not Speckle")))

    ft <- fisher.test(mat_2x2, alternative = "greater")

    speckle_enrichment_results[[comp_name]] <- data.frame(
      comparison   = comparison_labels[comp_name],
      n_sig        = length(sig_genes),
      n_speckle_sig = in_speckle_sig,
      pct_speckle   = round(100 * in_speckle_sig / max(length(sig_genes), 1), 1),
      bg_pct_speckle = round(100 * in_speckle_bg / length(all_dataset_genes), 1),
      odds_ratio    = round(ft$estimate, 3),
      p_value       = ft$p.value,
      stringsAsFactors = FALSE
    )

    cat(sprintf("  %s: %d/%d sig genes in speckle list (%.1f%% vs %.1f%% bg), OR=%.2f, p=%.4g\n",
                comparison_labels[comp_name], in_speckle_sig, length(sig_genes),
                100 * in_speckle_sig / max(length(sig_genes), 1),
                100 * in_speckle_bg / length(all_dataset_genes),
                ft$estimate, ft$p.value))
  }

  speckle_enrich_df <- bind_rows(speckle_enrichment_results)

  # Save results
  write.csv(speckle_enrich_df,
            file.path(RESULTS_DIR, "nuclear_speckle_enrichment.csv"),
            row.names = FALSE)

  # 13.2 Barplot: % speckle genes in each comparison vs background
  speckle_plot_df <- speckle_enrich_df %>%
    dplyr::select(comparison, pct_speckle, bg_pct_speckle, p_value, odds_ratio) %>%
    pivot_longer(cols = c(pct_speckle, bg_pct_speckle),
                 names_to = "group", values_to = "pct") %>%
    mutate(
      group = ifelse(group == "pct_speckle", "Significant Events", "Background"),
      sig_label = ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                  ifelse(p_value < 0.05, "*", "ns")))
    )

  p_speckle <- ggplot(speckle_plot_df,
                      aes(x = comparison, y = pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6,
             color = "grey30", linewidth = 0.3) +
    scale_fill_manual(values = c("Significant Events" = "#D62839",
                                 "Background" = "grey70")) +
    labs(
      title = "Nuclear Speckle Gene Enrichment",
      subtitle = "Are taxol-affected genes enriched for speckle association?",
      x = NULL, y = "% Genes in Speckle List", fill = NULL
    ) +
    theme_taxol() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = rel(0.85)),
          legend.position = "top") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  ggsave(file.path(PLOTS_DIR, "21_nuclear_speckle_enrichment.pdf"), p_speckle,
         width = 8, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 21_nuclear_speckle_enrichment.pdf\n")

  # 13.3 Odds ratio forest plot
  speckle_forest <- speckle_enrich_df %>%
    mutate(
      log2OR = log2(pmax(odds_ratio, 0.01)),
      sig = p_value < 0.05,
      comparison = factor(comparison,
                          levels = rev(comparison))
    )

  p_forest <- ggplot(speckle_forest,
                     aes(x = log2OR, y = comparison, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 4) +
    geom_text(aes(label = sprintf("p=%.3g", p_value)),
              hjust = -0.2, size = 3, family = FONT_FAMILY) +
    scale_color_manual(values = c(`TRUE` = "#D62839", `FALSE` = "grey60"),
                       labels = c("p < 0.05", "ns"), name = NULL) +
    labs(
      title = "Nuclear Speckle Enrichment (Odds Ratios)",
      subtitle = "Fisher's exact test, one-sided (greater)",
      x = expression(log[2](Odds~Ratio)), y = NULL
    ) +
    theme_taxol() +
    theme(legend.position = "top")

  ggsave(file.path(PLOTS_DIR, "22_speckle_odds_ratio.pdf"), p_forest,
         width = 7, height = 4, device = SAVE_DEVICE)
  cat("  Saved: 22_speckle_odds_ratio.pdf\n")

  # 13.4 dPSI distribution: speckle vs non-speckle genes (all 4 comparisons, 2×2)
  speckle_dpsi_list <- list()
  for (tc in all_comps) {
    if (!(tc %in% names(fdr_results))) next
    df_tc <- fdr_results[[tc]] %>%
      filter(!is.na(deltapsi)) %>%
      mutate(
        speckle = ifelse(GENE %in% speckle_genes_all,
                         "Speckle-Associated", "Other"),
        comparison = comparison_labels[tc]
      )
    speckle_dpsi_list[[tc]] <- df_tc
  }

  if (length(speckle_dpsi_list) > 0) {
    speckle_dpsi_all <- bind_rows(speckle_dpsi_list) %>%
      mutate(comparison = factor(comparison, levels = unname(comparison_labels)))

    p_speckle_dpsi <- ggplot(speckle_dpsi_all,
                             aes(x = speckle, y = deltapsi, fill = speckle)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.1, size = 0.5) +
      stat_compare_means(method = "wilcox.test", size = 3, family = FONT_FAMILY) +
      scale_fill_manual(values = c("Speckle-Associated" = "#762A83",
                                   "Other" = "grey75")) +
      facet_wrap(~comparison, ncol = 2) +
      labs(
        title = "\u0394PSI: Speckle vs Non-Speckle Genes",
        x = NULL, y = "\u0394PSI"
      ) +
      theme_taxol() +
      theme(legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"))

    ggsave(file.path(PLOTS_DIR, "23_speckle_dpsi_all.pdf"), p_speckle_dpsi,
           width = 9, height = 9, device = SAVE_DEVICE)
    cat("  Saved: 23_speckle_dpsi_all.pdf\n")
  }

  # 13.5 Among significant events: are speckle genes more skipped or included? (all 4 comps, 2×2)
  speckle_dir_list <- list()
  for (tc in all_comps) {
    if (!(tc %in% names(fdr_results))) next
    sig_tc <- fdr_results[[tc]] %>%
      filter(!is.na(FDR), FDR <= 0.05, abs(deltapsi) >= 0.1) %>%
      mutate(
        speckle    = ifelse(GENE %in% speckle_genes_all, "Speckle", "Other"),
        direction  = ifelse(deltapsi > 0, "Included", "Skipped"),
        comparison = comparison_labels[tc]
      )
    if (nrow(sig_tc) > 5) speckle_dir_list[[tc]] <- sig_tc
  }

  if (length(speckle_dir_list) > 0) {
    speckle_dir_all <- bind_rows(speckle_dir_list) %>%
      mutate(comparison = factor(comparison, levels = unname(comparison_labels)))

    speckle_dir_summary <- speckle_dir_all %>%
      dplyr::count(comparison, speckle, direction) %>%
      group_by(comparison, speckle) %>%
      mutate(pct = 100 * n / sum(n)) %>%
      ungroup()

    p_dir <- ggplot(speckle_dir_summary,
                    aes(x = speckle, y = pct, fill = direction)) +
      geom_col(position = "stack", width = 0.6,
               color = "grey30", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
                position = position_stack(vjust = 0.5),
                size = 2.5, family = FONT_FAMILY) +
      scale_fill_manual(values = sig_colors[1:2]) +
      facet_wrap(~comparison, ncol = 2) +
      labs(
        title = "Inclusion vs Skipping in Speckle Genes",
        x = NULL, y = "% of Significant Events", fill = NULL
      ) +
      theme_taxol() +
      theme(legend.position = "top",
            strip.text = element_text(size = 10, face = "bold")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 105))

    ggsave(file.path(PLOTS_DIR, "24_speckle_direction_all.pdf"), p_dir,
           width = 9, height = 9, device = SAVE_DEVICE)
    cat("  Saved: 24_speckle_direction_all.pdf\n")
  }
} else {
  cat("  nuclear_speckle_associated_genes.txt not found — skipping\n")
}

# ============================================================================
# ============================================================================
#  SECTION 14: FORCE-PREDICTED GENES OVERLAP
# ============================================================================
# ============================================================================
#
# Test overlap of taxol-affected splicing genes with genes previously
# predicted to be affected by cytoplasmic/nuclear forces (from mechanical
# perturbation analysis).
#
# final_all_layers.csv       -> genes affected by forces (all contexts)
# final_myo_forces_layers.csv -> genes affected by force + differentiation
#
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 14: FORCE-PREDICTED GENES OVERLAP\n")
cat(strrep("=", 70), "\n\n")

force_all_file <- file.path(getwd(), "final_all_layers.csv")
force_myo_file <- file.path(getwd(), "final_myo_forces_layers.csv")

if (file.exists(force_all_file) && file.exists(force_myo_file)) {
  force_all <- read.csv(force_all_file, stringsAsFactors = FALSE)
  force_myo <- read.csv(force_myo_file, stringsAsFactors = FALSE)

  force_all_genes <- unique(force_all$GENE)
  force_myo_genes <- unique(force_myo$GENE)

  cat(sprintf("  Force-all genes:    %d unique genes (%d events)\n",
              length(force_all_genes), nrow(force_all)))
  cat(sprintf("  Force-myo genes:    %d unique genes (%d events)\n",
              length(force_myo_genes), nrow(force_myo)))

  # 14.1 Hypergeometric overlap tests
  all_dataset_genes_force <- unique(all_events$PSI$GENE)
  N_bg <- length(all_dataset_genes_force)

  overlap_results <- list()

  for (comp_name in names(comparison_labels)) {
    if (!(comp_name %in% names(fdr_results))) next
    sig_genes_comp <- get_sig_genes(comp_name)
    n_sig <- length(sig_genes_comp)

    for (force_label in c("Force-All", "Force-Myo")) {
      force_genes <- if (force_label == "Force-All") force_all_genes else force_myo_genes
      n_force <- sum(force_genes %in% all_dataset_genes_force)
      n_overlap <- sum(sig_genes_comp %in% force_genes)

      # Hypergeometric test (phyper: one-sided, enrichment)
      p_hyper <- phyper(n_overlap - 1, n_force, N_bg - n_force, n_sig,
                        lower.tail = FALSE)

      expected <- n_sig * n_force / N_bg
      fold_enrich <- n_overlap / max(expected, 0.01)

      overlap_results[[paste(comp_name, force_label)]] <- data.frame(
        comparison     = comparison_labels[comp_name],
        force_set      = force_label,
        n_sig          = n_sig,
        n_force_in_bg  = n_force,
        n_overlap      = n_overlap,
        expected       = round(expected, 2),
        fold_enrichment = round(fold_enrich, 2),
        p_value        = p_hyper,
        stringsAsFactors = FALSE
      )

      cat(sprintf("  %s x %s: %d overlap (expected %.1f), fold=%.2f, p=%.4g\n",
                  comparison_labels[comp_name], force_label,
                  n_overlap, expected, fold_enrich, p_hyper))
    }
  }

  overlap_df <- bind_rows(overlap_results)
  write.csv(overlap_df,
            file.path(RESULTS_DIR, "force_gene_overlap.csv"),
            row.names = FALSE)

  # 14.2 Visualize fold enrichment
  overlap_df$neg_log_p <- -log10(overlap_df$p_value)
  overlap_df$sig <- overlap_df$p_value < 0.05

  p_overlap <- ggplot(overlap_df,
                      aes(x = fold_enrichment, y = comparison,
                          color = force_set, size = neg_log_p, shape = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c("Force-All" = "#E66101", "Force-Myo" = "#5E3C99")) +
    scale_size_continuous(range = c(2, 8), name = expression(-log[10](p))) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c("p < 0.05", "ns"), name = NULL) +
    labs(
      title = "Overlap with Force-Predicted Genes",
      subtitle = "Hypergeometric test for enrichment",
      x = "Fold Enrichment over Expected", y = NULL,
      color = "Force Gene Set"
    ) +
    theme_taxol() +
    theme(legend.position = "right")

  ggsave(file.path(PLOTS_DIR, "25_force_gene_overlap.pdf"), p_overlap,
         width = 8, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 25_force_gene_overlap.pdf\n")

  # 14.3 Venn-style overlap for taxol comparisons
  for (tc in c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO")) {
    if (!(tc %in% names(fdr_results))) next
    sig_genes_tc <- get_sig_genes(tc)

    overlap_list <- list(
      "Taxol Sig"   = sig_genes_tc,
      "Force-All"   = force_all_genes[force_all_genes %in% all_dataset_genes_force],
      "Force-Myo"   = force_myo_genes[force_myo_genes %in% all_dataset_genes_force]
    )

    nonempty_ov <- sapply(overlap_list, length) > 0
    if (sum(nonempty_ov) >= 2) {
      open_pdf(file.path(PLOTS_DIR, sprintf("26_force_upset_%s.pdf", tc)),
               width = 8, height = 5)
      binary_ov <- fromList(overlap_list[nonempty_ov])
      UpSetR::upset(
        binary_ov,
        sets = names(overlap_list[nonempty_ov]),
        sets.bar.color = c("#D62839", "#E66101", "#5E3C99")[nonempty_ov],
        order.by = "freq",
        mainbar.y.label = "Gene Overlap",
        sets.x.label = "Total Genes",
        keep.order = TRUE,
        text.scale = c(1.3, 1.2, 1.0, 1.0, 1.3, 1.2),
        mb.ratio = c(0.65, 0.35)
      )
      grid.text(sprintf("Taxol Sig Genes vs Force-Predicted (%s)",
                        comparison_labels[tc]),
                x = 0.5, y = 0.97,
                gp = gpar(fontsize = 11, fontface = "bold", fontfamily = FONT_FAMILY))
      dev.off()
      cat(sprintf("  Saved: 26_force_upset_%s.pdf\n", tc))
    }
  }

  # 14.4 dPSI comparison: force-predicted vs other genes (all 4 comparisons, 2×2)
  force_dpsi_list <- list()
  for (tc in names(comparison_labels)) {
    if (!(tc %in% names(fdr_results))) next
    df_force <- fdr_results[[tc]] %>%
      filter(!is.na(deltapsi)) %>%
      mutate(
        force_category = case_when(
          GENE %in% force_all_genes & GENE %in% force_myo_genes ~ "Force-All + Myo",
          GENE %in% force_all_genes ~ "Force-All Only",
          GENE %in% force_myo_genes ~ "Force-Myo Only",
          TRUE ~ "Not Force-Predicted"
        ),
        force_category = factor(force_category,
                                levels = c("Not Force-Predicted", "Force-All Only",
                                           "Force-Myo Only", "Force-All + Myo")),
        comparison = comparison_labels[tc]
      )
    force_dpsi_list[[tc]] <- df_force
  }

  if (length(force_dpsi_list) > 0) {
    force_dpsi_all <- bind_rows(force_dpsi_list) %>%
      mutate(comparison = factor(comparison, levels = unname(comparison_labels)))

    force_cat_cols <- c(
      "Not Force-Predicted" = "grey75",
      "Force-All Only"      = "#E66101",
      "Force-Myo Only"      = "#5E3C99",
      "Force-All + Myo"     = "#B2182B"
    )

    p_force_dpsi <- ggplot(force_dpsi_all,
                           aes(x = force_category, y = abs(deltapsi),
                               fill = force_category)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.65) +
      geom_jitter(width = 0.15, alpha = 0.08, size = 0.4) +
      stat_compare_means(ref.group = "Not Force-Predicted",
                         method = "wilcox.test", label = "p.signif",
                         hide.ns = TRUE, size = 3) +
      scale_fill_manual(values = force_cat_cols) +
      facet_wrap(~comparison, ncol = 2) +
      labs(
        title = "|\u0394PSI| by Force-Prediction Status",
        subtitle = "Do force-predicted genes show larger splicing changes?",
        x = NULL, y = "|\u0394PSI|"
      ) +
      theme_taxol() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 25, hjust = 1, size = rel(0.75)),
            strip.text = element_text(size = 10, face = "bold"))

    ggsave(file.path(PLOTS_DIR, "27_force_dpsi_all.pdf"), p_force_dpsi,
           width = 10, height = 10, device = SAVE_DEVICE)
    cat("  Saved: 27_force_dpsi_all.pdf\n")
  }
} else {
  cat("  Force gene files not found — skipping\n")
}

# ============================================================================
# ============================================================================
#  SECTION 15: NUCLEOPLASMIC AGITATION — FEATURE COMPARISON
# ============================================================================
# ============================================================================
#
# RATIONALE: If taxol alters splicing through nucleoplasmic agitation,
#   we expect taxol-affected exons to differ from differentiation-only
#   exons in structural features. Weak splice sites, shorter exons, or
#   specific sequence features may make splicing more sensitive to
#   mechanical perturbation of the spliceosome/speckle environment.
#
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 15: NUCLEOPLASMIC AGITATION FEATURES\n")
cat(strrep("=", 70), "\n\n")

event_info_file_s15 <- file.path(getwd(), "EVENT_INFO-mm10.tab")

if (file.exists(event_info_file_s15)) {
  events_info_full <- read.table(event_info_file_s15, header = TRUE, sep = "\t")

  # Define event classes:
  # (A) Taxol-only: significant in taxol, NOT in differentiation
  # (B) Differentiation-only: significant in differentiation, NOT in taxol
  # (C) Shared: significant in both
  # (D) Background: not significant in either

  taxol_sig_events <- unique(c(get_sig_events("Myoblast_Taxol_vs_DMSO"),
                                get_sig_events("Myotube_Taxol_vs_DMSO")))
  diff_sig_events  <- unique(c(get_sig_events("Differentiation_DMSO"),
                                get_sig_events("Differentiation_Taxol")))

  all_tested <- unique(all_events$PSI$EVENT)

  event_class_df <- data.frame(
    EVENT = all_tested,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      in_taxol = EVENT %in% taxol_sig_events,
      in_diff  = EVENT %in% diff_sig_events,
      event_class = case_when(
        in_taxol & !in_diff ~ "Taxol-Specific",
        !in_taxol & in_diff ~ "Differentiation-Specific",
        in_taxol & in_diff  ~ "Shared (Taxol + Diff)",
        TRUE                ~ "Background"
      ),
      event_class = factor(event_class,
                           levels = c("Background", "Differentiation-Specific",
                                      "Taxol-Specific", "Shared (Taxol + Diff)"))
    )

  cat(sprintf("  Taxol-specific events:           %d\n",
              sum(event_class_df$event_class == "Taxol-Specific")))
  cat(sprintf("  Differentiation-specific events: %d\n",
              sum(event_class_df$event_class == "Differentiation-Specific")))
  cat(sprintf("  Shared events:                   %d\n",
              sum(event_class_df$event_class == "Shared (Taxol + Diff)")))
  cat(sprintf("  Background events:               %d\n",
              sum(event_class_df$event_class == "Background")))

  # Merge with event info (exon sequences, lengths)
  feature_df <- event_class_df %>%
    left_join(events_info_full %>%
                dplyr::select(EVENT, Seq_A, Seq_C1, Seq_C2, LE_o),
              by = "EVENT") %>%
    filter(!is.na(Seq_A), Seq_A != "") %>%
    mutate(
      exon_length = nchar(Seq_A),
      GC_content  = sapply(Seq_A, function(s) {
        s <- toupper(s)
        gc <- sum(strsplit(s, "")[[1]] %in% c("G", "C"))
        100 * gc / nchar(s)
      })
    )

  class_colors <- c(
    "Background"                = "grey80",
    "Differentiation-Specific"  = "#2166AC",
    "Taxol-Specific"            = "#D62839",
    "Shared (Taxol + Diff)"     = "#762A83"
  )

  # 15.1 Exon length by event class
  cat("--- 15.1 Exon length by event class ---\n")

  p_len_class <- ggplot(feature_df %>% filter(event_class != "Background"),
                        aes(x = event_class, y = log10(exon_length + 1),
                            fill = event_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.2, size = 0.8) +
    stat_compare_means(
      comparisons = list(
        c("Taxol-Specific", "Differentiation-Specific"),
        c("Taxol-Specific", "Shared (Taxol + Diff)")
      ),
      method = "wilcox.test", label = "p.signif",
      hide.ns = FALSE, tip.length = 0.02, size = 3.5
    ) +
    scale_fill_manual(values = class_colors) +
    labs(
      title = "Exon Length by Splicing Response Class",
      subtitle = "Are taxol-sensitive exons structurally different?",
      x = NULL, y = expression(log[10](Exon~Length~+~1))
    ) +
    theme_taxol() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1, size = rel(0.85)))

  ggsave(file.path(PLOTS_DIR, "28_exon_length_by_class.pdf"), p_len_class,
         width = 6, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 28_exon_length_by_class.pdf\n")

  # 15.2 GC content by event class
  cat("--- 15.2 GC content by event class ---\n")

  p_gc_class <- ggplot(feature_df %>% filter(event_class != "Background"),
                       aes(x = event_class, y = GC_content, fill = event_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.2, size = 0.8) +
    stat_compare_means(
      comparisons = list(
        c("Taxol-Specific", "Differentiation-Specific"),
        c("Taxol-Specific", "Shared (Taxol + Diff)")
      ),
      method = "wilcox.test", label = "p.signif",
      hide.ns = FALSE, tip.length = 0.02, size = 3.5
    ) +
    scale_fill_manual(values = class_colors) +
    labs(
      title = "GC Content by Splicing Response Class",
      x = NULL, y = "GC Content (%)"
    ) +
    theme_taxol() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1, size = rel(0.85)))

  ggsave(file.path(PLOTS_DIR, "29_gc_by_class.pdf"), p_gc_class,
         width = 6, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 29_gc_by_class.pdf\n")

  # 15.3 Event type proportions across classes
  cat("--- 15.3 Event type proportions by class ---\n")

  feature_df_typed <- feature_df %>%
    filter(event_class != "Background") %>%
    mutate(
      event_type = case_when(
        grepl("^MmuEX", EVENT)   ~ "Exon (EX)",
        grepl("^MmuINT", EVENT)  ~ "Intron Retention (IR)",
        grepl("^MmuALT", EVENT)  ~ "Alt Splice Site",
        TRUE ~ "Other"
      )
    )

  type_class_summary <- feature_df_typed %>%
    dplyr::count(event_class, event_type) %>%
    group_by(event_class) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  type_class_colors <- c(
    "Exon (EX)"               = "#4BA3C3",
    "Intron Retention (IR)"   = "#D62839",
    "Alt Splice Site"         = "#762A83",
    "Other"                   = "grey60"
  )

  p_type_class <- ggplot(type_class_summary,
                         aes(x = event_class, y = pct, fill = event_type)) +
    geom_col(position = "stack", width = 0.65,
             color = "grey30", linewidth = 0.3) +
    scale_fill_manual(values = type_class_colors) +
    labs(
      title = "Event Type Composition by Response Class",
      subtitle = "Is intron retention enriched in taxol-specific events?",
      x = NULL, y = "% of Events", fill = "Event Type"
    ) +
    theme_taxol() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = rel(0.85)),
          legend.position = "right") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105))

  ggsave(file.path(PLOTS_DIR, "30_event_type_by_class.pdf"), p_type_class,
         width = 7, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 30_event_type_by_class.pdf\n")

  # 15.4 Summary feature table
  feature_summary <- feature_df %>%
    filter(event_class != "Background") %>%
    group_by(event_class) %>%
    summarise(
      n_events       = n(),
      median_length  = median(exon_length, na.rm = TRUE),
      mean_length    = round(mean(exon_length, na.rm = TRUE), 1),
      median_GC      = round(median(GC_content, na.rm = TRUE), 1),
      mean_GC        = round(mean(GC_content, na.rm = TRUE), 1),
      pct_microexon  = round(100 * mean(exon_length <= 27, na.rm = TRUE), 1),
      pct_short      = round(100 * mean(exon_length <= 50, na.rm = TRUE), 1),
      .groups = "drop"
    )

  write.csv(feature_summary,
            file.path(RESULTS_DIR, "event_class_features.csv"),
            row.names = FALSE)
  cat("  Feature summary saved\n")
  print(feature_summary)

} else {
  cat("  EVENT_INFO-mm10.tab not found — skipping feature analysis\n")
}

# ============================================================================
# ============================================================================
#  SECTION 16: TAXOL-SPECIFIC DEEP CHARACTERIZATION
# ============================================================================
# ============================================================================
#
# Deep-dive into the ~200 taxol-specific splicing events:
#   - Which genes & pathways
#   - Chromosomal distribution
#   - Nuclear speckle & force overlap statistics
#   - Focused GO/pathway enrichment
#
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 16: TAXOL-SPECIFIC DEEP CHARACTERIZATION\n")
cat(strrep("=", 70), "\n\n")

# Combine taxol-significant genes from both myoblast and myotube comparisons
taxol_sig_genes_myo  <- get_sig_genes("Myoblast_Taxol_vs_DMSO")
taxol_sig_genes_tube <- get_sig_genes("Myotube_Taxol_vs_DMSO")
all_taxol_genes <- unique(c(taxol_sig_genes_myo, taxol_sig_genes_tube))

cat(sprintf("  Taxol-responsive genes (myoblast): %d\n", length(taxol_sig_genes_myo)))
cat(sprintf("  Taxol-responsive genes (myotube):  %d\n", length(taxol_sig_genes_tube)))
cat(sprintf("  Union of taxol-responsive genes:   %d\n", length(all_taxol_genes)))

# 16.1 Focused enrichment: nuclear/cytoskeletal terms
cat("--- 16.1 Focused GO enrichment for taxol genes ---\n")

if (length(all_taxol_genes) >= 5) {
  entrez_taxol <- tryCatch({
    bitr(all_taxol_genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Mm.eg.db)$ENTREZID
  }, error = function(e) character(0))

  entrez_all_bg <- tryCatch({
    bitr(unique(all_events$PSI$GENE), fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Mm.eg.db)$ENTREZID
  }, error = function(e) character(0))

  if (length(entrez_taxol) >= 3) {
    # Run enrichment for BP, CC, MF
    for (ont in c("BP", "CC", "MF")) {
      ego_taxol <- tryCatch({
        enrichGO(
          gene          = entrez_taxol,
          universe      = entrez_all_bg,
          OrgDb         = org.Mm.eg.db,
          keyType       = "ENTREZID",
          ont           = ont,
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.1,
          qvalueCutoff  = 0.5
        )
      }, error = function(e) NULL)

      if (!is.null(ego_taxol) && nrow(as.data.frame(ego_taxol)) > 0) {
        df_ego <- as.data.frame(ego_taxol) %>%
          mutate(log_pval = -log10(pvalue)) %>%
          arrange(pvalue) %>%
          slice_head(n = 20) %>%
          mutate(Description = fct_reorder(Description, log_pval))

        # Highlight nuclear/speckle/cytoskeleton terms
        nuc_keywords <- c("nucle", "speckle", "chromatin", "chromosome",
                          "cytoskelet", "microtub", "actin", "myosin",
                          "envelope", "lamina", "splicing", "spliceosom",
                          "mRNA process", "RNA splicing")
        df_ego$is_nuclear <- grepl(
          paste(nuc_keywords, collapse = "|"),
          df_ego$Description, ignore.case = TRUE
        )

        p_ego <- ggplot(df_ego, aes(x = log_pval, y = Description,
                                    size = Count, color = p.adjust,
                                    shape = is_nuclear)) +
          geom_point(alpha = 0.85) +
          scale_color_gradient(low = "#2166AC", high = "#B2182B", trans = "reverse") +
          scale_size_continuous(range = c(3, 8)) +
          scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                             labels = c("Nuclear/Cytoskeletal", "Other"),
                             name = "Category") +
          labs(
            x = expression(-log[10](p-value)), y = NULL,
            color = "Adj. P", size = "Count",
            title = sprintf("GO %s: All Taxol-Responsive Genes", ont),
            subtitle = sprintf("%d genes, triangles = nuclear/cytoskeletal terms",
                               length(all_taxol_genes))
          ) +
          theme_taxol() +
          theme(axis.text.y = element_text(size = rel(0.75)),
                legend.position = "right")

        fname <- sprintf("31_taxol_GO_%s.pdf", ont)
        ggsave(file.path(PLOTS_DIR, fname), p_ego,
               width = 10, height = 7, device = SAVE_DEVICE)
        cat(sprintf("  Saved: %s\n", fname))
      }
    }
  }
}

# 16.2 Chromosomal distribution of taxol events
cat("--- 16.2 Chromosomal distribution ---\n")

all_taxol_events_combined <- unique(c(
  get_sig_events("Myoblast_Taxol_vs_DMSO"),
  get_sig_events("Myotube_Taxol_vs_DMSO")
))

if (length(all_taxol_events_combined) > 10 && exists("events_info_full")) {
  # Extract chromosome from COORD column
  taxol_coord_df <- bind_rows(
    if ("Myoblast_Taxol_vs_DMSO" %in% names(fdr_results))
      fdr_results[["Myoblast_Taxol_vs_DMSO"]] %>%
        filter(EVENT %in% all_taxol_events_combined) %>%
        dplyr::select(EVENT, GENE, COORD),
    if ("Myotube_Taxol_vs_DMSO" %in% names(fdr_results))
      fdr_results[["Myotube_Taxol_vs_DMSO"]] %>%
        filter(EVENT %in% all_taxol_events_combined) %>%
        dplyr::select(EVENT, GENE, COORD)
  ) %>%
    distinct(EVENT, .keep_all = TRUE) %>%
    mutate(
      chr = str_extract(COORD, "chr[0-9XYM]+"),
      chr = factor(chr, levels = paste0("chr", c(1:19, "X", "Y", "M")))
    ) %>%
    filter(!is.na(chr))

  # Compare to background
  bg_coord <- all_events$PSI %>%
    left_join(events_info_full %>% dplyr::select(EVENT, COORD_o), by = "EVENT") %>%
    mutate(chr = str_extract(COORD_o, "chr[0-9XYM]+")) %>%
    filter(!is.na(chr))

  chr_taxol <- taxol_coord_df %>% dplyr::count(chr, name = "taxol_n")
  chr_bg    <- bg_coord %>% dplyr::count(chr, name = "bg_n")

  chr_comparison <- full_join(chr_taxol, chr_bg, by = "chr") %>%
    replace_na(list(taxol_n = 0, bg_n = 0)) %>%
    mutate(
      taxol_pct = 100 * taxol_n / sum(taxol_n),
      bg_pct    = 100 * bg_n    / sum(bg_n),
      log2_ratio = log2((taxol_pct + 0.1) / (bg_pct + 0.1))
    )

  p_chr <- ggplot(chr_comparison, aes(x = chr, y = log2_ratio)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
    geom_col(aes(fill = log2_ratio > 0), width = 0.7,
             color = "grey30", linewidth = 0.3) +
    scale_fill_manual(values = c(`TRUE` = "#D62839", `FALSE` = "#4BA3C3"),
                      guide = "none") +
    labs(
      title = "Chromosomal Enrichment of Taxol-Responsive Events",
      subtitle = "log2(% taxol / % background)",
      x = NULL, y = expression(log[2](Enrichment))
    ) +
    theme_taxol() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.8)))

  ggsave(file.path(PLOTS_DIR, "32_chromosomal_distribution.pdf"), p_chr,
         width = 9, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 32_chromosomal_distribution.pdf\n")
}

# 16.3 Summary table of all taxol-responsive genes with annotations
cat("--- 16.3 Annotated taxol gene table ---\n")

if (exists("speckle_genes_all") && exists("force_all_genes")) {
  taxol_gene_table <- data.frame(
    gene = all_taxol_genes,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      in_myoblast    = gene %in% taxol_sig_genes_myo,
      in_myotube     = gene %in% taxol_sig_genes_tube,
      is_speckle     = gene %in% speckle_genes_all,
      is_force_all   = gene %in% force_all_genes,
      is_force_myo   = gene %in% force_myo_genes,
      is_splicing_factor = tolower(gene) %in% tolower(known_sfs)
    )

  write.csv(taxol_gene_table,
            file.path(RESULTS_DIR, "taxol_responsive_genes_annotated.csv"),
            row.names = FALSE)
  cat(sprintf("  Annotated table: %d genes saved\n", nrow(taxol_gene_table)))
  cat(sprintf("    - Speckle genes: %d (%.1f%%)\n",
              sum(taxol_gene_table$is_speckle),
              100 * mean(taxol_gene_table$is_speckle)))
  cat(sprintf("    - Force-All:     %d (%.1f%%)\n",
              sum(taxol_gene_table$is_force_all),
              100 * mean(taxol_gene_table$is_force_all)))
  cat(sprintf("    - Force-Myo:     %d (%.1f%%)\n",
              sum(taxol_gene_table$is_force_myo),
              100 * mean(taxol_gene_table$is_force_myo)))
  cat(sprintf("    - Known SFs:     %d (%.1f%%)\n",
              sum(taxol_gene_table$is_splicing_factor),
              100 * mean(taxol_gene_table$is_splicing_factor)))
}

# ============================================================================
# ============================================================================
#  SECTION 17: INTEGRATED NUCLEAR AGITATION MODEL
# ============================================================================
# ============================================================================
#
# SYNTHESIS: Bring together speckle, force, and feature analyses to
#   build an integrated picture of whether taxol acts through
#   nucleoplasmic agitation.
#
#   Key predictions of the nuclear agitation hypothesis:
#   1. Taxol-affected genes should be enriched in speckle-associated genes
#   2. Taxol-affected events should overlap with force-predicted genes
#   3. Taxol-specific exons may have distinct structural features
#      (weaker splice sites -> more sensitive to mechanical perturbation)
#   4. Effect should be independent of differentiation program
#
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 17: INTEGRATED NUCLEAR AGITATION MODEL\n")
cat(strrep("=", 70), "\n\n")

# 17.1 Combined annotation heatmap for taxol events
cat("--- 17.1 Annotation heatmap for taxol-significant events ---\n")

if (exists("taxol_gene_table") && nrow(taxol_gene_table) > 10) {
  # Build annotation matrix
  annot_mat <- taxol_gene_table %>%
    column_to_rownames("gene") %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

  # Only show genes with at least one annotation
  has_annot <- rowSums(annot_mat[, c("is_speckle", "is_force_all",
                                      "is_force_myo", "is_splicing_factor")]) > 0
  annot_sub <- annot_mat[has_annot, , drop = FALSE]

  if (nrow(annot_sub) > 3) {
    colnames(annot_sub) <- c("Myoblast", "Myotube", "Speckle",
                             "Force-All", "Force-Myo", "Spl. Factor")

    ht_annot <- Heatmap(
      annot_sub,
      name = "Present",
      col  = colorRamp2(c(0, 1), c("grey95", "#D62839")),
      cluster_rows  = TRUE,
      cluster_columns = FALSE,
      show_row_names = nrow(annot_sub) <= 60,
      row_names_gp   = gpar(fontsize = 7),
      column_names_rot = 45,
      column_title   = "Taxol-Responsive Genes: Multi-Layer Annotations",
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      rect_gp        = gpar(col = "grey80", lwd = 0.3),
      heatmap_legend_param = list(
        title    = "Present",
        at       = c(0, 1),
        labels   = c("No", "Yes")
      )
    )

    annot_hm_file <- file.path(PLOTS_DIR, "33_taxol_annotation_heatmap.pdf")
    pdf(annot_hm_file, width = 7,
        height = max(5, min(20, nrow(annot_sub) * 0.2)))
    draw(ht_annot)
    dev.off()
    cat("  Saved: 33_taxol_annotation_heatmap.pdf\n")
  }
}

# 17.2 Three-way overlap summary
cat("--- 17.2 Three-way overlap summary ---\n")

if (exists("speckle_genes_all") && exists("force_all_genes")) {
  threeway_df <- data.frame(
    category = c("Taxol Only",
                 "Taxol + Speckle",
                 "Taxol + Force",
                 "Taxol + Speckle + Force",
                 "Total Taxol Genes"),
    n = c(
      sum(!all_taxol_genes %in% speckle_genes_all &
            !all_taxol_genes %in% force_all_genes),
      sum(all_taxol_genes %in% speckle_genes_all &
            !all_taxol_genes %in% force_all_genes),
      sum(!all_taxol_genes %in% speckle_genes_all &
            all_taxol_genes %in% force_all_genes),
      sum(all_taxol_genes %in% speckle_genes_all &
            all_taxol_genes %in% force_all_genes),
      length(all_taxol_genes)
    ),
    stringsAsFactors = FALSE
  ) %>%
    filter(category != "Total Taxol Genes") %>%
    mutate(pct = round(100 * n / length(all_taxol_genes), 1))

  threeway_colors <- c(
    "Taxol Only"                = "grey70",
    "Taxol + Speckle"           = "#762A83",
    "Taxol + Force"             = "#E66101",
    "Taxol + Speckle + Force"   = "#D62839"
  )

  p_threeway <- ggplot(threeway_df,
                       aes(x = reorder(category, -n), y = n, fill = category)) +
    geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
              vjust = -0.3, family = FONT_FAMILY, size = 3.5) +
    scale_fill_manual(values = threeway_colors) +
    labs(
      title = "Taxol-Responsive Genes: Multi-Layer Annotation",
      subtitle = "Speckle and force-prediction overlaps",
      x = NULL, y = "Number of Genes"
    ) +
    theme_taxol() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1, size = rel(0.85))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)))

  ggsave(file.path(PLOTS_DIR, "34_threeway_overlap.pdf"), p_threeway,
         width = 7, height = 5, device = SAVE_DEVICE)
  cat("  Saved: 34_threeway_overlap.pdf\n")
}

# 17.3 dPSI magnitude comparison: speckle + force vs others
cat("--- 17.3 dPSI magnitude by annotation layers ---\n")

layer_dpsi_list <- list()
for (tc in names(comparison_labels)) {
  if (!(tc %in% names(fdr_results))) next
  if (!exists("speckle_genes_all") || !exists("force_all_genes")) next

  df_layers <- fdr_results[[tc]] %>%
    filter(!is.na(deltapsi)) %>%
    mutate(
      annotation = case_when(
        GENE %in% speckle_genes_all & GENE %in% force_all_genes ~ "Speckle + Force",
        GENE %in% speckle_genes_all ~ "Speckle Only",
        GENE %in% force_all_genes   ~ "Force Only",
        TRUE ~ "Neither"
      ),
      annotation = factor(annotation,
                          levels = c("Neither", "Speckle Only",
                                     "Force Only", "Speckle + Force")),
      comparison = comparison_labels[tc]
    )
  layer_dpsi_list[[tc]] <- df_layers
}

if (length(layer_dpsi_list) > 0) {
  layer_dpsi_all <- bind_rows(layer_dpsi_list) %>%
    mutate(comparison = factor(comparison, levels = unname(comparison_labels)))

  layer_cols <- c("Neither" = "grey75", "Speckle Only" = "#762A83",
                  "Force Only" = "#E66101", "Speckle + Force" = "#D62839")

  p_layer_dpsi <- ggplot(layer_dpsi_all,
                         aes(x = annotation, y = deltapsi, fill = annotation)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.65) +
    geom_jitter(width = 0.15, alpha = 0.08, size = 0.4) +
    stat_compare_means(ref.group = "Neither",
                       method = "wilcox.test", label = "p.signif",
                       hide.ns = TRUE, size = 3) +
    scale_fill_manual(values = layer_cols) +
    facet_wrap(~comparison, ncol = 2) +
    labs(
      title = "\u0394PSI by Annotation Layer",
      subtitle = "Nuclear agitation hypothesis: speckle + force genes most affected?",
      x = NULL, y = "\u0394PSI"
    ) +
    theme_taxol() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1, size = rel(0.75)),
          strip.text = element_text(size = 10, face = "bold"))

  ggsave(file.path(PLOTS_DIR, "35_layer_dpsi_all.pdf"), p_layer_dpsi,
         width = 10, height = 10, device = SAVE_DEVICE)
  cat("  Saved: 35_layer_dpsi_all.pdf\n")
}

# 17.4 Taxol dPSI direction: do speckle genes show biased skipping?
cat("--- 17.4 Taxol effect direction analysis ---\n")

for (tc in names(comparison_labels)) {
  if (!(tc %in% names(fdr_results))) next
  if (!exists("speckle_genes_all")) next

  sig_events_tc <- fdr_results[[tc]] %>%
    filter(!is.na(FDR), FDR <= 0.05, abs(deltapsi) >= 0.1) %>%
    mutate(
      is_speckle = GENE %in% speckle_genes_all,
      layer = ifelse(is_speckle, "Speckle Gene", "Non-Speckle Gene"),
      direction = ifelse(deltapsi > 0, "Included", "Skipped")
    )

  if (nrow(sig_events_tc) > 10) {
    direction_test <- table(sig_events_tc$layer, sig_events_tc$direction)
    if (all(dim(direction_test) == c(2, 2))) {
      ft_dir <- fisher.test(direction_test)
      cat(sprintf("  %s — Direction bias (Fisher's test): OR=%.2f, p=%.4g\n",
                  comparison_labels[tc], ft_dir$estimate, ft_dir$p.value))
    }
  }
}

# 17.5 PSI heatmap: taxol events in force-predicted genes
cat("--- 17.5 Taxol x Force PSI heatmap ---\n")

if (exists("force_all_genes")) {
  # Get taxol-significant events that are in force-predicted genes
  taxol_force_events <- all_taxol_events_combined[
    all_taxol_events_combined %in% all_events$PSI$EVENT
  ]

  force_event_genes <- all_events$PSI %>%
    filter(EVENT %in% taxol_force_events, GENE %in% force_all_genes)

  if (nrow(force_event_genes) > 3) {
    hm_force_mat <- force_event_genes %>%
      dplyr::select(EVENT, GENE, all_of(metadata$sample_id))

    # Label rows with gene name
    hm_labels <- paste0(hm_force_mat$GENE, " (", hm_force_mat$EVENT, ")")
    hm_mat_num <- hm_force_mat %>%
      dplyr::select(-EVENT, -GENE) %>%
      mutate(across(everything(), as.numeric)) %>%
      as.matrix()
    rownames(hm_mat_num) <- hm_labels
    hm_mat_num <- hm_mat_num[complete.cases(hm_mat_num), , drop = FALSE]

    if (nrow(hm_mat_num) > 2) {
      col_fun_force <- colorRamp2(c(0, 50, 100), viridis(3))

      ha_force <- HeatmapAnnotation(
        Condition = metadata$condition[match(colnames(hm_mat_num), metadata$sample_id)],
        col = list(Condition = condition_colors),
        annotation_name_gp = gpar(fontsize = 9)
      )

      ht_force <- Heatmap(
        hm_mat_num,
        name = "PSI",
        col  = col_fun_force,
        clustering_distance_rows = "pearson",
        top_annotation = ha_force,
        column_title   = "Taxol-Responsive Events in Force-Predicted Genes",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        column_labels  = metadata$condition[match(colnames(hm_mat_num),
                                                  metadata$sample_id)],
        column_names_rot = 45,
        column_names_gp  = gpar(fontsize = 9),
        show_row_names   = nrow(hm_mat_num) <= 50,
        row_names_gp     = gpar(fontsize = 7),
        rect_gp          = gpar(col = "grey90", lwd = 0.3),
        use_raster       = FALSE
      )

      hm_force_file <- file.path(PLOTS_DIR, "36_taxol_force_heatmap.pdf")
      pdf(hm_force_file, width = 10,
          height = max(6, min(18, nrow(hm_mat_num) * 0.25)))
      draw(ht_force, heatmap_legend_side = "right")
      dev.off()
      cat("  Saved: 36_taxol_force_heatmap.pdf\n")
    }
  }
}

# 17.6 Hypothesis Summary Statistics
cat("\n--- Nuclear Agitation Hypothesis Summary ---\n")
cat("  Evidence assessment:\n")

if (exists("speckle_enrich_df")) {
  taxol_speckle <- speckle_enrich_df %>%
    filter(grepl("Taxol", comparison))
  for (i in seq_len(nrow(taxol_speckle))) {
    cat(sprintf("  [Speckle] %s: OR=%.2f, p=%.4g %s\n",
                taxol_speckle$comparison[i],
                taxol_speckle$odds_ratio[i],
                taxol_speckle$p_value[i],
                ifelse(taxol_speckle$p_value[i] < 0.05, "\u2713", "\u2717")))
  }
}

if (exists("overlap_df")) {
  taxol_force_ov <- overlap_df %>%
    filter(grepl("Taxol", comparison))
  for (i in seq_len(nrow(taxol_force_ov))) {
    cat(sprintf("  [Force]   %s x %s: fold=%.2f, p=%.4g %s\n",
                taxol_force_ov$comparison[i],
                taxol_force_ov$force_set[i],
                taxol_force_ov$fold_enrichment[i],
                taxol_force_ov$p_value[i],
                ifelse(taxol_force_ov$p_value[i] < 0.05, "\u2713", "\u2717")))
  }
}

cat("\n")

# ============================================================================
# FINAL OUTPUT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Output directories:\n")
cat(sprintf("  Results: %s\n", RESULTS_DIR))
cat(sprintf("  Plots:   %s\n", PLOTS_DIR))
cat(sprintf("  Plots generated: %d\n", length(list.files(PLOTS_DIR, pattern = "\\.pdf$"))))
cat(sprintf("  Result files: %d\n", length(list.files(RESULTS_DIR, pattern = "\\.csv$"))))

cat("\nSession info:\n")
sessionInfo()
