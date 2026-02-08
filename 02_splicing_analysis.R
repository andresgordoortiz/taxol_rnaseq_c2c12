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

font_add_google(name = "Courier Prime", family = "Courier Prime")
showtext_auto()
showtext_opts(dpi = 300)
FONT_FAMILY <- "Courier Prime"

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
INCLUSION_TABLE <- file.path(getwd(), "inclusion_tables",
                             "INCLUSION_LEVELS_FULL-mm10.tab")

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
result_files <- result_files[!grepl("summary", result_files)]

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

ggsave(file.path(PLOTS_DIR, "01_pca_psi.pdf"), p_pca, width = 7, height = 5)
cat("  Saved: 01_pca_psi.pdf\n")

# ----------------------------------------------------------------------------
# 1.2 PSI Distribution (U-shape check for exons)
# ----------------------------------------------------------------------------

cat("--- 1.2 PSI distributions ---\n")

# Flatten PSI values for each condition
psi_long <- exons$PSI %>%
  select(EVENT, all_of(metadata$sample_id)) %>%
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
       width = 8, height = 6)
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
    ggsave(file.path(PLOTS_DIR, fname), p, width = 5.5, height = 4)
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
    mutate(direction = ifelse(deltapsi > 0, "Included/Retained", "Skipped/Spliced Out"))
  data.frame(
    comparison = name,
    direction  = sig$direction,
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
  labs(
    title = "Differential Splicing Events per Comparison",
    subtitle = "FDR \u2264 0.05, |\u0394PSI| \u2265 0.1",
    x = NULL, y = "Number of Events", fill = NULL
  ) +
  theme_taxol() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = rel(0.8)),
    legend.position = "top"
  )

ggsave(file.path(PLOTS_DIR, "04_de_events_summary.pdf"), p_summary,
       width = 10, height = 5)
cat("  Saved: 04_de_events_summary.pdf\n")

# ----------------------------------------------------------------------------
# 3.2 Taxol effect comparison: scatter of dPSI in myoblasts vs myotubes
# ----------------------------------------------------------------------------

cat("--- 3.2 Taxol effect scatter (myoblasts vs myotubes) ---\n")

if (all(c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO") %in%
        names(fdr_results))) {

  taxol_myo <- fdr_results[["Myoblast_Taxol_vs_DMSO"]] %>%
    select(EVENT, dpsi_myoblast = deltapsi, FDR_myoblast = FDR)

  taxol_tube <- fdr_results[["Myotube_Taxol_vs_DMSO"]] %>%
    select(EVENT, dpsi_myotube = deltapsi, FDR_myotube = FDR)

  scatter_data <- inner_join(taxol_myo, taxol_tube, by = "EVENT") %>%
    mutate(
      sig_myoblast = !is.na(FDR_myoblast) & FDR_myoblast <= 0.05 & abs(dpsi_myoblast) >= 0.1,
      sig_myotube  = !is.na(FDR_myotube)  & FDR_myotube  <= 0.05 & abs(dpsi_myotube)  >= 0.1,
      category = case_when(
        sig_myoblast & sig_myotube  ~ "Both",
        sig_myoblast & !sig_myotube ~ "Myoblast only",
        !sig_myoblast & sig_myotube ~ "Myotube only",
        TRUE ~ "NS"
      )
    )

  cat_colors_scatter <- c("Both" = "#762A83", "Myoblast only" = "#D62839",
                          "Myotube only" = "#2166AC", "NS" = "grey85")

  p_taxol_scatter <- ggplot(scatter_data,
                            aes(x = dpsi_myoblast, y = dpsi_myotube, color = category)) +
    geom_point(data = filter(scatter_data, category == "NS"), alpha = 0.1, size = 0.5) +
    geom_point(data = filter(scatter_data, category != "NS"), alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted",
                linewidth = 0.6) +
    scale_color_manual(values = cat_colors_scatter,
                       breaks = c("Both", "Myoblast only", "Myotube only")) +
    labs(
      title = "Taxol Effect on Exon Splicing",
      subtitle = "Comparing Taxol-induced \u0394PSI in myoblasts vs myotubes",
      x = "\u0394PSI (Taxol vs DMSO in Myoblasts)",
      y = "\u0394PSI (Taxol vs DMSO in Myotubes)",
      color = "Significant in:"
    ) +
    coord_fixed(ratio = 1, xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_taxol() +
    theme(legend.position = "right")

  ggsave(file.path(PLOTS_DIR, "05_taxol_scatter.pdf"), p_taxol_scatter,
         width = 6, height = 6)
  cat("  Saved: 05_taxol_scatter.pdf\n")

  # Blocking-style summary barplot
  scatter_summary <- scatter_data %>%
    filter(category != "NS") %>%
    count(category) %>%
    mutate(pct = 100 * n / sum(n))

  p_taxol_bar <- ggplot(scatter_summary, aes(x = category, y = n, fill = category)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)), vjust = -0.3,
              family = FONT_FAMILY, size = 3.5) +
    scale_fill_manual(values = cat_colors_scatter) +
    labs(
      title = "Distribution of Taxol-Responsive Events",
      x = NULL, y = "Number of Events"
    ) +
    theme_taxol() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  ggsave(file.path(PLOTS_DIR, "06_taxol_category_bar.pdf"), p_taxol_bar,
         width = 5, height = 4)
  cat("  Saved: 06_taxol_category_bar.pdf\n")
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
  pdf(file.path(PLOTS_DIR, "07_upset_events.pdf"), width = 8, height = 5)
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
  pdf(file.path(PLOTS_DIR, "08_upset_genes.pdf"), width = 8, height = 5)
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

  # 5.1 dPSI histogram
  cat("--- 5.1 dPSI distribution ---\n")

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
         width = 8, height = 5)
  cat("  Saved: 10_dpsi_distribution.pdf\n")

  # 5.2 Skewness by exon count per gene (bootstrapped)
  cat("--- 5.2 Bootstrapped skewness by gene complexity ---\n")

  # Get exon counts from VastDB EVENT_INFO if available
  event_info_file <- file.path(getwd(), "EVENT_INFO-mm10.tab")
  if (file.exists(event_info_file)) {
    vastdb_events <- read.table(event_info_file, header = TRUE, sep = "\t")
    vastdb_exons <- vastdb_events[grep("EX", vastdb_events$EVENT), ]
    gene_exon_counts <- table(vastdb_exons$GENE)

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
      mutate(
        exon_bin = cut(numb_exons,
                       breaks = c(0, 10, 50, Inf),
                       labels = c("1\u201310", "11\u201350", "51+"),
                       right = TRUE)
      ) %>%
      filter(!is.na(exon_bin))

    set.seed(123)
    boot_results <- binned %>%
      group_by(exon_bin) %>%
      summarise(
        boot_skews = list(replicate(1000, skewness_manual(sample(deltapsi, replace = TRUE)))),
        .groups = "drop"
      ) %>%
      mutate(
        skew_estimate = map_dbl(boot_skews, ~ mean(.x, na.rm = TRUE)),
        ci_lower      = map_dbl(boot_skews, ~ quantile(.x, 0.025, na.rm = TRUE)),
        ci_upper      = map_dbl(boot_skews, ~ quantile(.x, 0.975, na.rm = TRUE)),
        significant   = (ci_lower > 0) | (ci_upper < 0)
      )

    boot_long <- boot_results %>%
      mutate(exon_bin = factor(exon_bin, levels = c("1\u201310", "11\u201350", "51+"))) %>%
      unnest(boot_skews)

    p_skew <- ggplot(boot_results, aes(x = exon_bin, y = skew_estimate)) +
      geom_jitter(data = boot_long,
                  aes(fill = exon_bin, y = boot_skews),
                  position = position_jitter(width = 0.15),
                  shape = 21, color = "gray30", size = 2,
                  alpha = 0.05, stroke = 0.8) +
      geom_smooth(data = boot_long,
                  aes(x = as.numeric(exon_bin), y = boot_skews),
                  method = "lm", formula = y ~ poly(x, 2),
                  se = FALSE, linetype = "dashed", color = "grey40",
                  linewidth = 1.2) +
      geom_point(aes(fill = exon_bin, color = significant, shape = significant),
                 size = 4, alpha = 0.9, stroke = 1) +
      scale_fill_manual(values = c(
        "1\u201310"  = "#d1cbe5",
        "11\u201350" = "#9671bd",
        "51+"        = "#6a408d"
      ), name = "Exon Bin") +
      scale_color_manual(values = c(`TRUE` = "#378d94", `FALSE` = "grey60"),
                         name = "Significant") +
      scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                         name = "Significant") +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
      scale_y_continuous(limits = c(-2.5, 2.5)) +
      labs(
        title = "Bootstrapped Skewness of Taxol Splicing",
        subtitle = "Does gene complexity affect inclusion/skipping bias?",
        x = "Number of exons per gene",
        y = "Skewness of \u0394PSI Distribution"
      ) +
      theme_taxol() +
      theme(aspect.ratio = 1, legend.position = "top", legend.direction = "horizontal")

    ggsave(file.path(PLOTS_DIR, "11_skewness_by_exon_count.pdf"), p_skew,
           width = 6, height = 5)
    cat("  Saved: 11_skewness_by_exon_count.pdf\n")
  } else {
    cat("  Skipping skewness analysis: EVENT_INFO-mm10.tab not found\n")
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
if (file.exists(event_info_file) && primary_comp %in% names(fdr_results)) {

  events_info <- read.table(event_info_file, header = TRUE, sep = "\t")

  diff_exons <- fdr_results[[primary_comp]] %>%
    filter(!is.na(FDR), FDR <= 0.05, abs(deltapsi) >= 0.1)

  diff_exons <- diff_exons %>%
    left_join(events_info %>% select(EVENT, Seq_A, LE_o), by = "EVENT") %>%
    mutate(
      LENGTH = nchar(Seq_A),
      direction = ifelse(deltapsi > 0, "Included", "Skipped")
    )

  # 6.1 Exon length by direction
  cat("--- 6.1 Exon length by inclusion direction ---\n")

  diff_exons$length_bin <- cut(
    diff_exons$LENGTH,
    breaks = c(-Inf, 27, 50, Inf),
    labels = c("\u226427", "28\u201350", ">50")
  )

  seq_cols_len <- c("\u226427" = "#c6dbef", "28\u201350" = "#6baed6", ">50" = "#2171b5")

  p_len <- ggplot(diff_exons, aes(x = length_bin, y = deltapsi, fill = length_bin)) +
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
    labs(
      title = "\u0394PSI by Exon Length Category",
      subtitle = "Taxol vs DMSO in Myoblasts",
      x = "Length bin", y = "\u0394PSI"
    ) +
    theme_taxol() +
    theme(aspect.ratio = 1, legend.position = "none")

  ggsave(file.path(PLOTS_DIR, "12_exon_length_dpsi.pdf"), p_len,
         width = 5, height = 5)
  cat("  Saved: 12_exon_length_dpsi.pdf\n")

  # 6.2 GC content by direction
  cat("--- 6.2 GC content by splicing direction ---\n")

  gc_content <- function(seqs) {
    sapply(seqs, function(seq) {
      seq <- toupper(seq)
      gc <- sum(strsplit(seq, "")[[1]] %in% c("G", "C"))
      total <- nchar(seq)
      if (total == 0) return(NA)
      return(100 * gc / total)
    })
  }

  diff_exons$GC <- gc_content(diff_exons$Seq_A)

  p_gc <- ggplot(diff_exons, aes(x = direction, y = GC, fill = direction)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 1.2) +
    stat_compare_means(method = "wilcox.test", label.y = max(diff_exons$GC, na.rm = TRUE) + 3,
                       size = 4, family = FONT_FAMILY) +
    scale_fill_manual(values = c("Included" = "#D62839", "Skipped" = "#4BA3C3")) +
    labs(
      title = "GC Content of Differentially Spliced Exons",
      x = "Splicing Direction", y = "GC Content (%)"
    ) +
    theme_taxol() +
    theme(legend.position = "none")

  ggsave(file.path(PLOTS_DIR, "13_gc_content_by_direction.pdf"), p_gc,
         width = 5, height = 5)
  cat("  Saved: 13_gc_content_by_direction.pdf\n")
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
  # Build PSI matrix
  all_filtered <- filterEvents(all_events, types = c("C1","C2","C3","S","MIC","IR"), N = 10)
  heatmap_mat <- all_filtered$PSI[, c(1, 7:ncol(all_filtered$PSI))] %>%
    filter(EVENT %in% all_sig_events)

  # Use metadata sample IDs as columns
  heatmap_mat_num <- heatmap_mat %>%
    select(EVENT, all_of(metadata$sample_id)) %>%
    column_to_rownames("EVENT") %>%
    as.matrix()

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
      annotation_name_gp = gpar(fontsize = 10, fontfamily = FONT_FAMILY)
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
      column_title_gp  = gpar(fontsize = 16, fontface = "bold", fontfamily = FONT_FAMILY),
      row_title        = "Splicing Events",
      row_title_gp     = gpar(fontsize = 12, fontface = "italic", fontfamily = FONT_FAMILY),
      column_labels    = metadata$condition[match(colnames(heatmap_mat_num), metadata$sample_id)],
      column_names_rot = 45,
      column_names_gp  = gpar(fontsize = 10, fontfamily = FONT_FAMILY),
      row_km           = n_clust,
      show_row_names   = FALSE,
      rect_gp          = gpar(col = "grey90", lwd = 0.3),
      row_gap          = unit(3, "mm"),
      column_gap       = unit(2, "mm"),
      heatmap_legend_param = list(
        title      = "PSI",
        title_gp   = gpar(fontsize = 12, fontface = "bold", fontfamily = FONT_FAMILY),
        labels_gp  = gpar(fontsize = 9, fontfamily = FONT_FAMILY)
      )
    )

    pdf(file.path(PLOTS_DIR, "14_heatmap_splicing_dynamics.pdf"),
        width = 10, height = 12)
    draw(ht, heatmap_legend_side = "right")
    dev.off()
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
        ggsave(file.path(PLOTS_DIR, fname), p, width = 9, height = 6)
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
    select(EVENT, GENE, dpsi_dmso = deltapsi, FDR_dmso = FDR)
  diff_taxol <- fdr_results[["Differentiation_Taxol"]] %>%
    select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)

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
      x = "\u0394PSI Differentiation (DMSO)",
      y = "\u0394PSI Differentiation (Taxol)",
      color = NULL
    ) +
    coord_fixed(ratio = 1, xlim = c(-1, 1), ylim = c(-1, 1)) +
    theme_taxol() +
    theme(legend.position = "right")

  ggsave(file.path(PLOTS_DIR, "16_interaction_scatter.pdf"), p_interaction,
         width = 7, height = 6)
  cat("  Saved: 16_interaction_scatter.pdf\n")

  # Stacked bar of categories
  interaction_summary <- interaction %>%
    filter(category != "NS") %>%
    count(category) %>%
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
         width = 5, height = 5)
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
    count(comparison, event_subtype) %>%
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
         width = 9, height = 5)
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
          count(SF_Gene, name = "n_targets") %>%
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

        ggsave(file.path(PLOTS_DIR, "19_sf_interactions.pdf"), p_sf,
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
         width = 10, height = 10)
  cat("  Saved: 20_combined_volcano.pdf\n")
}

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
