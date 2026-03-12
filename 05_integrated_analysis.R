# ============================================================================
# 05_integrated_analysis.R
# Integrated Splicing + Gene Expression Analysis
# Taxol Effect in C2C12 Myoblasts and Myotubes
#
# PURPOSE: Integrate results from:
#   - 01_fdr_calculation.R & 02_splicing_analysis.R  (splicing: betAS/vast-tools)
#   - 03_gene_expression.R   (gene expression: DESeq2 with LFC shrinkage)
#   - 04_matt_feature_analysis.R (exon features: Matt)
#
# KEY QUESTION: Are genes with differential splicing also differentially
#   expressed, or do splicing and expression changes affect distinct genes?
#   This dissection is critical for understanding the taxol mechanism:
#   if taxol primarily affects splicing through nucleoplasmic agitation,
#   we expect splicing changes largely independent of expression changes.
#
# EXPERIMENTAL DESIGN:
#   12 samples — 4 conditions × 3 replicates
#   Comparisons: Myoblast Taxol vs DMSO, Myotube Taxol vs DMSO,
#                Differentiation DMSO, Differentiation Taxol
#
# OUTPUT:
#   plots/int_*                              — integration figures
#   results/integrated_splicing_expression_* — merged tables
#
# Author: Andrés Gordo Ortiz
# ============================================================================

cat(strrep("=", 70), "\n")
cat("05_integrated_analysis.R — Splicing × Expression Integration\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
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
  library(UpSetR)
  library(grid)
  library(viridis)
  library(scales)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(ggpubr)
})

# ============================================================================
# FONT SETUP
# ============================================================================

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
  cat("  Font backend: base pdf (Courier)\n")
}

open_pdf <- function(file, ...) {
  if (HAS_CAIRO) cairo_pdf(file, ...) else pdf(file, ...)
}

# ============================================================================
# SHARED THEME
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

condition_colors <- c(
  "Myoblast_DMSO"  = "#4BA3C3",
  "Myoblast_Taxol" = "#D62839",
  "Myotube_DMSO"   = "#2166AC",
  "Myotube_Taxol"  = "#B2182B"
)

# Integration categories
int_colors <- c(
  "Both"             = "#762A83",
  "Splicing only"    = "#D62839",
  "Expression only"  = "#4BA3C3",
  "Neither"          = "grey85"
)

# ============================================================================
# PATHS
# ============================================================================

RESULTS_DIR   <- file.path(getwd(), "results")
PLOTS_DIR     <- file.path(getwd(), "plots")
METADATA_FILE <- file.path(getwd(), "metadata", "metadata.csv")

if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("--- Loading splicing FDR results ---\n")

comparison_names <- c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO",
                       "Differentiation_DMSO", "Differentiation_Taxol")

comparison_labels <- c(
  "Myoblast_Taxol_vs_DMSO" = "Myoblast: Taxol vs DMSO",
  "Myotube_Taxol_vs_DMSO"  = "Myotube: Taxol vs DMSO",
  "Differentiation_DMSO"   = "Differentiation (DMSO)",
  "Differentiation_Taxol"  = "Differentiation (Taxol)"
)

# Load splicing results
splicing_results <- list()
for (comp_name in comparison_names) {
  f <- file.path(RESULTS_DIR, paste0(comp_name, ".csv"))
  if (file.exists(f)) {
    splicing_results[[comp_name]] <- read.csv(f)
    cat(sprintf("  Splicing: %s (%d events)\n", comp_name,
                nrow(splicing_results[[comp_name]])))
  }
}

# Load gene expression results
cat("\n--- Loading DESeq2 results ---\n")
deseq_results <- list()
for (comp_name in comparison_names) {
  f <- file.path(RESULTS_DIR, paste0("deseq2_", comp_name, ".csv"))
  if (file.exists(f)) {
    deseq_results[[comp_name]] <- read.csv(f)
    cat(sprintf("  Expression: %s (%d genes)\n", comp_name,
                nrow(deseq_results[[comp_name]])))
  }
}

# Load VST matrix (for heatmaps)
vst_file <- file.path(RESULTS_DIR, "deseq2_vst_matrix.csv")
if (file.exists(vst_file)) {
  vst_df <- read.csv(vst_file)
  cat(sprintf("  VST matrix: %d genes × %d samples\n",
              nrow(vst_df), ncol(vst_df) - 2))  # minus ensembl_id and gene_symbol
}

# Metadata
metadata <- read.csv(METADATA_FILE, stringsAsFactors = FALSE) %>%
  arrange(condition)

# Check that both result types are available
if (length(splicing_results) == 0) stop("No splicing results found in ", RESULTS_DIR)
if (length(deseq_results)   == 0) stop("No DESeq2 results found in ", RESULTS_DIR,
                                        "\nRun 03_gene_expression.R first.")

# ============================================================================
# GENE SYMBOL MAPPING
# ============================================================================

# Build ensembl → symbol map from DESeq2 results
gene_map <- deseq_results[[1]] %>%
  dplyr::select(ensembl_id, gene_symbol) %>%
  distinct() %>%
  filter(!is.na(gene_symbol))

# ============================================================================
# THRESHOLDS
# ============================================================================

# Splicing significance (betAS FDR)
SPLICING_FDR_THRESH <- 0.05
SPLICING_DPSI_THRESH <- 0.1

# Gene expression significance (DESeq2)
EXPR_PADJ_THRESH <- 0.05
EXPR_LFC_THRESH  <- 1.0

cat(sprintf("\n  Splicing:   FDR ≤ %.2f, |dPSI| ≥ %.2f\n",
            SPLICING_FDR_THRESH, SPLICING_DPSI_THRESH))
cat(sprintf("  Expression: padj ≤ %.2f, |LFC| ≥ %.1f\n",
            EXPR_PADJ_THRESH, EXPR_LFC_THRESH))

# ============================================================================
# ============================================================================
#  SECTION 1: GENE-LEVEL INTEGRATION (per comparison)
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 1: GENE-LEVEL INTEGRATION\n")
cat(strrep("=", 70), "\n\n")

# Helper functions
get_splicing_genes <- function(comp_name) {
  if (comp_name %in% names(splicing_results)) {
    df <- splicing_results[[comp_name]]
    unique(na.omit(df$GENE[!is.na(df$FDR) & df$FDR <= SPLICING_FDR_THRESH &
                             !is.na(df$deltapsi) & abs(df$deltapsi) >= SPLICING_DPSI_THRESH]))
  } else character(0)
}

get_expression_genes <- function(comp_name) {
  if (comp_name %in% names(deseq_results)) {
    df <- deseq_results[[comp_name]]
    unique(na.omit(df$gene_symbol[!is.na(df$significance) &
                                    df$significance != "Not Significant"]))
  } else character(0)
}

# Build integration table per comparison
integration_summary <- list()

for (comp_name in comparison_names) {
  spl_genes  <- get_splicing_genes(comp_name)
  expr_genes <- get_expression_genes(comp_name)

  # All genes in the universe
  all_genes <- unique(c(spl_genes, expr_genes))

  if (length(all_genes) == 0) next

  int_df <- data.frame(
    gene = all_genes,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      splicing_sig   = gene %in% spl_genes,
      expression_sig = gene %in% expr_genes,
      category = case_when(
        splicing_sig & expression_sig  ~ "Both",
        splicing_sig & !expression_sig ~ "Splicing only",
        !splicing_sig & expression_sig ~ "Expression only",
        TRUE                           ~ "Neither"
      ),
      comparison = comp_name
    )

  integration_summary[[comp_name]] <- int_df

  n_both <- sum(int_df$category == "Both")
  n_spl  <- sum(int_df$category == "Splicing only")
  n_expr <- sum(int_df$category == "Expression only")

  cat(sprintf("  %s:\n", comp_name))
  cat(sprintf("    Splicing significant:   %d genes\n", length(spl_genes)))
  cat(sprintf("    Expression significant: %d genes\n", length(expr_genes)))
  cat(sprintf("    Both:                   %d genes\n", n_both))
  cat(sprintf("    Splicing only:          %d genes\n", n_spl))
  cat(sprintf("    Expression only:        %d genes\n", n_expr))

  # Fisher's exact test for enrichment/depletion of overlap
  if (comp_name %in% names(splicing_results) &&
      comp_name %in% names(deseq_results)) {

    # Background: all genes present in both analyses
    bg_spl  <- unique(na.omit(splicing_results[[comp_name]]$GENE))
    bg_expr <- unique(na.omit(deseq_results[[comp_name]]$gene_symbol))
    bg_all  <- intersect(bg_spl, bg_expr)

    n_bg     <- length(bg_all)
    n_s      <- sum(bg_all %in% spl_genes)
    n_e      <- sum(bg_all %in% expr_genes)
    n_se     <- sum(bg_all %in% spl_genes & bg_all %in% expr_genes)

    cont_mat <- matrix(c(n_se,
                          n_s - n_se,
                          n_e - n_se,
                          n_bg - n_s - n_e + n_se),
                        nrow = 2)
    ft <- fisher.test(cont_mat, alternative = "greater")

    cat(sprintf("    Fisher's test (overlap enrichment): OR = %.2f, p = %.2e\n",
                ft$estimate, ft$p.value))
  }
  cat("\n")
}

# Save integration summary
int_combined <- bind_rows(integration_summary)
write.csv(int_combined, file.path(RESULTS_DIR, "integrated_splicing_expression_genes.csv"),
          row.names = FALSE)
cat("  Saved: integrated_splicing_expression_genes.csv\n")

# ============================================================================
# ============================================================================
#  SECTION 2: INTEGRATION SCATTER — dPSI vs LFC
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 2: dPSI vs LFC SCATTER\n")
cat(strrep("=", 70), "\n\n")

scatter_plots <- list()

for (comp_name in comparison_names) {
  if (!(comp_name %in% names(splicing_results)) ||
      !(comp_name %in% names(deseq_results))) next

  # Get per-gene summary for splicing (max |dPSI| event per gene)
  spl_per_gene <- splicing_results[[comp_name]] %>%
    filter(!is.na(deltapsi) & !is.na(FDR)) %>%
    group_by(GENE) %>%
    summarise(
      max_abs_dpsi = max(abs(deltapsi), na.rm = TRUE),
      best_dpsi    = deltapsi[which.max(abs(deltapsi))],
      best_fdr     = FDR[which.max(abs(deltapsi))],
      n_sig_events = sum(FDR <= SPLICING_FDR_THRESH & abs(deltapsi) >= SPLICING_DPSI_THRESH),
      .groups = "drop"
    ) %>%
    rename(gene = GENE) %>%
    mutate(
      splicing_sig = best_fdr <= SPLICING_FDR_THRESH & max_abs_dpsi >= SPLICING_DPSI_THRESH
    )

  # Get expression data
  expr_per_gene <- deseq_results[[comp_name]] %>%
    filter(!is.na(log2FoldChange) & !is.na(padj)) %>%
    dplyr::select(gene = gene_symbol, lfc = log2FoldChange, padj,
                  significance) %>%
    mutate(
      expression_sig = significance != "Not Significant"
    )

  # Merge
  merged <- inner_join(spl_per_gene, expr_per_gene, by = "gene") %>%
    mutate(
      category = case_when(
        splicing_sig & expression_sig  ~ "Both",
        splicing_sig & !expression_sig ~ "Splicing only",
        !splicing_sig & expression_sig ~ "Expression only",
        TRUE                           ~ "Neither"
      )
    )

  # Correlation for significant genes
  sig_merged <- merged %>% filter(category != "Neither")
  cor_val <- if (nrow(sig_merged) > 5) {
    cor(sig_merged$best_dpsi, sig_merged$lfc,
        use = "complete.obs", method = "spearman")
  } else NA

  # Label top genes (in "Both")
  top_both <- merged %>%
    filter(category == "Both") %>%
    arrange(desc(abs(best_dpsi) + abs(lfc))) %>%
    head(10)

  p <- ggplot(merged, aes(x = best_dpsi, y = lfc, color = category)) +
    geom_point(data = filter(merged, category == "Neither"),
               alpha = 0.1, size = 0.5) +
    geom_point(data = filter(merged, category != "Neither"),
               alpha = 0.7, size = 1.5) +
    geom_text_repel(
      data = top_both,
      aes(label = gene),
      size = 2.5, fontface = "italic", family = FONT_FAMILY,
      max.overlaps = 15, segment.size = 0.2, color = "black"
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_hline(yintercept = c(-EXPR_LFC_THRESH, EXPR_LFC_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-SPLICING_DPSI_THRESH, SPLICING_DPSI_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_color_manual(values = int_colors,
                       breaks = c("Both", "Splicing only", "Expression only")) +
    labs(
      title    = comparison_labels[comp_name],
      subtitle = ifelse(is.na(cor_val), "",
                          sprintf("Spearman r = %.2f (significant genes)", cor_val)),
      x = expression(Delta*PSI~"(best event per gene)"),
      y = expression(log[2]~"Fold Change (shrunken)"),
      color = NULL
    ) +
    theme_taxol(base_size = 11) +
    theme(legend.position = "right")

  scatter_plots[[comp_name]] <- p

  ggsave(file.path(PLOTS_DIR, sprintf("int_02_scatter_%s.pdf", comp_name)),
         p, width = 7, height = 6, device = SAVE_DEVICE)
  cat(sprintf("  Saved: int_02_scatter_%s.pdf\n", comp_name))
}

# ============================================================================
# ============================================================================
#  SECTION 3: UPSET PLOT — GENE OVERLAP ACROSS MODALITIES
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 3: UPSET PLOT — GENE OVERLAP\n")
cat(strrep("=", 70), "\n\n")

# For the taxol comparisons, make upset of splicing vs expression gene sets
taxol_comps <- intersect(c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO"),
                          names(splicing_results))
taxol_comps <- intersect(taxol_comps, names(deseq_results))

if (length(taxol_comps) >= 1) {
  upset_list <- list()

  for (comp_name in taxol_comps) {
    spl_g  <- get_splicing_genes(comp_name)
    expr_g <- get_expression_genes(comp_name)

    short <- gsub("_Taxol_vs_DMSO", "", comp_name)
    upset_list[[paste0(short, " Splicing")]]   <- spl_g
    upset_list[[paste0(short, " Expression")]] <- expr_g
  }

  # Only plot if at least 2 non-empty sets
  nonempty <- sapply(upset_list, length) > 0
  if (sum(nonempty) >= 2) {
    open_pdf(file.path(PLOTS_DIR, "int_03_upset_splicing_expression.pdf"),
             width = 9, height = 6)
    binary_mat <- fromList(upset_list[nonempty])
    UpSetR::upset(
      binary_mat,
      nsets        = length(upset_list[nonempty]),
      order.by     = "freq",
      point.size   = 3,
      line.size    = 1,
      text.scale   = c(1.3, 1.0, 1.0, 1.0, 1.3, 1.0),
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label    = "Genes per Set",
      mb.ratio = c(0.6, 0.4),
      sets.bar.color = c("#D62839", "#4BA3C3", "#B2182B", "#2166AC")[1:sum(nonempty)]
    )
    grid.text("Splicing vs Expression: Gene Overlap (Taxol Effect)",
              x = 0.5, y = 0.97,
              gp = gpar(fontsize = 12, fontface = "bold", fontfamily = FONT_FAMILY))
    dev.off()
    cat("  Saved: int_03_upset_splicing_expression.pdf\n")
  }
}

# ============================================================================
# ============================================================================
#  SECTION 4: STACKED BAR CHART — INTEGRATION CATEGORIES
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 4: INTEGRATION CATEGORY BAR CHART\n")
cat(strrep("=", 70), "\n\n")

if (nrow(int_combined) > 0) {

  bar_data <- int_combined %>%
    filter(category != "Neither") %>%
    group_by(comparison, category) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      comparison_label = comparison_labels[comparison],
      category = factor(category,
                         levels = c("Splicing only", "Both", "Expression only"))
    )

  p_bar <- ggplot(bar_data, aes(x = comparison_label, y = n, fill = category)) +
    geom_col(position = "dodge", width = 0.7, color = "grey30", linewidth = 0.2) +
    geom_text(aes(label = n), position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3, family = FONT_FAMILY) +
    scale_fill_manual(values = int_colors[c("Splicing only", "Both", "Expression only")]) +
    labs(
      title    = "Genes Affected by Splicing, Expression, or Both",
      subtitle = sprintf("Splicing: |dPSI| ≥ %.2f & FDR ≤ %.2f  |  Expression: |LFC| ≥ %.1f & padj ≤ %.2f",
                          SPLICING_DPSI_THRESH, SPLICING_FDR_THRESH,
                          EXPR_LFC_THRESH, EXPR_PADJ_THRESH),
      x = NULL, y = "Number of Genes", fill = NULL
    ) +
    theme_taxol() +
    theme(
      axis.text.x = element_text(angle = 15, hjust = 1, size = 9),
      legend.position = "top"
    )

  ggsave(file.path(PLOTS_DIR, "int_04_category_bar.pdf"), p_bar,
         width = 10, height = 5, device = SAVE_DEVICE)
  cat("  Saved: int_04_category_bar.pdf\n")

  # Proportional version
  prop_data <- bar_data %>%
    group_by(comparison_label) %>%
    mutate(total = sum(n), pct = 100 * n / total) %>%
    ungroup()

  p_prop <- ggplot(prop_data, aes(x = comparison_label, y = pct, fill = category)) +
    geom_col(position = "stack", width = 0.7, color = "grey30", linewidth = 0.2) +
    geom_text(
      aes(label = ifelse(pct > 8, sprintf("%d\n(%.0f%%)", n, pct), "")),
      position = position_stack(vjust = 0.5),
      size = 2.8, color = "white", fontface = "bold", family = FONT_FAMILY
    ) +
    scale_fill_manual(values = int_colors[c("Splicing only", "Both", "Expression only")]) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(
      title    = "Proportion: Splicing vs Expression",
      x = NULL, y = "% of Affected Genes", fill = NULL
    ) +
    theme_taxol() +
    theme(
      axis.text.x = element_text(angle = 15, hjust = 1, size = 9),
      legend.position = "bottom"
    )

  ggsave(file.path(PLOTS_DIR, "int_04b_category_proportions.pdf"), p_prop,
         width = 8, height = 5, device = SAVE_DEVICE)
  cat("  Saved: int_04b_category_proportions.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 5: HEATMAP — GENES WITH BOTH SPLICING & EXPRESSION CHANGES
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 5: DUAL-CHANGE GENES HEATMAP\n")
cat(strrep("=", 70), "\n\n")

# Identify genes that are significant in BOTH splicing and expression
# in at least one comparison

both_genes <- int_combined %>%
  filter(category == "Both") %>%
  pull(gene) %>%
  unique()

cat(sprintf("  Genes with both splicing & expression changes: %d\n",
            length(both_genes)))

if (length(both_genes) >= 5 && exists("vst_df")) {

  # Map gene symbols to ensembl IDs
  both_ensembl <- gene_map %>%
    filter(gene_symbol %in% both_genes) %>%
    pull(ensembl_id)

  # Get VST values
  vst_for_ht <- vst_df %>%
    filter(ensembl_id %in% both_ensembl) %>%
    column_to_rownames("ensembl_id")

  if ("gene_symbol" %in% colnames(vst_for_ht)) {
    row_labels <- vst_for_ht$gene_symbol
    vst_for_ht$gene_symbol <- NULL
  } else {
    row_labels <- rownames(vst_for_ht)
  }

  vst_mat_ht <- as.matrix(vst_for_ht[, metadata$sample_id])
  rownames(vst_mat_ht) <- row_labels

  # Remove rows with NA
  vst_mat_ht <- vst_mat_ht[complete.cases(vst_mat_ht), ]

  if (nrow(vst_mat_ht) > 3) {

    # Cap at 60 genes for readability
    if (nrow(vst_mat_ht) > 60) {
      # Prioritise genes with strongest combined effects
      var_order <- apply(vst_mat_ht, 1, var)
      vst_mat_ht <- vst_mat_ht[order(-var_order)[1:60], ]
    }

    # Scale rows
    vst_scaled <- t(scale(t(vst_mat_ht)))

    ha_top <- HeatmapAnnotation(
      Condition = metadata$condition,
      col = list(Condition = condition_colors),
      annotation_name_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
      annotation_legend_param = list(
        Condition = list(
          title_gp  = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
          labels_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY)
        )
      )
    )

    col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

    open_pdf(file.path(PLOTS_DIR, "int_05_both_heatmap.pdf"),
             width = 8, height = max(8, nrow(vst_scaled) * 0.15 + 3))
    ht <- Heatmap(
      vst_scaled,
      name = "Z-score",
      col  = col_fun,
      top_annotation = ha_top,
      column_title = "Genes with Both Splicing & Expression Changes (VST z-score)",
      column_title_gp = gpar(fontsize = 11, fontface = "bold",
                              fontfamily = FONT_FAMILY),
      row_names_gp = gpar(fontsize = 6, fontfamily = FONT_FAMILY,
                            fontface = "italic"),
      column_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
      show_column_names = TRUE,
      cluster_columns = FALSE,
      row_names_max_width = unit(8, "cm"),
      heatmap_legend_param = list(
        title_gp  = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
        labels_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY)
      )
    )
    draw(ht)
    dev.off()
    cat("  Saved: int_05_both_heatmap.pdf\n")
  }
}

# ============================================================================
# ============================================================================
#  SECTION 6: GO ENRICHMENT — CATEGORY-SPECIFIC
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 6: GO ENRICHMENT BY CATEGORY\n")
cat(strrep("=", 70), "\n\n")

# Enrichment for "Splicing only", "Expression only", "Both" genes
# Focus on the taxol comparisons

for (comp_name in c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO")) {
  if (!(comp_name %in% names(integration_summary))) next

  int_df <- integration_summary[[comp_name]]

  for (cat_label in c("Splicing only", "Expression only", "Both")) {
    cat_genes_sym <- int_df %>%
      filter(category == cat_label) %>%
      pull(gene)

    if (length(cat_genes_sym) < 5) {
      cat(sprintf("  %s — %s: only %d genes, skipping GO\n",
                  comp_name, cat_label, length(cat_genes_sym)))
      next
    }

    # Map to Ensembl for enrichGO
    cat_genes_ens <- gene_map %>%
      filter(gene_symbol %in% cat_genes_sym) %>%
      pull(ensembl_id)

    bg_genes_ens <- gene_map$ensembl_id

    ego <- tryCatch({
      enrichGO(
        gene          = cat_genes_ens,
        universe      = bg_genes_ens,
        OrgDb         = org.Mm.eg.db,
        keyType       = "ENSEMBL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.1,
        readable      = TRUE
      )
    }, error = function(e) NULL)

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      ego_df <- as.data.frame(ego)
      cat(sprintf("  %s — %s: %d GO terms\n", comp_name, cat_label, nrow(ego_df)))

      write.csv(ego_df,
                file.path(RESULTS_DIR,
                          sprintf("integrated_GO_%s_%s.csv",
                                  comp_name, gsub(" ", "_", cat_label))),
                row.names = FALSE)

      p_go <- dotplot(ego, showCategory = 12) +
        labs(title = sprintf("%s\n%s genes (%d)",
                              comparison_labels[comp_name], cat_label,
                              length(cat_genes_sym))) +
        theme_taxol(base_size = 9) +
        theme(axis.text.y = element_text(size = 7))

      ggsave(file.path(PLOTS_DIR,
                        sprintf("int_06_go_%s_%s.pdf",
                                comp_name, gsub(" ", "_", cat_label))),
             p_go, width = 8, height = 6, device = SAVE_DEVICE)
    } else {
      cat(sprintf("  %s — %s: no significant GO terms\n", comp_name, cat_label))
    }
  }
}

# ============================================================================
# ============================================================================
#  SECTION 7: COMBINED FIGURE — 2×2 INTEGRATION SUMMARY
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 7: COMBINED INTEGRATION FIGURE\n")
cat(strrep("=", 70), "\n\n")

# Build a combined figure showing the key integration results
if (length(scatter_plots) >= 2) {

  # Get 2 taxol scatter plots
  p1 <- scatter_plots[["Myoblast_Taxol_vs_DMSO"]] +
    theme(legend.position = "none", plot.title = element_text(size = 10))
  p2 <- scatter_plots[["Myotube_Taxol_vs_DMSO"]] +
    theme(legend.position = "none", plot.title = element_text(size = 10))

  # Category bar
  if (exists("p_bar")) {
    p3 <- p_bar + theme(plot.title = element_text(size = 10))
  } else {
    p3 <- ggplot() + theme_void()
  }

  # Proportion bar
  if (exists("p_prop")) {
    p4 <- p_prop + theme(plot.title = element_text(size = 10),
                          legend.position = "bottom")
  } else {
    p4 <- ggplot() + theme_void()
  }

  combined <- (p1 + p2) / (p3 + p4) +
    plot_layout(widths = c(1, 1), heights = c(1, 0.8)) +
    plot_annotation(
      title    = "Integrated Splicing × Expression Analysis — Taxol in C2C12",
      subtitle = "Upper: dPSI vs LFC per gene  |  Lower: category breakdown",
      tag_levels = "A",
      theme = theme(
        plot.title    = element_text(size = 14, face = "bold", hjust = 0.5,
                                      family = FONT_FAMILY),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40",
                                      family = FONT_FAMILY)
      )
    )

  ggsave(file.path(PLOTS_DIR, "int_07_combined_figure.pdf"), combined,
         width = 13, height = 10, device = SAVE_DEVICE)
  cat("  Saved: int_07_combined_figure.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 8: SPLICING vs EXPRESSION INDEPENDENCE TEST
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 8: INDEPENDENCE ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

# For genes present in both analyses, test whether |dPSI| is correlated
# with expression level or expression change

for (comp_name in c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO")) {
  if (!(comp_name %in% names(splicing_results)) ||
      !(comp_name %in% names(deseq_results))) next

  cat(sprintf("--- %s ---\n", comp_name))

  # Per-gene splicing summary
  spl_df <- splicing_results[[comp_name]] %>%
    filter(!is.na(deltapsi) & !is.na(FDR)) %>%
    group_by(GENE) %>%
    summarise(
      max_abs_dpsi = max(abs(deltapsi), na.rm = TRUE),
      n_events     = n(),
      n_sig        = sum(FDR <= SPLICING_FDR_THRESH & abs(deltapsi) >= SPLICING_DPSI_THRESH),
      .groups = "drop"
    ) %>%
    rename(gene = GENE)

  # Expression data
  expr_df <- deseq_results[[comp_name]] %>%
    dplyr::select(gene = gene_symbol, lfc = log2FoldChange,
                  baseMean, padj) %>%
    filter(!is.na(gene))

  merged <- inner_join(spl_df, expr_df, by = "gene")

  if (nrow(merged) > 10) {
    # Correlation: |dPSI| vs |LFC|
    cor_dpsi_lfc <- cor.test(merged$max_abs_dpsi, abs(merged$lfc),
                              method = "spearman")
    cat(sprintf("  |dPSI| vs |LFC|:    rho = %.3f, p = %.2e\n",
                cor_dpsi_lfc$estimate, cor_dpsi_lfc$p.value))

    # Correlation: |dPSI| vs expression level
    cor_dpsi_expr <- cor.test(merged$max_abs_dpsi, log10(merged$baseMean + 1),
                               method = "spearman")
    cat(sprintf("  |dPSI| vs baseMean: rho = %.3f, p = %.2e\n",
                cor_dpsi_expr$estimate, cor_dpsi_expr$p.value))

    # Plot: |dPSI| vs |LFC|
    p_ind <- ggplot(merged, aes(x = max_abs_dpsi, y = abs(lfc))) +
      geom_point(alpha = 0.15, size = 0.5, color = "grey50") +
      geom_smooth(method = "lm", color = "#762A83", linewidth = 0.8,
                  fill = alpha("#762A83", 0.2)) +
      labs(
        title    = comparison_labels[comp_name],
        subtitle = sprintf("Spearman rho = %.3f, p = %.2e",
                            cor_dpsi_lfc$estimate, cor_dpsi_lfc$p.value),
        x = expression("|"*Delta*PSI*"|"~"(best event per gene)"),
        y = expression("|log"[2]~"FC|")
      ) +
      theme_taxol(base_size = 11)

    ggsave(file.path(PLOTS_DIR, sprintf("int_08_independence_%s.pdf", comp_name)),
           p_ind, width = 6, height = 5, device = SAVE_DEVICE)
    cat(sprintf("  Saved: int_08_independence_%s.pdf\n", comp_name))
  }
}

# ============================================================================
# ============================================================================
#  SECTION 9: EXPRESSION OF DIFFERENTIALLY SPLICED GENES
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 9: EXPRESSION LEVELS OF SPLICED GENES\n")
cat(strrep("=", 70), "\n\n")

# Do differentially spliced genes tend to be highly expressed?

for (comp_name in c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO")) {
  if (!(comp_name %in% names(splicing_results)) ||
      !(comp_name %in% names(deseq_results))) next

  spl_genes <- get_splicing_genes(comp_name)

  expr_df <- deseq_results[[comp_name]] %>%
    filter(!is.na(baseMean) & !is.na(gene_symbol)) %>%
    mutate(
      is_spliced = gene_symbol %in% spl_genes,
      group = ifelse(is_spliced, "Differentially\nspliced", "Not\nspliced")
    )

  if (sum(expr_df$is_spliced) >= 5) {
    # Wilcoxon test
    wt <- wilcox.test(
      log10(expr_df$baseMean[expr_df$is_spliced] + 1),
      log10(expr_df$baseMean[!expr_df$is_spliced] + 1)
    )

    p_expr <- ggplot(expr_df, aes(x = group, y = log10(baseMean + 1), fill = group)) +
      geom_violin(alpha = 0.7, color = "grey30", linewidth = 0.3) +
      geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
      scale_fill_manual(values = c("Differentially\nspliced" = "#D62839",
                                     "Not\nspliced" = "grey70")) +
      labs(
        title    = sprintf("%s\nExpression of Spliced Genes",
                            comparison_labels[comp_name]),
        subtitle = sprintf("Wilcoxon p = %.2e", wt$p.value),
        x = NULL, y = expression(log[10]~"(baseMean + 1)")
      ) +
      theme_taxol(base_size = 11) +
      theme(legend.position = "none")

    ggsave(file.path(PLOTS_DIR, sprintf("int_09_expr_of_spliced_%s.pdf", comp_name)),
           p_expr, width = 4.5, height = 5, device = SAVE_DEVICE)
    cat(sprintf("  Saved: int_09_expr_of_spliced_%s.pdf\n", comp_name))
  }
}

# ============================================================================
# ============================================================================
#  SECTION 10: NUCLEAR AGITATION MODEL — EXPRESSION CHECK
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 10: NUCLEAR AGITATION — EXPRESSION CONTEXT\n")
cat(strrep("=", 70), "\n\n")

# If nuclear speckle and force gene lists exist, check their expression

speckle_file <- file.path(getwd(), "nuclear_speckle_associated_genes.txt")
force_file   <- file.path(getwd(), "final_all_layers.csv")

if (file.exists(speckle_file)) {
  speckle_genes <- scan(speckle_file, what = "character", quiet = TRUE)
  cat(sprintf("  Loaded %d nuclear speckle genes\n", length(speckle_genes)))

  for (comp_name in c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO")) {
    if (!(comp_name %in% names(deseq_results))) next

    expr_df <- deseq_results[[comp_name]] %>%
      filter(!is.na(gene_symbol)) %>%
      mutate(is_speckle = gene_symbol %in% speckle_genes)

    n_speckle_de <- sum(expr_df$is_speckle & expr_df$significance != "Not Significant",
                         na.rm = TRUE)
    cat(sprintf("  %s: %d speckle genes DE (out of %d tested)\n",
                comp_name, n_speckle_de, sum(expr_df$is_speckle)))
  }
}

if (file.exists(force_file)) {
  force_df <- read.csv(force_file)
  force_genes <- unique(force_df$GENE)
  cat(sprintf("  Loaded %d force-predicted genes\n", length(force_genes)))

  for (comp_name in c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO")) {
    if (!(comp_name %in% names(deseq_results))) next

    expr_df <- deseq_results[[comp_name]] %>%
      filter(!is.na(gene_symbol)) %>%
      mutate(is_force = gene_symbol %in% force_genes)

    n_force_de <- sum(expr_df$is_force & expr_df$significance != "Not Significant",
                       na.rm = TRUE)
    cat(sprintf("  %s: %d force genes DE (out of %d tested)\n",
                comp_name, n_force_de, sum(expr_df$is_force)))
  }
}

# ============================================================================
# FINAL OUTPUT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("INTEGRATED ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Output:\n")
cat(sprintf("  Results: %s\n", RESULTS_DIR))
cat(sprintf("  Plots:   %s\n", PLOTS_DIR))

if (exists("int_combined")) {
  for (comp_name in unique(int_combined$comparison)) {
    ic <- int_combined %>% filter(comparison == comp_name)
    cat(sprintf("  %s: Both=%d, Splicing=%d, Expression=%d\n",
                comp_name,
                sum(ic$category == "Both"),
                sum(ic$category == "Splicing only"),
                sum(ic$category == "Expression only")))
  }
}

cat("\nKey insight: If taxol acts primarily through nucleoplasmic agitation\n")
cat("(mechanical perturbation of splicing), we expect:\n")
cat("  - Splicing changes largely INDEPENDENT of expression changes\n")
cat("  - Low correlation between |dPSI| and |LFC|\n")
cat("  - Many 'Splicing only' genes relative to 'Both'\n")

cat("\nSession info:\n")
sessionInfo()
