# ============================================================================
# 03_gene_expression.R
# Taxol Effect on Gene Expression in C2C12 Myoblasts and Myotubes
# DESeq2 Analysis with Variance Stabilisation & LFC Shrinkage
#
# PURPOSE: Perform differential gene expression analysis on Salmon counts
#          using DESeq2. Apply variance-stabilising transformation (VST) for
#          QC visualisations and apeglm LFC shrinkage for accurate effect
#          sizes. Output tables and publication-ready figures that feed into
#          the integrated splicing + expression analysis (05_integrated.R).
#
# EXPERIMENTAL DESIGN:
#   12 samples — 4 conditions × 3 replicates:
#     Myoblast_DMSO, Myoblast_Taxol, Myotube_DMSO, Myotube_Taxol
#
# COMPARISONS:
#   1. Myoblast_Taxol_vs_DMSO   — Taxol effect in undifferentiated cells
#   2. Myotube_Taxol_vs_DMSO    — Taxol effect in differentiated cells
#   3. Differentiation_DMSO     — Differentiation without Taxol
#   4. Differentiation_Taxol    — Differentiation with Taxol
#
# INPUT:
#   Salmon_Gene_table_counts.csv   (gene-level raw counts from Salmon)
#   metadata/metadata.csv
#
# OUTPUT:
#   results/deseq2_<comparison>.csv     — full shrunken results per comparison
#   results/deseq2_vst_matrix.csv       — VST-normalised expression matrix
#   plots/ge_*                          — publication-ready figures
#
# Author: Andrés Gordo Ortiz
# ============================================================================

cat(strrep("=", 70), "\n")
cat("03_gene_expression.R — DESeq2 Gene Expression Analysis\n")
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
  library(pheatmap)
  library(viridis)
  library(scales)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichR)
  library(apeglm)
  library(ggpubr)
  library(EnhancedVolcano)
})

# ============================================================================
# FONT SETUP (consistent with 02_splicing_analysis.R)
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

de_colors <- c(
  "Up"              = "#D62839",
  "Down"            = "#4BA3C3",
  "Not Significant" = "grey80"
)

# ============================================================================
# PATHS
# ============================================================================

RESULTS_DIR   <- file.path(getwd(), "results")
PLOTS_DIR     <- file.path(getwd(), "plots")
METADATA_FILE <- file.path(getwd(), "metadata", "metadata.csv")
COUNTS_FILE   <- file.path(getwd(), "Salmon_Gene_table_counts.csv")

if (!dir.exists(PLOTS_DIR))   dir.create(PLOTS_DIR, recursive = TRUE)
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ============================================================================
# LOAD & PREPARE DATA
# ============================================================================

cat("--- Loading Salmon gene counts ---\n")
counts_raw <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
counts_int <- round(counts_raw)

metadata <- read.csv(METADATA_FILE, stringsAsFactors = FALSE) %>%
  arrange(condition)

# Ensure column order matches metadata
stopifnot(all(metadata$sample_id %in% colnames(counts_int)))
counts_int <- counts_int[, metadata$sample_id]

cat(sprintf("  Loaded: %d genes × %d samples\n", nrow(counts_int), ncol(counts_int)))

# Pre-filter: keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts_int >= 10) >= 3
counts_filtered <- counts_int[keep, ]
cat(sprintf("  After filtering: %d genes\n", nrow(counts_filtered)))

# ============================================================================
# GENE SYMBOL ANNOTATION (Ensembl → Symbol)
# ============================================================================

cat("--- Annotating Ensembl IDs → gene symbols ---\n")

gene_symbols <- tryCatch({
  AnnotationDbi::mapIds(org.Mm.eg.db,
                        keys     = rownames(counts_filtered),
                        column   = "SYMBOL",
                        keytype  = "ENSEMBL",
                        multiVals = "first")
}, error = function(e) {
  cat("  Warning: org.Mm.eg.db annotation failed, using Ensembl IDs\n")
  setNames(rownames(counts_filtered), rownames(counts_filtered))
})

gene_map <- data.frame(
  ensembl_id  = names(gene_symbols),
  gene_symbol = unname(gene_symbols),
  stringsAsFactors = FALSE
) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_id, gene_symbol))

# ============================================================================
# ============================================================================
#  DESeq2: FULL MODEL
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("BUILDING DESeq2 DATASET\n")
cat(strrep("=", 70), "\n\n")

# Use condition as the single factor — 4 levels
col_data <- data.frame(
  row.names = metadata$sample_id,
  condition = factor(metadata$condition,
                     levels = c("Myoblast_DMSO", "Myoblast_Taxol",
                                "Myotube_DMSO", "Myotube_Taxol")),
  cell_type = factor(metadata$cell_type),
  treatment = factor(metadata$treatment)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = col_data,
  design    = ~ condition
)

# Run DESeq2
dds <- DESeq(dds, quiet = FALSE)

cat("\n  Dispersion estimation complete.\n")
cat(sprintf("  Genes tested: %d\n", nrow(dds)))

# ============================================================================
# VARIANCE STABILISING TRANSFORMATION (for QC)
# ============================================================================

cat("\n--- Variance-stabilising transformation (VST) ---\n")
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

# Save VST matrix
vst_out <- as.data.frame(vst_mat) %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_map, by = "ensembl_id")
write.csv(vst_out, file.path(RESULTS_DIR, "deseq2_vst_matrix.csv"),
          row.names = FALSE)
cat("  Saved: deseq2_vst_matrix.csv\n")

# ============================================================================
# ============================================================================
#  SECTION 1: QC PLOTS
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 1: QUALITY CONTROL\n")
cat(strrep("=", 70), "\n\n")

# ----------------------------------------------------------------------------
# 1.1 PCA on VST-normalised counts
# ----------------------------------------------------------------------------

cat("--- 1.1 PCA on VST expression ---\n")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_data$Condition <- pca_data$condition
percentVar <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = condition_colors) +
  labs(
    title    = "Gene Expression PCA (VST)",
    subtitle = "12 samples — 4 conditions × 3 replicates",
    x = sprintf("PC1 (%s%% variance)", percentVar[1]),
    y = sprintf("PC2 (%s%% variance)", percentVar[2]),
    color = NULL
  ) +
  coord_fixed() +
  theme_taxol(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(PLOTS_DIR, "ge_01_pca_vst.pdf"), p_pca,
       width = 7, height = 5, device = SAVE_DEVICE)
cat("  Saved: ge_01_pca_vst.pdf\n")

# ----------------------------------------------------------------------------
# 1.2 Sample-to-sample distance heatmap
# ----------------------------------------------------------------------------

cat("--- 1.2 Sample distance heatmap ---\n")

sample_dists <- dist(t(vst_mat))
dist_mat <- as.matrix(sample_dists)

# Annotation
ha <- HeatmapAnnotation(
  Condition = col_data$condition,
  col = list(Condition = condition_colors),
  annotation_name_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
  annotation_legend_param = list(
    Condition = list(title_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
                     labels_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY))
  )
)

col_fun <- colorRamp2(
  c(min(dist_mat), median(dist_mat), max(dist_mat)),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

open_pdf(file.path(PLOTS_DIR, "ge_02_sample_distances.pdf"), width = 7, height = 6)
ht <- Heatmap(dist_mat,
              name = "Distance",
              col  = col_fun,
              top_annotation = ha,
              column_title = "Sample-to-Sample Euclidean Distance (VST)",
              column_title_gp = gpar(fontsize = 12, fontface = "bold",
                                      fontfamily = FONT_FAMILY),
              row_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
              column_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
                labels_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY)
              ))
draw(ht)
dev.off()
cat("  Saved: ge_02_sample_distances.pdf\n")

# ============================================================================
# ============================================================================
#  SECTION 2: DIFFERENTIAL EXPRESSION — ALL 4 COMPARISONS
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 2: DIFFERENTIAL EXPRESSION (LFC SHRINKAGE)\n")
cat(strrep("=", 70), "\n\n")

# Define contrasts (test vs reference)
comparisons <- list(
  Myoblast_Taxol_vs_DMSO = c("condition", "Myoblast_Taxol", "Myoblast_DMSO"),
  Myotube_Taxol_vs_DMSO  = c("condition", "Myotube_Taxol",  "Myotube_DMSO"),
  Differentiation_DMSO   = c("condition", "Myotube_DMSO",   "Myoblast_DMSO"),
  Differentiation_Taxol  = c("condition", "Myotube_Taxol",  "Myoblast_Taxol")
)

# Store results
deseq_results <- list()

LFC_THRESH  <- 1.0
PADJ_THRESH <- 0.05

for (comp_name in names(comparisons)) {
  cat(sprintf("--- %s ---\n", comp_name))

  # Unshrunken results (for p-values)
  res_raw <- results(dds, contrast = comparisons[[comp_name]], alpha = PADJ_THRESH)

  # LFC shrinkage with apeglm requires a coefficient name
  # Build the coefficient name from the contrast
  coef_name <- paste0("condition_",
                       comparisons[[comp_name]][2], "_vs_",
                       comparisons[[comp_name]][3])

  # If this exact coefficient does not exist we need to use the manual approach
  if (coef_name %in% resultsNames(dds)) {
    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  } else {
    # Fallback: use ashr shrinkage which accepts contrast vectors
    res_shrunk <- lfcShrink(dds, contrast = comparisons[[comp_name]],
                            type = "ashr", quiet = TRUE)
  }

  # Build output table
  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("ensembl_id") %>%
    left_join(gene_map, by = "ensembl_id") %>%
    mutate(
      significance = case_when(
        padj < PADJ_THRESH & log2FoldChange >=  LFC_THRESH ~ "Up",
        padj < PADJ_THRESH & log2FoldChange <= -LFC_THRESH ~ "Down",
        TRUE ~ "Not Significant"
      )
    ) %>%
    arrange(padj)

  # Summary
  n_up   <- sum(res_df$significance == "Up",   na.rm = TRUE)
  n_down <- sum(res_df$significance == "Down", na.rm = TRUE)
  cat(sprintf("  DE genes: %d up, %d down (|LFC| >= %.1f, padj < %.2f)\n",
              n_up, n_down, LFC_THRESH, PADJ_THRESH))

  # Save
  write.csv(res_df, file.path(RESULTS_DIR, paste0("deseq2_", comp_name, ".csv")),
            row.names = FALSE)
  cat(sprintf("  Saved: deseq2_%s.csv\n", comp_name))

  deseq_results[[comp_name]] <- res_df
}

# ============================================================================
# ============================================================================
#  SECTION 3: VOLCANO PLOTS
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 3: VOLCANO PLOTS\n")
cat(strrep("=", 70), "\n\n")

comparison_labels <- c(
  "Myoblast_Taxol_vs_DMSO" = "Myoblast: Taxol vs DMSO",
  "Myotube_Taxol_vs_DMSO"  = "Myotube: Taxol vs DMSO",
  "Differentiation_DMSO"   = "Differentiation (DMSO)",
  "Differentiation_Taxol"  = "Differentiation (Taxol)"
)

make_volcano_ge <- function(df, title_str) {
  df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange))

  # Cap -log10(padj)
  max_y <- 50
  df$negLogPadj <- pmin(-log10(df$padj), max_y)

  # Label top genes
  top_genes <- df %>%
    filter(significance != "Not Significant") %>%
    arrange(padj) %>%
    head(15)

  n_up   <- sum(df$significance == "Up")
  n_down <- sum(df$significance == "Down")

  ggplot(df, aes(x = log2FoldChange, y = negLogPadj, color = significance)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_point(data = filter(df, significance != "Not Significant"),
               alpha = 0.7, size = 1.3) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene_symbol),
      size         = 2.5,
      fontface     = "italic",
      family       = FONT_FAMILY,
      max.overlaps = 20,
      segment.size = 0.25,
      box.padding  = 0.3,
      color        = "black"
    ) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(PADJ_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_color_manual(values = de_colors,
                       breaks = c("Up", "Down")) +
    labs(
      title    = title_str,
      subtitle = sprintf("Up: %d  |  Down: %d  (|LFC| ≥ %.1f, padj < %.2f)",
                          n_up, n_down, LFC_THRESH, PADJ_THRESH),
      x = expression(log[2]~"Fold Change (shrunken)"),
      y = expression(-log[10]~"adjusted p-value"),
      color = NULL
    ) +
    theme_taxol(base_size = 11) +
    theme(legend.position = "top")
}

for (comp_name in names(comparison_labels)) {
  if (comp_name %in% names(deseq_results)) {
    p <- make_volcano_ge(deseq_results[[comp_name]], comparison_labels[comp_name])
    fname <- sprintf("ge_03_volcano_%s.pdf", comp_name)
    ggsave(file.path(PLOTS_DIR, fname), p, width = 7, height = 6, device = SAVE_DEVICE)
    cat(sprintf("  Saved: %s\n", fname))
  }
}

# ============================================================================
# ============================================================================
#  SECTION 4: MA PLOTS
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 4: MA PLOTS\n")
cat(strrep("=", 70), "\n\n")

make_ma_plot <- function(df, title_str) {
  df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange) & !is.na(baseMean))
  df$log10BaseMean <- log10(df$baseMean + 1)

  ggplot(df, aes(x = log10BaseMean, y = log2FoldChange, color = significance)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_point(data = filter(df, significance != "Not Significant"),
               alpha = 0.6, size = 1.0) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_hline(yintercept = c(-LFC_THRESH, LFC_THRESH),
               linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_color_manual(values = de_colors, breaks = c("Up", "Down")) +
    labs(
      title    = title_str,
      x = expression(log[10]~"Mean Expression"),
      y = expression(log[2]~"Fold Change (shrunken)"),
      color = NULL
    ) +
    theme_taxol(base_size = 11) +
    theme(legend.position = "top")
}

for (comp_name in names(comparison_labels)) {
  if (comp_name %in% names(deseq_results)) {
    p <- make_ma_plot(deseq_results[[comp_name]], comparison_labels[comp_name])
    fname <- sprintf("ge_04_ma_%s.pdf", comp_name)
    ggsave(file.path(PLOTS_DIR, fname), p, width = 7, height = 5, device = SAVE_DEVICE)
    cat(sprintf("  Saved: %s\n", fname))
  }
}

# ============================================================================
# ============================================================================
#  SECTION 5: DE SUMMARY BAR CHART
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 5: DE GENE SUMMARY\n")
cat(strrep("=", 70), "\n\n")

de_summary <- bind_rows(
  lapply(names(deseq_results), function(n) {
    df <- deseq_results[[n]]
    data.frame(
      comparison = n,
      comparison_label = comparison_labels[n],
      direction = c("Up", "Down"),
      count = c(
        sum(df$significance == "Up", na.rm = TRUE),
        sum(df$significance == "Down", na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    )
  })
)

p_summary <- ggplot(de_summary, aes(x = comparison_label, y = count, fill = direction)) +
  geom_col(position = "dodge", width = 0.7, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3, family = FONT_FAMILY) +
  scale_fill_manual(values = c("Up" = "#D62839", "Down" = "#4BA3C3")) +
  labs(
    title    = "Differentially Expressed Genes per Comparison",
    subtitle = sprintf("|LFC| ≥ %.1f, padj < %.2f (apeglm/ashr shrinkage)",
                       LFC_THRESH, PADJ_THRESH),
    x = NULL, y = "Number of DE Genes", fill = NULL
  ) +
  theme_taxol() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9))

ggsave(file.path(PLOTS_DIR, "ge_05_de_summary.pdf"), p_summary,
       width = 9, height = 5, device = SAVE_DEVICE)
cat("  Saved: ge_05_de_summary.pdf\n")

# ============================================================================
# ============================================================================
#  SECTION 6: TAXOL EFFECT SCATTER (MYOBLAST vs MYOTUBE)
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 6: TAXOL EFFECT SCATTER\n")
cat(strrep("=", 70), "\n\n")

if (all(c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO") %in%
        names(deseq_results))) {

  myo_df <- deseq_results[["Myoblast_Taxol_vs_DMSO"]] %>%
    dplyr::select(ensembl_id, gene_symbol,
                  lfc_myo = log2FoldChange, padj_myo = padj,
                  sig_myo = significance)

  tube_df <- deseq_results[["Myotube_Taxol_vs_DMSO"]] %>%
    dplyr::select(ensembl_id,
                  lfc_tube = log2FoldChange, padj_tube = padj,
                  sig_tube = significance)

  scatter_df <- inner_join(myo_df, tube_df, by = "ensembl_id") %>%
    mutate(
      category = case_when(
        sig_myo != "Not Significant" & sig_tube != "Not Significant" ~ "Both",
        sig_myo != "Not Significant" ~ "Myoblast only",
        sig_tube != "Not Significant" ~ "Myotube only",
        TRUE ~ "NS"
      )
    )

  cat_colors <- c("Both" = "#762A83", "Myoblast only" = "#D62839",
                   "Myotube only" = "#2166AC", "NS" = "grey85")

  cor_val <- cor(scatter_df$lfc_myo, scatter_df$lfc_tube,
                 use = "complete.obs", method = "pearson")

  p_scatter <- ggplot(scatter_df, aes(x = lfc_myo, y = lfc_tube, color = category)) +
    geom_point(data = filter(scatter_df, category == "NS"),
               alpha = 0.1, size = 0.5) +
    geom_point(data = filter(scatter_df, category != "NS"),
               alpha = 0.7, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, color = "red",
                linetype = "dotted", linewidth = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = cat_colors,
                       breaks = c("Both", "Myoblast only", "Myotube only")) +
    labs(
      title    = "Taxol Effect on Gene Expression",
      subtitle = sprintf("Myoblast vs Myotube LFC (r = %.2f)", cor_val),
      x = expression(log[2]~"FC — Myoblast (Taxol vs DMSO)"),
      y = expression(log[2]~"FC — Myotube (Taxol vs DMSO)"),
      color = NULL
    ) +
    coord_fixed(ratio = 1) +
    theme_taxol(base_size = 11) +
    theme(legend.position = "right")

  ggsave(file.path(PLOTS_DIR, "ge_06_taxol_scatter.pdf"), p_scatter,
         width = 7, height = 7, device = SAVE_DEVICE)
  cat("  Saved: ge_06_taxol_scatter.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 7: TOP DE GENES HEATMAP
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 7: TOP DE GENES HEATMAP\n")
cat(strrep("=", 70), "\n\n")

# Collect top 50 DE genes across all comparisons (by padj)
top_genes_all <- bind_rows(
  lapply(names(deseq_results), function(n) {
    deseq_results[[n]] %>%
      filter(significance != "Not Significant") %>%
      arrange(padj) %>%
      head(30) %>%
      dplyr::select(ensembl_id, gene_symbol)
  })
) %>%
  distinct(ensembl_id, .keep_all = TRUE)

if (nrow(top_genes_all) > 5) {
  # Extract VST values for these genes
  ht_mat <- vst_mat[top_genes_all$ensembl_id[top_genes_all$ensembl_id %in%
                                                rownames(vst_mat)], ]

  # Label rows with gene symbols
  row_labels <- gene_map$gene_symbol[match(rownames(ht_mat), gene_map$ensembl_id)]
  row_labels[is.na(row_labels)] <- rownames(ht_mat)[is.na(row_labels)]

  # Scale rows (z-score)
  ht_scaled <- t(scale(t(ht_mat)))

  # Annotation
  ha_top <- HeatmapAnnotation(
    Condition = col_data[colnames(ht_scaled), "condition"],
    col = list(Condition = condition_colors),
    annotation_name_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
    annotation_legend_param = list(
      Condition = list(title_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
                       labels_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY))
    )
  )

  col_fun_z <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

  open_pdf(file.path(PLOTS_DIR, "ge_07_top_de_heatmap.pdf"), width = 8, height = 12)
  ht_de <- Heatmap(ht_scaled,
                   name = "Z-score",
                   col  = col_fun_z,
                   top_annotation = ha_top,
                   row_labels = row_labels,
                   column_title = "Top Differentially Expressed Genes (VST z-score)",
                   column_title_gp = gpar(fontsize = 12, fontface = "bold",
                                           fontfamily = FONT_FAMILY),
                   row_names_gp = gpar(fontsize = 6, fontfamily = FONT_FAMILY,
                                        fontface = "italic"),
                   column_names_gp = gpar(fontsize = 7, fontfamily = FONT_FAMILY),
                   show_column_names = TRUE,
                   cluster_columns = FALSE,
                   row_names_max_width = unit(8, "cm"),
                   heatmap_legend_param = list(
                     title_gp = gpar(fontsize = 9, fontfamily = FONT_FAMILY),
                     labels_gp = gpar(fontsize = 8, fontfamily = FONT_FAMILY)
                   ))
  draw(ht_de)
  dev.off()
  cat("  Saved: ge_07_top_de_heatmap.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 8: GO ENRICHMENT FOR TAXOL-RESPONSIVE GENES
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 8: GO ENRICHMENT — TAXOL COMPARISONS\n")
cat(strrep("=", 70), "\n\n")

run_go_enrichment <- function(de_genes_ensembl, bg_genes_ensembl, comp_label, direction_label) {
  ego <- tryCatch({
    enrichGO(
      gene          = de_genes_ensembl,
      universe      = bg_genes_ensembl,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENSEMBL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1,
      readable      = TRUE
    )
  }, error = function(e) {
    cat(sprintf("    GO enrichment failed for %s %s: %s\n",
                comp_label, direction_label, e$message))
    return(NULL)
  })

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_df <- as.data.frame(ego) %>%
      mutate(comparison = comp_label, direction = direction_label)
    write.csv(ego_df,
              file.path(RESULTS_DIR,
                        sprintf("deseq2_GO_BP_%s_%s.csv", comp_label, direction_label)),
              row.names = FALSE)
    cat(sprintf("    %s %s: %d significant GO terms\n",
                comp_label, direction_label, nrow(ego_df)))
    return(ego)
  } else {
    cat(sprintf("    %s %s: no significant GO terms\n", comp_label, direction_label))
    return(NULL)
  }
}

bg_genes <- rownames(counts_filtered)

go_plots <- list()
for (comp_name in c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO")) {
  df <- deseq_results[[comp_name]]

  for (dir in c("Up", "Down")) {
    de_genes <- df %>% filter(significance == dir) %>% pull(ensembl_id)
    if (length(de_genes) >= 5) {
      ego <- run_go_enrichment(de_genes, bg_genes, comp_name, dir)
      if (!is.null(ego)) {
        p <- dotplot(ego, showCategory = 15) +
          labs(title = sprintf("%s — %s-regulated genes",
                               comparison_labels[comp_name], dir)) +
          theme_taxol(base_size = 9) +
          theme(axis.text.y = element_text(size = 7))
        go_plots[[paste0(comp_name, "_", dir)]] <- p
      }
    }
  }
}

# Combine GO plots if we have any
if (length(go_plots) > 0) {
  combined_go <- wrap_plots(go_plots, ncol = 2)
  ggsave(file.path(PLOTS_DIR, "ge_08_go_enrichment.pdf"), combined_go,
         width = 14, height = 7 * ceiling(length(go_plots) / 2),
         device = SAVE_DEVICE, limitsize = FALSE)
  cat("  Saved: ge_08_go_enrichment.pdf\n")
}

# ============================================================================
# ============================================================================
#  SECTION 9: COMBINED FIGURE (2×2)
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 9: COMBINED FIGURE\n")
cat(strrep("=", 70), "\n\n")

if (all(c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO") %in%
        names(deseq_results))) {

  v_myo  <- make_volcano_ge(deseq_results[["Myoblast_Taxol_vs_DMSO"]],
                             "Myoblast: Taxol vs DMSO") +
    theme(legend.position = "none")
  v_tube <- make_volcano_ge(deseq_results[["Myotube_Taxol_vs_DMSO"]],
                             "Myotube: Taxol vs DMSO") +
    theme(legend.position = "none")

  combined <- (v_myo + v_tube) / (p_pca + p_summary) +
    plot_layout(widths = c(1, 1), heights = c(1, 1)) +
    plot_annotation(
      title    = "Taxol Effect on Gene Expression — C2C12",
      subtitle = "DESeq2 with apeglm/ashr LFC shrinkage",
      theme = theme(
        plot.title    = element_text(size = 16, face = "bold", hjust = 0.5,
                                      family = FONT_FAMILY),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40",
                                      family = FONT_FAMILY)
      )
    )

  ggsave(file.path(PLOTS_DIR, "ge_09_combined_figure.pdf"), combined,
         width = 12, height = 10, device = SAVE_DEVICE)
  cat("  Saved: ge_09_combined_figure.pdf\n")
}

# ============================================================================
# SAVE SESSION & SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("GENE EXPRESSION ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Output:\n")
cat(sprintf("  Results:  %s\n", RESULTS_DIR))
cat(sprintf("  Plots:    %s\n", PLOTS_DIR))

for (comp_name in names(deseq_results)) {
  df <- deseq_results[[comp_name]]
  n_up   <- sum(df$significance == "Up",   na.rm = TRUE)
  n_down <- sum(df$significance == "Down", na.rm = TRUE)
  cat(sprintf("  %s: %d up, %d down\n", comp_name, n_up, n_down))
}

cat("\nSession info:\n")
sessionInfo()
