# ============================================================================
# 04_matt_feature_analysis.R
# Taxol Effect on Exon Features — Matt Output Analysis
# Publication-Ready Dot Matrix Plots
#
# PURPOSE: Analyse exon features from Matt output files. Matt compares
#          structural/sequence features of included vs skipped vs unchanged
#          exons. This script reads the 4 comparison-level Matt outputs,
#          performs Wilcoxon tests (up/down vs ndiff), and generates
#          publication-ready dot matrix plots showing which features
#          discriminate taxol-affected and differentiation-affected exons.
#
# INPUT FILES (place in matt_out/ directory):
#   matt_out/Differentiation_DMSO_exons_with_efeatures.tab
#   matt_out/Differentiation_Taxol_exons_with_efeatures.tab
#   matt_out/Myoblast_Taxol_vs_DMSO_exons_with_efeatures.tab
#   matt_out/Myotube_Taxol_vs_DMSO_exons_with_efeatures.tab
#
# OUTPUT:
#   plots/matt_dot_matrix_exons.pdf      — main dot matrix figure
#   plots/matt_feature_boxplots.pdf      — detailed boxplots for key features
#   results/matt_feature_statistics.csv  — full statistical results
#   results/matt_effect_directions.csv   — effect directions
#
# Author: Andrés Gordo Ortiz
# ============================================================================

cat(strrep("=", 70), "\n")
cat("04_matt_feature_analysis.R — Matt Exon Feature Analysis\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(showtext)
  library(sysfonts)
  library(forcats)
  library(patchwork)
  library(scales)
  library(ggpubr)
  library(viridis)
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

# Publication-ready matrix theme (for dot matrix plot)
theme_matrix <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = FONT_FAMILY) %+replace%
    theme(
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.background   = element_rect(fill = "white", color = NA),
      panel.border       = element_blank(),
      axis.text.x        = element_blank(),
      axis.text.y        = element_text(size = rel(1.0), hjust = 1,
                                         color = "grey20"),
      axis.title         = element_blank(),
      axis.ticks         = element_line(color = "grey50", linewidth = 0.3),
      axis.ticks.x       = element_blank(),
      axis.ticks.length  = unit(0.12, "cm"),
      strip.text         = element_text(face = "bold", size = rel(1.1)),
      strip.background   = element_rect(fill = "grey95", color = "grey70",
                                          linewidth = 0.3),
      legend.position    = "none",
      plot.title         = element_text(face = "bold", size = rel(1.35),
                                         hjust = 0.5,
                                         margin = margin(b = 8, t = 2)),
      plot.subtitle      = element_text(size = rel(1.0), hjust = 0.5,
                                         color = "grey40",
                                         margin = margin(b = 10)),
      plot.margin        = margin(25, 20, 10, 12)
    )
}

# ============================================================================
# COLOR PALETTES
# ============================================================================

# Comparison colors
comparison_colors <- c(
  "Myoblast: Taxol vs DMSO" = "#D62839",
  "Myotube: Taxol vs DMSO"  = "#B2182B",
  "Differentiation (DMSO)"  = "#4BA3C3",
  "Differentiation (Taxol)" = "#2166AC"
)

# Direction colors
direction_colors <- c(
  "Increase" = "#C85450",
  "Decrease" = "#5B8FA3"
)

# Event direction colors
event_colors <- c(
  "Included" = "#D4A574",
  "Skipped"  = "#B85450"
)

# ============================================================================
# PATHS
# ============================================================================

RESULTS_DIR  <- file.path(getwd(), "results")
PLOTS_DIR    <- file.path(getwd(), "plots")
MATT_DIR     <- file.path(getwd())

if (!dir.exists(PLOTS_DIR))   dir.create(PLOTS_DIR, recursive = TRUE)
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ============================================================================
# LOAD MATT OUTPUT FILES
# ============================================================================

cat("--- Loading Matt exon feature files ---\n")

matt_files <- list(
  "Myoblast_Taxol_vs_DMSO" = file.path(MATT_DIR, "Myoblast_Taxol_vs_DMSO_exons_with_efeatures.tab"),
  "Myotube_Taxol_vs_DMSO"  = file.path(MATT_DIR, "Myotube_Taxol_vs_DMSO_exons_with_efeatures.tab"),
  "Differentiation_DMSO"   = file.path(MATT_DIR, "Differentiation_DMSO_exons_with_efeatures.tab"),
  "Differentiation_Taxol"  = file.path(MATT_DIR, "Differentiation_Taxol_exons_with_efeatures.tab")
)

# Nice labels for plots
comparison_labels <- c(
  "Myoblast_Taxol_vs_DMSO" = "Myoblast:\nTaxol vs DMSO",
  "Myotube_Taxol_vs_DMSO"  = "Myotube:\nTaxol vs DMSO",
  "Differentiation_DMSO"   = "Differentiation\n(DMSO)",
  "Differentiation_Taxol"  = "Differentiation\n(Taxol)"
)

# Check which files exist
files_exist <- sapply(matt_files, file.exists)
if (!any(files_exist)) {
  stop("No Matt output files found in ", MATT_DIR,
       "\nExpected files: ", paste(basename(unlist(matt_files)), collapse = ", "))
}

cat(sprintf("  Found %d / %d Matt files\n", sum(files_exist), length(matt_files)))

# Load available files
matt_data <- list()
for (name in names(matt_files)[files_exist]) {
  matt_data[[name]] <- read_tsv(matt_files[[name]], show_col_types = FALSE)
  cat(sprintf("  Loaded: %s (%d events)\n", basename(matt_files[[name]]),
              nrow(matt_data[[name]])))
}

# ============================================================================
# ============================================================================
#  SECTION 1: FEATURE STATISTICS (Wilcoxon: up/down vs ndiff)
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 1: FEATURE STATISTICAL TESTS\n")
cat(strrep("=", 70), "\n\n")

# Columns to exclude from feature analysis (identifiers, sequences, etc.)
exclude_cols <- c("GENE", "EVENT", "START", "END", "SCAFFOLD", "STRAND",
                  "LENGTH", "COMPLEX", "GROUP", "DATASET", "FDR", "DELTAPSI",
                  "GENEID_ENSEMBL", "ENSEMBL_GENEID",
                  "EXON_ID", "GENE_BIOTYPE", "EXON_FOUND_IN_GTF",
                  "EXON_FOUND_IN_THESE_TRS", "EXON_COOCCURS_WITH_THESE_EXONS",
                  "EXON_COOCCURS_WITH_OTHER_EXONS",
                  grep("^SEQ_", colnames(matt_data[[1]]), value = TRUE))

# Function to compute statistics for one comparison
compute_feature_stats <- function(df, comp_name) {
  # Identify numeric feature columns
  numeric_cols <- colnames(df)[sapply(df, is.numeric)]
  feature_cols <- setdiff(numeric_cols, exclude_cols)

  results_list <- list()

  for (feat in feature_cols) {
    for (direction in c("up", "down")) {
      vals_dir  <- df %>% filter(GROUP == direction) %>% pull(!!sym(feat))
      vals_bg   <- df %>% filter(GROUP == "ndiff")    %>% pull(!!sym(feat))

      # Remove NAs
      vals_dir <- vals_dir[!is.na(vals_dir)]
      vals_bg  <- vals_bg[!is.na(vals_bg)]

      if (length(vals_dir) >= 3 && length(vals_bg) >= 3) {
        wt <- tryCatch(
          wilcox.test(vals_dir, vals_bg, exact = FALSE),
          error = function(e) list(p.value = NA, statistic = NA)
        )

        med_dir <- median(vals_dir, na.rm = TRUE)
        med_bg  <- median(vals_bg,  na.rm = TRUE)

        results_list[[paste(feat, direction, sep = "_")]] <- data.frame(
          feature    = feat,
          comparison = direction,
          treatment  = comp_name,
          n_dir      = length(vals_dir),
          n_bg       = length(vals_bg),
          median_dir = med_dir,
          median_bg  = med_bg,
          diff       = med_dir - med_bg,
          log2fc     = log2((med_dir + 0.01) / (med_bg + 0.01)),
          p.value    = wt$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  bind_rows(results_list)
}

# Run for all comparisons
all_stats <- bind_rows(
  lapply(names(matt_data), function(n) {
    cat(sprintf("  Computing stats for %s...\n", n))
    compute_feature_stats(matt_data[[n]], n)
  })
)

# Multiple testing correction (BH within each comparison)
all_stats <- all_stats %>%
  group_by(treatment) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    sig_level = case_when(
      is.na(p.value) ~ "ns",
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    is_significant = !is.na(p.value) & p.value <= 0.05,
    direction_label = case_when(
      diff > 0 ~ "Increase",
      diff < 0 ~ "Decrease",
      TRUE     ~ "No change"
    )
  ) %>%
  ungroup()

write_csv(all_stats, file.path(RESULTS_DIR, "matt_feature_statistics.csv"))
cat(sprintf("  Saved: matt_feature_statistics.csv (%d tests)\n", nrow(all_stats)))
cat(sprintf("  Significant: %d (p < 0.05)\n", sum(all_stats$is_significant, na.rm = TRUE)))

# ============================================================================
# ============================================================================
#  SECTION 2: DOT MATRIX PLOT — EXON FEATURES
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 2: DOT MATRIX PLOT\n")
cat(strrep("=", 70), "\n\n")

# Select key exon features for the dot matrix
# (adapt based on what columns are available in Matt output)
exon_features_priority <- c(
  # Exon properties
  "EXON_LENGTH", "MEDIAN_EXON_LENGTH", "EXON_GCC",
  # Flanking intron lengths
  "UPINTRON_MEDIANLENGTH", "DOINTRON_MEDIANLENGTH",
  "UPINTRON_LENGTH", "DOINTRON_LENGTH",
  # Splice site scores
  "MAXENTSCR_HSAMODEL_3SS", "MAXENTSCR_HSAMODEL_5SS",
  "MAXENTSCR_MMUMODEL_3SS", "MAXENTSCR_MMUMODEL_5SS",
  # Exon position
  "MEDIAN_EXON_NUMBER", "EXON_NUMBER",
  # GC content
  "UPINTRON_GCC", "DOINTRON_GCC",
  # Other features
  "TOTAL_EXON_NUMBER", "GENE_SIZE",
  "NMD_CODING_POTENTIAL"
)

# Check which features are actually present in the data
available_features <- unique(all_stats$feature)
features_to_plot <- intersect(exon_features_priority, available_features)

# If no priority features found, take most significant ones
if (length(features_to_plot) < 3) {
  features_to_plot <- all_stats %>%
    filter(is_significant) %>%
    group_by(feature) %>%
    summarise(min_p = min(p.value, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_p) %>%
    head(15) %>%
    pull(feature)
}

cat(sprintf("  Features selected for dot matrix: %d\n", length(features_to_plot)))

# Display name mapping
feature_display_names <- c(
  "EXON_LENGTH"              = "Exon Length",
  "MEDIAN_EXON_LENGTH"       = "Median Exon Length",
  "EXON_GCC"                 = "Exon GC Content",
  "UPINTRON_MEDIANLENGTH"    = "Upstream Intron Length",
  "DOINTRON_MEDIANLENGTH"    = "Downstream Intron Length",
  "UPINTRON_LENGTH"          = "Upstream Intron Length",
  "DOINTRON_LENGTH"          = "Downstream Intron Length",
  "MAXENTSCR_HSAMODEL_3SS"   = "MaxEnt 3'SS (human)",
  "MAXENTSCR_HSAMODEL_5SS"   = "MaxEnt 5'SS (human)",
  "MAXENTSCR_MMUMODEL_3SS"   = "MaxEnt 3'SS (mouse)",
  "MAXENTSCR_MMUMODEL_5SS"   = "MaxEnt 5'SS (mouse)",
  "MEDIAN_EXON_NUMBER"       = "Exon Number",
  "EXON_NUMBER"              = "Exon Number",
  "UPINTRON_GCC"             = "Upstream Intron GC",
  "DOINTRON_GCC"             = "Downstream Intron GC",
  "TOTAL_EXON_NUMBER"        = "Total Exon Count",
  "GENE_SIZE"                = "Gene Size",
  "NMD_CODING_POTENTIAL"     = "NMD Coding Potential"
)

# Prepare significant data for plotting
sig_data <- all_stats %>%
  filter(is_significant, feature %in% features_to_plot) %>%
  mutate(
    event_label = case_when(
      comparison == "up"   ~ "Included",
      comparison == "down" ~ "Skipped"
    ),
    treatment_label = case_when(
      treatment == "Myoblast_Taxol_vs_DMSO" ~ "Myoblast:\nTaxol vs DMSO",
      treatment == "Myotube_Taxol_vs_DMSO"  ~ "Myotube:\nTaxol vs DMSO",
      treatment == "Differentiation_DMSO"   ~ "Differentiation\n(DMSO)",
      treatment == "Differentiation_Taxol"  ~ "Differentiation\n(Taxol)"
    ),
    feature_display = ifelse(feature %in% names(feature_display_names),
                              feature_display_names[feature],
                              gsub("_", " ", feature)),
    dot_size = case_when(
      sig_level == "***" ~ 6,
      sig_level == "**"  ~ 4.5,
      sig_level == "*"   ~ 3,
      TRUE               ~ 0
    ),
    direction_color = direction_label
  )

if (nrow(sig_data) > 0) {

  # Order features by number of significant comparisons (most shared → top)
  feat_order <- sig_data %>%
    group_by(feature_display) %>%
    summarise(n_sig = n(), .groups = "drop") %>%
    arrange(n_sig) %>%
    pull(feature_display)

  sig_data$feature_display <- factor(sig_data$feature_display,
                                      levels = feat_order)

  # Treatment ordering
  treatment_order <- c("Myoblast:\nTaxol vs DMSO", "Myotube:\nTaxol vs DMSO",
                        "Differentiation\n(DMSO)", "Differentiation\n(Taxol)")
  sig_data$treatment_label <- factor(sig_data$treatment_label,
                                      levels = treatment_order)

  sig_data$event_label <- factor(sig_data$event_label,
                                  levels = c("Included", "Skipped"))

  n_features <- length(unique(sig_data$feature_display))
  n_treatments <- length(unique(sig_data$treatment_label))

  # Create interaction column for x-axis
  sig_data$x_group <- interaction(sig_data$treatment_label,
                                   sig_data$event_label, sep = "\n")

  # Build annotation bars (event type and comparison)
  n_treat <- length(treatment_order)
  event_bars <- data.frame(
    xmin  = c(0.5, n_treat + 0.5),
    xmax  = c(n_treat + 0.5, 2 * n_treat + 0.5),
    ymin  = n_features + 0.85,
    ymax  = n_features + 1.05,
    event = c("Included", "Skipped"),
    stringsAsFactors = FALSE
  )

  # Comparison annotation bars
  comp_bar_colors <- c(
    "Myoblast:\nTaxol vs DMSO" = "#D62839",
    "Myotube:\nTaxol vs DMSO"  = "#B2182B",
    "Differentiation\n(DMSO)"  = "#4BA3C3",
    "Differentiation\n(Taxol)" = "#2166AC"
  )

  comp_bars <- data.frame(
    xmin = seq(0.5, by = 1, length.out = 2 * n_treat),
    xmax = seq(1.5, by = 1, length.out = 2 * n_treat),
    ymin = n_features + 0.55,
    ymax = n_features + 0.75,
    comp = rep(treatment_order, 2),
    stringsAsFactors = FALSE
  )

  # ---- DOT MATRIX PLOT ----

  p_matrix <- ggplot(sig_data,
                      aes(x = x_group, y = feature_display)) +
    # Comparison colour bars
    geom_rect(
      data = comp_bars,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = comp),
      inherit.aes = FALSE, color = "white", linewidth = 0.5
    ) +
    # Event direction bars
    geom_rect(
      data = event_bars,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = event),
      inherit.aes = FALSE, color = "white", linewidth = 0.5
    ) +
    geom_text(
      data = event_bars,
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = event),
      inherit.aes = FALSE, size = 3.2, fontface = "bold", color = "white",
      family = FONT_FAMILY
    ) +
    # Dots: size = significance, colour = direction
    geom_point(
      aes(color = direction_color, size = dot_size),
      alpha = 0.85
    ) +
    scale_color_manual(
      name   = "Direction",
      values = direction_colors,
      breaks = c("Increase", "Decrease"),
      labels = c("Higher", "Lower")
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c(comp_bar_colors, event_colors),
      breaks = c("Included", "Skipped"),
      labels = c("Included", "Skipped")
    ) +
    scale_size_identity() +
    scale_x_discrete(expand = expansion(add = c(0.08, 0.08))) +
    scale_y_discrete(expand = expansion(add = c(0.1, 1.5))) +
    labs(
      title = "Exon Splicing Features — Taxol & Differentiation",
      subtitle = "Significant features (Wilcoxon p < 0.05, included/skipped vs unchanged)",
      x = NULL, y = NULL
    ) +
    theme_matrix() +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 22, size = 3),
        nrow = 1, title = NULL, order = 1
      ),
      color = guide_legend(
        override.aes = list(size = 3.5),
        nrow = 1, order = 2, title = "Direction"
      )
    ) +
    # Invisible points for size legend
    geom_point(
      data = data.frame(x = rep(Inf, 3), y = rep(Inf, 3),
                          size = c(3, 4.5, 6),
                          label = c("*", "**", "***")),
      aes(x = x, y = y, size = size),
      inherit.aes = FALSE, alpha = 0
    ) +
    scale_size_identity(
      name   = "p-value",
      breaks = c(3, 4.5, 6),
      labels = c("* (p < 0.05)", "** (p < 0.01)", "*** (p < 0.001)"),
      guide  = guide_legend(
        override.aes = list(alpha = 0.85, color = "grey40"),
        nrow = 1, order = 3
      )
    )

  plot_height <- max(5.0, n_features * 0.35 + 4.0)

  ggsave(file.path(PLOTS_DIR, "matt_dot_matrix_exons.pdf"),
         plot = p_matrix, width = 6.5, height = plot_height,
         device = SAVE_DEVICE, bg = "white")
  ggsave(file.path(PLOTS_DIR, "matt_dot_matrix_exons.png"),
         plot = p_matrix, width = 6.5, height = plot_height,
         dpi = 300, bg = "white")
  cat("  Saved: matt_dot_matrix_exons.pdf/.png\n")

} else {
  cat("  No significant features found — skipping dot matrix plot.\n")
}

# ============================================================================
# ============================================================================
#  SECTION 3: FEATURE BOXPLOTS (DETAILED)
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 3: FEATURE BOXPLOTS\n")
cat(strrep("=", 70), "\n\n")

# Select top features for detailed boxplots
top_features_for_boxplot <- all_stats %>%
  filter(is_significant) %>%
  group_by(feature) %>%
  summarise(min_p = min(p.value, na.rm = TRUE),
            n_sig = n(), .groups = "drop") %>%
  arrange(desc(n_sig), min_p) %>%
  head(6) %>%
  pull(feature)

cat(sprintf("  Plotting %d top features as boxplots\n", length(top_features_for_boxplot)))

if (length(top_features_for_boxplot) > 0) {

  boxplot_list <- list()

  for (feat in top_features_for_boxplot) {
    # Combine data from all comparisons
    feat_data <- bind_rows(
      lapply(names(matt_data), function(comp_name) {
        df <- matt_data[[comp_name]]
        if (feat %in% colnames(df) && "GROUP" %in% colnames(df)) {
          df %>%
            filter(GROUP %in% c("up", "down", "ndiff")) %>%
            dplyr::select(GROUP, value = !!sym(feat)) %>%
            mutate(treatment = comp_name) %>%
            filter(!is.na(value))
        }
      })
    ) %>%
      mutate(
        DATASET = factor(GROUP, levels = c("up", "ndiff", "down"),
                          labels = c("Included", "Unchanged", "Skipped")),
        treatment_label = case_when(
          treatment == "Myoblast_Taxol_vs_DMSO" ~ "Myoblast:\nTaxol vs DMSO",
          treatment == "Myotube_Taxol_vs_DMSO"  ~ "Myotube:\nTaxol vs DMSO",
          treatment == "Differentiation_DMSO"   ~ "Differentiation\n(DMSO)",
          treatment == "Differentiation_Taxol"  ~ "Differentiation\n(Taxol)"
        )
      )

    # Get display name
    disp_name <- ifelse(feat %in% names(feature_display_names),
                         feature_display_names[feat],
                         gsub("_", " ", feat))

    dataset_cols <- c("Included" = "#D62839", "Unchanged" = "grey70", "Skipped" = "#4BA3C3")

    p_box <- ggplot(feat_data, aes(x = DATASET, y = value, fill = DATASET)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3,
                   color = "grey30") +
      coord_cartesian(ylim = quantile(feat_data$value, c(0.01, 0.99), na.rm = TRUE)) +
      facet_wrap(~treatment_label, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = dataset_cols) +
      labs(title = disp_name, x = NULL, y = NULL) +
      theme_taxol(base_size = 9) +
      theme(
        legend.position = "none",
        axis.text.x     = element_text(angle = 30, hjust = 1, size = 7),
        strip.text       = element_text(size = 7)
      )

    boxplot_list[[feat]] <- p_box
  }

  # Combine boxplots
  if (length(boxplot_list) >= 2) {
    combined_box <- wrap_plots(boxplot_list, ncol = 2) +
      plot_annotation(
        title = "Key Exon Features: Included vs Skipped vs Unchanged",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                                     family = FONT_FAMILY)
        )
      )

    ggsave(file.path(PLOTS_DIR, "matt_feature_boxplots.pdf"), combined_box,
           width = 12, height = 4 * ceiling(length(boxplot_list) / 2),
           device = SAVE_DEVICE)
    cat("  Saved: matt_feature_boxplots.pdf\n")
  }
}

# ============================================================================
# ============================================================================
#  SECTION 4: TAXOL vs DIFFERENTIATION FEATURE COMPARISON
# ============================================================================
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SECTION 4: TAXOL vs DIFFERENTIATION FEATURE COMPARISON\n")
cat(strrep("=", 70), "\n\n")

# Compare which features distinguish taxol vs differentiation

if (all(c("Myoblast_Taxol_vs_DMSO", "Differentiation_DMSO") %in% names(matt_data))) {

  # Identify features significant in taxol but not differentiation (and vice versa)
  feature_comparison <- all_stats %>%
    filter(is_significant) %>%
    mutate(
      is_taxol = treatment %in% c("Myoblast_Taxol_vs_DMSO", "Myotube_Taxol_vs_DMSO"),
      is_diff  = treatment %in% c("Differentiation_DMSO", "Differentiation_Taxol")
    ) %>%
    group_by(feature, comparison) %>%
    summarise(
      taxol_sig = any(is_taxol),
      diff_sig  = any(is_diff),
      .groups = "drop"
    ) %>%
    mutate(
      category = case_when(
        taxol_sig & diff_sig  ~ "Both",
        taxol_sig & !diff_sig ~ "Taxol-specific",
        !taxol_sig & diff_sig ~ "Differentiation-specific",
        TRUE                  ~ "Neither"
      )
    )

  write_csv(feature_comparison,
            file.path(RESULTS_DIR, "matt_taxol_vs_diff_features.csv"))

  cat("  Feature comparison:\n")
  cat(sprintf("    Both: %d\n", sum(feature_comparison$category == "Both")))
  cat(sprintf("    Taxol-specific: %d\n",
              sum(feature_comparison$category == "Taxol-specific")))
  cat(sprintf("    Differentiation-specific: %d\n",
              sum(feature_comparison$category == "Differentiation-specific")))

  # Effect size comparison: log2FC heatmap
  effect_wide <- all_stats %>%
    filter(feature %in% features_to_plot) %>%
    mutate(
      label = paste0(treatment, "\n", comparison),
      feature_display = ifelse(feature %in% names(feature_display_names),
                                feature_display_names[feature],
                                gsub("_", " ", feature))
    ) %>%
    dplyr::select(feature_display, label, log2fc, is_significant) %>%
    distinct()

  if (nrow(effect_wide) > 0) {
    p_effect <- ggplot(effect_wide,
                        aes(x = label, y = feature_display, fill = log2fc)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = ifelse(is_significant, "*", "")),
                size = 4, color = "black", vjust = 0.8) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, name = expression(log[2]~"FC")
      ) +
      labs(
        title    = "Effect Size of Exon Features Across Comparisons",
        subtitle = "* = significant (Wilcoxon p < 0.05)",
        x = NULL, y = NULL
      ) +
      theme_taxol(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        legend.position = "right"
      )

    ggsave(file.path(PLOTS_DIR, "matt_effect_size_heatmap.pdf"), p_effect,
           width = 10, height = max(4, length(features_to_plot) * 0.4 + 2),
           device = SAVE_DEVICE)
    cat("  Saved: matt_effect_size_heatmap.pdf\n")
  }
}

# ============================================================================
# FINAL OUTPUT
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("MATT FEATURE ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Output:\n")
cat(sprintf("  Results: %s\n", RESULTS_DIR))
cat(sprintf("  Plots:   %s\n", PLOTS_DIR))
cat(sprintf("  Significant feature-comparison pairs: %d\n",
            sum(all_stats$is_significant, na.rm = TRUE)))
cat(sprintf("  Unique significant features: %d\n",
            length(unique(all_stats$feature[all_stats$is_significant]))))

cat("\nSession info:\n")
sessionInfo()
