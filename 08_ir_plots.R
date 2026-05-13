#!/usr/bin/env Rscript
# ============================================================================
# 08_ir_plots.R  [DEPRECATED — superseded by 10_final_figures.R]
# These intermediate IR plots were used during exploratory analysis.
# The final publication figures are now in 10_final_figures.R (Figures 1–4).
# Keeping this script for reference but plot saving is disabled.
# ============================================================================

cat("08_ir_plots.R is DEPRECATED. Use 10_final_figures.R instead.\n")
cat("Exiting without generating plots.\n")
# To re-enable for debugging, comment out the next line:
quit(save = "no")

cat(strrep("=", 70), "\n")
cat("08_ir_plots.R — IR-specific visualisations\n")
cat(strrep("=", 70), "\n\n")

setwd("/Users/andres.ortiz/Projects/CRG/taxol_rnaseq_c2c12")
RESULTS_DIR <- "results"
PLOTS_DIR   <- "plots"

theme_clean <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        strip.background = element_rect(fill = "grey95"))

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading data...\n")
diff_dmso  <- readRDS(file.path(RESULTS_DIR, "Differentiation_DMSO.rds"))
diff_taxol <- readRDS(file.path(RESULTS_DIR, "Differentiation_Taxol.rds"))
myo_taxol  <- readRDS(file.path(RESULTS_DIR, "Myoblast_Taxol_vs_DMSO.rds"))
tube_taxol <- readRDS(file.path(RESULTS_DIR, "Myotube_Taxol_vs_DMSO.rds"))

cat("Loading precomputed IR features...\n")
ir_features_all <- read.csv("results/ir_features_precomputed.csv", stringsAsFactors = FALSE)
cat(sprintf("  %d IR events with features\n", nrow(ir_features_all)))

speckle_genes <- scan("nuclear_speckle_associated_genes.txt", what = "character", quiet = TRUE)
speckle_upper <- toupper(speckle_genes)

# Thresholds
FDR_THRESH <- 0.05; DPSI_THRESH <- 0.1
is_sig <- function(fdr, dpsi) !is.na(fdr) & fdr <= FDR_THRESH & !is.na(dpsi) & abs(dpsi) >= DPSI_THRESH

# Split by event type
filter_ir <- function(df) df %>% filter(grepl("MmuINT", EVENT))
filter_ex <- function(df) df %>% filter(grepl("MmuEX", EVENT))

# ============================================================================
# 2. BUILD SCATTER DATA
# ============================================================================

cat("Building scatter data...\n")

# IR scatter
dd_ir <- filter_ir(diff_dmso) %>% select(EVENT, GENE, dpsi_dmso = deltapsi, FDR_dmso = FDR)
dt_ir <- filter_ir(diff_taxol) %>% select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)
ir_scatter <- inner_join(dd_ir, dt_ir, by = "EVENT") %>%
  mutate(sig_dmso = is_sig(FDR_dmso, dpsi_dmso),
         sig_taxol = is_sig(FDR_taxol, dpsi_taxol),
         quadrant = case_when(
           sig_dmso & sig_taxol & sign(dpsi_dmso) == sign(dpsi_taxol) ~ "Concordant",
           sig_dmso & sig_taxol & sign(dpsi_dmso) != sign(dpsi_taxol) ~ "Reversed",
           sig_dmso & !sig_taxol ~ "Blocked",
           !sig_dmso & sig_taxol ~ "Enabled",
           TRUE ~ "NS"),
         event_type = "IR")

# EX scatter
dd_ex <- filter_ex(diff_dmso) %>% select(EVENT, GENE, dpsi_dmso = deltapsi, FDR_dmso = FDR)
dt_ex <- filter_ex(diff_taxol) %>% select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)
ex_scatter <- inner_join(dd_ex, dt_ex, by = "EVENT") %>%
  mutate(sig_dmso = is_sig(FDR_dmso, dpsi_dmso),
         sig_taxol = is_sig(FDR_taxol, dpsi_taxol),
         quadrant = case_when(
           sig_dmso & sig_taxol & sign(dpsi_dmso) == sign(dpsi_taxol) ~ "Concordant",
           sig_dmso & sig_taxol & sign(dpsi_dmso) != sign(dpsi_taxol) ~ "Reversed",
           sig_dmso & !sig_taxol ~ "Blocked",
           !sig_dmso & sig_taxol ~ "Enabled",
           TRUE ~ "NS"),
         event_type = "EX")

# ============================================================================
# PLOT 1: IR-specific differentiation scatter with quadrant colouring
# ============================================================================

cat("Plot 1: IR differentiation scatter...\n")

quad_colours <- c("NS" = "grey85", "Concordant" = "#2166AC",
                  "Blocked" = "#B2182B", "Enabled" = "#4DAF4A", "Reversed" = "orange")

p1 <- ggplot(ir_scatter, aes(dpsi_dmso, dpsi_taxol, colour = quadrant)) +
  geom_point(data = ir_scatter %>% filter(quadrant == "NS"), size = 0.3, alpha = 0.15) +
  geom_point(data = ir_scatter %>% filter(quadrant != "NS"), size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = c(-DPSI_THRESH, DPSI_THRESH), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-DPSI_THRESH, DPSI_THRESH), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(values = quad_colours, name = "Quadrant",
                      breaks = c("Concordant", "Blocked", "Enabled", "Reversed")) +
  annotate("text", x = -0.5, y = 0.02, label = paste0("Blocked\n(n=", sum(ir_scatter$quadrant == "Blocked"), ")"),
           colour = "#B2182B", fontface = "bold", size = 3.5) +
  annotate("text", x = -0.5, y = -0.5, label = paste0("Concordant\n(n=", sum(ir_scatter$quadrant == "Concordant"), ")"),
           colour = "#2166AC", fontface = "bold", size = 3.5) +
  annotate("text", x = 0.02, y = -0.5, label = paste0("Enabled\n(n=", sum(ir_scatter$quadrant == "Enabled"), ")"),
           colour = "#4DAF4A", fontface = "bold", size = 3.5) +
  labs(title = "Intron Retention: Differentiation DMSO vs Taxol",
       x = expression(Delta*PIR ~ "(DMSO differentiation)"),
       y = expression(Delta*PIR ~ "(Taxol differentiation)")) +
  coord_fixed(xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7)) +
  theme_clean

ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_diff_scatter_quadrants.pdf"),
       p1, width = 7, height = 6.5)
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_diff_scatter_quadrants.png"),
       p1, width = 7, height = 6.5, dpi = 200)
cat("  Saved ir_diff_scatter_quadrants.pdf/png\n")

# ============================================================================
# PLOT 2: Side-by-side EX vs IR scatter
# ============================================================================

cat("Plot 2: EX vs IR scatter side-by-side...\n")

both_scatter <- bind_rows(ex_scatter, ir_scatter)

p2 <- ggplot(both_scatter, aes(dpsi_dmso, dpsi_taxol, colour = quadrant)) +
  geom_point(data = both_scatter %>% filter(quadrant == "NS"), size = 0.2, alpha = 0.1) +
  geom_point(data = both_scatter %>% filter(quadrant != "NS"), size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = c(-DPSI_THRESH, DPSI_THRESH), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-DPSI_THRESH, DPSI_THRESH), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(values = quad_colours, name = "Quadrant",
                      breaks = c("Concordant", "Blocked", "Enabled", "Reversed")) +
  facet_wrap(~factor(event_type, levels = c("EX", "IR")),
             labeller = labeller(event_type = c("EX" = "Cassette Exons (EX)", "IR" = "Intron Retention (IR)"))) +
  labs(title = "Differentiation Scatter: EX vs IR Events",
       x = expression(Delta*PSI ~ "or" ~ Delta*PIR ~ "(DMSO differentiation)"),
       y = expression(Delta*PSI ~ "or" ~ Delta*PIR ~ "(Taxol differentiation)")) +
  coord_fixed(xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7)) +
  theme_clean

ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_vs_ex_diff_scatter.pdf"),
       p2, width = 12, height = 6)
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_vs_ex_diff_scatter.png"),
       p2, width = 12, height = 6, dpi = 200)
cat("  Saved ir_vs_ex_diff_scatter.pdf/png\n")

# ============================================================================
# PLOT 3: Speckle enrichment — IR vs EX bar plot
# ============================================================================

cat("Plot 3: Speckle enrichment IR vs EX...\n")

# Compute enrichment for each event_type x quadrant
enrich_data <- list()
for (etype in c("EX", "IR")) {
  df <- if (etype == "EX") ex_scatter else ir_scatter
  df$GENE_upper <- toupper(df$GENE)
  bg <- unique(df$GENE_upper)
  bg_sp <- sum(bg %in% speckle_upper)
  bg_pct <- 100 * bg_sp / length(bg)

  enrich_data[[paste0(etype, "_Background")]] <- data.frame(
    event_type = etype, quadrant = "Background",
    n_genes = length(bg), n_speckle = bg_sp,
    pct_speckle = bg_pct, stringsAsFactors = FALSE)

  for (q in c("Concordant", "Blocked", "Enabled")) {
    qg <- unique(df$GENE_upper[df$quadrant == q])
    qsp <- sum(qg %in% speckle_upper)
    qpct <- 100 * qsp / max(length(qg), 1)

    # Fisher test
    n_bg <- length(bg); n_sp_bg <- bg_sp
    n_q <- length(qg); n_sp_q <- qsp
    cont <- matrix(c(n_sp_q, n_q - n_sp_q, n_sp_bg - n_sp_q,
                      max(0, n_bg - n_q - n_sp_bg + n_sp_q)), nrow = 2)
    ft <- tryCatch(fisher.test(cont), error = function(e) list(estimate = NA, p.value = 1))

    enrich_data[[paste0(etype, "_", q)]] <- data.frame(
      event_type = etype, quadrant = q,
      n_genes = n_q, n_speckle = qsp,
      pct_speckle = qpct, or = unname(ft$estimate), p = ft$p.value,
      stringsAsFactors = FALSE)
  }
}

enrich_df <- bind_rows(enrich_data) %>%
  mutate(quadrant = factor(quadrant, levels = c("Background", "Concordant", "Blocked", "Enabled")),
         event_type = factor(event_type, levels = c("EX", "IR")),
         sig_label = ifelse(!is.na(p) & p < 0.001, "***",
                     ifelse(!is.na(p) & p < 0.01, "**",
                     ifelse(!is.na(p) & p < 0.05, "*", "ns"))))

quad_fill <- c("Background" = "grey70", "Concordant" = "#2166AC",
               "Blocked" = "#B2182B", "Enabled" = "#4DAF4A")

p3 <- ggplot(enrich_df, aes(x = quadrant, y = pct_speckle, fill = quadrant)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(sig_label) | quadrant == "Background", "",
                                paste0(sig_label, "\nOR=", ifelse(is.na(or), "", sprintf("%.1f", or))))),
            vjust = -0.3, size = 3) +
  geom_text(aes(label = paste0("n=", n_genes)), vjust = 1.5, size = 2.8, colour = "white", fontface = "bold") +
  facet_wrap(~event_type, labeller = labeller(event_type = c("EX" = "Cassette Exons", "IR" = "Intron Retention"))) +
  scale_fill_manual(values = quad_fill, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Nuclear Speckle Gene Enrichment: EX vs IR",
       subtitle = "Blocked IR shows strongest enrichment (OR=4.33); Blocked EX is NS",
       x = "Quadrant", y = "% Speckle-associated genes") +
  theme_clean +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_vs_ex_speckle_enrichment.pdf"),
       p3, width = 9, height = 6)
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_vs_ex_speckle_enrichment.png"),
       p3, width = 9, height = 6, dpi = 200)
cat("  Saved ir_vs_ex_speckle_enrichment.pdf/png\n")

# ============================================================================
# PLOT 4: Intron length and GC distributions per IR quadrant
# ============================================================================

cat("Plot 4: Intron features per quadrant...\n")

ir_features <- ir_features_all

ir_feat_plot <- ir_scatter %>%
  left_join(ir_features, by = "EVENT") %>%
  filter(quadrant != "Reversed") %>%
  mutate(quadrant = factor(quadrant, levels = c("NS", "Concordant", "Blocked", "Enabled")))

# Intron length
p4a <- ggplot(ir_feat_plot %>% filter(quadrant != "NS"),
              aes(x = quadrant, y = log10(INTRON_LENGTH), fill = quadrant)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, outlier.alpha = 0.3) +
  geom_hline(yintercept = log10(median(ir_feat_plot$INTRON_LENGTH[ir_feat_plot$quadrant == "NS"], na.rm = TRUE)),
             linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.55, y = log10(median(ir_feat_plot$INTRON_LENGTH[ir_feat_plot$quadrant == "NS"], na.rm = TRUE)) + 0.05,
           label = "NS median", size = 3, colour = "grey40", hjust = 0) +
  scale_fill_manual(values = c("Concordant" = "#2166AC", "Blocked" = "#B2182B", "Enabled" = "#4DAF4A"), guide = "none") +
  labs(title = "IR: Intron Length by Quadrant",
       subtitle = "Blocked introns are shortest (median 752 nt)",
       x = "Quadrant", y = expression(log[10] ~ "(Intron length, nt)")) +
  theme_clean

# Intron GC
p4b <- ggplot(ir_feat_plot %>% filter(quadrant != "NS"),
              aes(x = quadrant, y = INTRON_GCC, fill = quadrant)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, outlier.alpha = 0.3) +
  geom_hline(yintercept = median(ir_feat_plot$INTRON_GCC[ir_feat_plot$quadrant == "NS"], na.rm = TRUE),
             linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.55, y = median(ir_feat_plot$INTRON_GCC[ir_feat_plot$quadrant == "NS"], na.rm = TRUE) + 0.003,
           label = "NS median", size = 3, colour = "grey40", hjust = 0) +
  scale_fill_manual(values = c("Concordant" = "#2166AC", "Blocked" = "#B2182B", "Enabled" = "#4DAF4A"), guide = "none") +
  labs(title = "IR: Intron GC Content by Quadrant",
       subtitle = "All regulated introns are GC-rich; blocked most extreme",
       x = "Quadrant", y = "Intron GC content") +
  theme_clean

p4 <- p4a | p4b
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_intron_features_by_quadrant.pdf"),
       p4, width = 10, height = 5)
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_intron_features_by_quadrant.png"),
       p4, width = 10, height = 5, dpi = 200)
cat("  Saved ir_intron_features_by_quadrant.pdf/png\n")

# ============================================================================
# PLOT 5: Blocked IR direction breakdown
# ============================================================================

cat("Plot 5: Blocked IR direction...\n")

blocked_dir <- ir_scatter %>%
  filter(quadrant == "Blocked") %>%
  mutate(direction = ifelse(dpsi_dmso > 0, "Retention\nincreased", "Retention\ndecreased"),
         is_speckle = toupper(GENE) %in% speckle_upper)

dir_counts <- blocked_dir %>%
  group_by(direction, is_speckle) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(direction) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup()

p5a <- ggplot(blocked_dir, aes(x = direction, fill = direction)) +
  geom_bar(width = 0.6, colour = "grey30", linewidth = 0.3) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Retention\nincreased" = "#E08214", "Retention\ndecreased" = "#542788"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Blocked IR Events: Direction in DMSO Differentiation",
       subtitle = "90% are retention-decreased (introns that should be removed stay retained)",
       x = "", y = "Number of events") +
  theme_clean

p5b <- ggplot(blocked_dir, aes(x = direction, fill = is_speckle)) +
  geom_bar(position = "fill", width = 0.6, colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#D6604D"),
                    name = "Speckle gene", labels = c("No", "Yes")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Speckle Enrichment by Direction",
       subtitle = "Retention-decreased events in speckle genes: 34%",
       x = "", y = "Proportion") +
  theme_clean

p5 <- p5a | p5b
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_blocked_direction.pdf"),
       p5, width = 10, height = 5)
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_blocked_direction.png"),
       p5, width = 10, height = 5, dpi = 200)
cat("  Saved ir_blocked_direction.pdf/png\n")

# ============================================================================
# PLOT 6: Block rate comparison (EX vs IR bar chart)
# ============================================================================

cat("Plot 6: Block rate comparison...\n")

block_rates <- data.frame(
  event_type = c("EX", "IR"),
  sig_dmso = c(sum(is_sig(filter_ex(diff_dmso)$FDR, filter_ex(diff_dmso)$deltapsi)),
                sum(is_sig(filter_ir(diff_dmso)$FDR, filter_ir(diff_dmso)$deltapsi))),
  blocked = c(sum(ex_scatter$quadrant == "Blocked"),
              sum(ir_scatter$quadrant == "Blocked"))
) %>%
  mutate(block_rate = 100 * blocked / sig_dmso,
         event_type = factor(event_type, levels = c("EX", "IR")))

p6 <- ggplot(block_rates, aes(x = event_type, y = block_rate, fill = event_type)) +
  geom_col(width = 0.5, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d/%d)", block_rate, blocked, sig_dmso)),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("EX" = "#4393C3", "IR" = "#D6604D"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)), limits = c(0, 65)) +
  labs(title = "Taxol Block Rate: Cassette Exons vs Intron Retention",
       subtitle = "IR events are 2.5× more likely to be blocked during differentiation",
       x = "Event type", y = "% of DMSO-significant events blocked by Taxol") +
  theme_clean

ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_vs_ex_block_rate.pdf"),
       p6, width = 5.5, height = 5)
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_vs_ex_block_rate.png"),
       p6, width = 5.5, height = 5, dpi = 200)
cat("  Saved ir_vs_ex_block_rate.pdf/png\n")

# ============================================================================
# PLOT 7: Combined summary figure
# ============================================================================

cat("Plot 7: Combined summary figure...\n")

p_combined <- (p6 | p3) / (p4) / (p5) +
  plot_annotation(
    title = "Intron Retention Analysis: Taxol Blocks Speckle-Dependent Intron Removal",
    tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_analysis_combined_figure.pdf"),
       p_combined, width = 14, height = 16)
ggsave(file.path(PLOTS_DIR, "09_speckle_analysis", "ir_analysis_combined_figure.png"),
       p_combined, width = 14, height = 16, dpi = 200)
cat("  Saved ir_analysis_combined_figure.pdf/png\n")

cat("\n", strrep("=", 70), "\n")
cat("All IR plots generated successfully.\n")
cat(strrep("=", 70), "\n")
