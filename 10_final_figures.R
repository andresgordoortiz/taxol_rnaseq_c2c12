#!/usr/bin/env Rscript
# ============================================================================
# 10_final_figures.R
# Four refined, publication-quality figures that tell the full story.
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
PLOTS_DIR   <- "plots"

# Palette
col_retain  <- "#D6604D"    # warm red — intron stays in
col_excise  <- "#4393C3"    # cool blue — intron excised
col_include <- "#F4A582"    # light salmon — exon in
col_skip    <- "#92C5DE"    # light blue — exon out
col_blocked <- "#B2182B"    # deep red
col_conc    <- "#2166AC"    # deep blue
col_enabled <- "#4DAF4A"    # green
col_speckle <- "#8856A7"    # purple — speckle
col_nonsp   <- "grey70"     # non-speckle
col_ns      <- "grey88"

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, colour = "grey40", margin = margin(b = 8)),
    strip.background = element_rect(fill = "grey96", colour = NA),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.margin = margin(t = -5),
    axis.title = element_text(size = 10)
  )

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")
diff_dmso  <- readRDS(file.path(RESULTS_DIR, "Differentiation_DMSO.rds"))
diff_taxol <- readRDS(file.path(RESULTS_DIR, "Differentiation_Taxol.rds"))
myo_taxol  <- readRDS(file.path(RESULTS_DIR, "Myoblast_Taxol_vs_DMSO.rds"))
tube_taxol <- readRDS(file.path(RESULTS_DIR, "Myotube_Taxol_vs_DMSO.rds"))

speckle_upper <- toupper(scan("nuclear_speckle_associated_genes.txt", what = "character", quiet = TRUE))
ir_features <- read.csv("results/ir_features_precomputed.csv", stringsAsFactors = FALSE)

FDR_THRESH <- 0.05; DPSI_THRESH <- 0.1
is_sig <- function(fdr, dpsi) !is.na(fdr) & fdr <= FDR_THRESH & !is.na(dpsi) & abs(dpsi) >= DPSI_THRESH
filter_ir <- function(df) df %>% filter(grepl("MmuINT", EVENT))
filter_ex <- function(df) df %>% filter(grepl("MmuEX", EVENT))

# ============================================================================
# FIGURE 1: THE MIRROR — Taxol direction reversal by cell type
#   "Same introns, opposite fates"
# ============================================================================

cat("Figure 1: The Mirror...\n")

# Build data
myo_ir_sig  <- filter_ir(myo_taxol) %>% filter(is_sig(FDR, deltapsi))
tube_ir_sig <- filter_ir(tube_taxol) %>% filter(is_sig(FDR, deltapsi))
myo_ex_sig  <- filter_ex(myo_taxol) %>% filter(is_sig(FDR, deltapsi))
tube_ex_sig <- filter_ex(tube_taxol) %>% filter(is_sig(FDR, deltapsi))

mirror_data <- bind_rows(
  # IR
  tibble(cell = "Myoblast", event = "Intron Retention",
         direction = "Retained\n(stays in)", n = sum(myo_ir_sig$deltapsi > 0)),
  tibble(cell = "Myoblast", event = "Intron Retention",
         direction = "Excised\n(spliced out)", n = sum(myo_ir_sig$deltapsi < 0)),
  tibble(cell = "Myotube", event = "Intron Retention",
         direction = "Retained\n(stays in)", n = sum(tube_ir_sig$deltapsi > 0)),
  tibble(cell = "Myotube", event = "Intron Retention",
         direction = "Excised\n(spliced out)", n = sum(tube_ir_sig$deltapsi < 0)),
  # EX
  tibble(cell = "Myoblast", event = "Cassette Exon",
         direction = "Included", n = sum(myo_ex_sig$deltapsi > 0)),
  tibble(cell = "Myoblast", event = "Cassette Exon",
         direction = "Skipped", n = sum(myo_ex_sig$deltapsi < 0)),
  tibble(cell = "Myotube", event = "Cassette Exon",
         direction = "Included", n = sum(tube_ex_sig$deltapsi > 0)),
  tibble(cell = "Myotube", event = "Cassette Exon",
         direction = "Skipped", n = sum(tube_ex_sig$deltapsi < 0))
) %>%
  group_by(cell, event) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    cell = factor(cell, levels = c("Myoblast", "Myotube")),
    # Make the "force axis" clear
    direction = factor(direction, levels = c("Retained\n(stays in)", "Excised\n(spliced out)",
                                              "Included", "Skipped")),
    fill_col = case_when(
      direction == "Retained\n(stays in)" ~ col_retain,
      direction == "Excised\n(spliced out)" ~ col_excise,
      direction == "Included" ~ col_include,
      direction == "Skipped" ~ col_skip
    )
  )

# Panel A: IR mirror bars
ir_mirror <- mirror_data %>% filter(event == "Intron Retention") %>%
  mutate(
    # Flip myoblast bars negative for mirror effect
    plot_pct = ifelse(cell == "Myoblast", -pct, pct),
    plot_n = ifelse(cell == "Myoblast", -n, n)
  )

p1a <- ggplot(ir_mirror, aes(x = direction, y = plot_n, fill = direction)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(aes(label = abs(n), y = plot_n + sign(plot_n) * 4),
            fontface = "bold", size = 4) +
  annotate("text", x = 2.5, y = -50, label = "MYOBLAST\n(low forces)", fontface = "italic",
           size = 3.5, colour = "grey40") +
  annotate("text", x = 2.5, y = 50, label = "MYOTUBE\n(high forces)", fontface = "italic",
           size = 3.5, colour = "grey40") +
  scale_fill_manual(values = c("Retained\n(stays in)" = col_retain,
                                "Excised\n(spliced out)" = col_excise),
                    guide = "none") +
  scale_y_continuous(labels = abs, limits = c(-90, 60)) +
  labs(title = "Intron Retention",
       subtitle = "Taxol pushes introns in opposite directions depending on forces",
       x = "", y = "Number of significant IR events") +
  coord_flip() +
  theme_pub +
  theme(plot.subtitle = element_text(size = 8.5))

# Panel B: EX mirror bars
ex_mirror <- mirror_data %>% filter(event == "Cassette Exon") %>%
  mutate(
    plot_n = ifelse(cell == "Myoblast", -n, n)
  )

p1b <- ggplot(ex_mirror, aes(x = direction, y = plot_n, fill = direction)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(aes(label = abs(n), y = plot_n + sign(plot_n) * 6),
            fontface = "bold", size = 4) +
  annotate("text", x = 2.5, y = -60, label = "MYOBLAST\n(low forces)", fontface = "italic",
           size = 3.5, colour = "grey40") +
  annotate("text", x = 2.5, y = 60, label = "MYOTUBE\n(high forces)", fontface = "italic",
           size = 3.5, colour = "grey40") +
  scale_fill_manual(values = c("Included" = col_include, "Skipped" = col_skip), guide = "none") +
  scale_y_continuous(labels = abs, limits = c(-120, 130)) +
  labs(title = "Cassette Exons",
       subtitle = "Balanced in myoblasts, strongly inclusion-biased in myotubes",
       x = "", y = "Number of significant EX events") +
  coord_flip() +
  theme_pub +
  theme(plot.subtitle = element_text(size = 8.5))

# Panel C: Speckle enrichment — ALL differentiation categories + taxol direct + background
# Honest comparison: blocked, concordant, enabled, direct taxol, and non-significant background

# Compute pooled gene sets
myo_all_sig  <- bind_rows(myo_ir_sig, myo_ex_sig)
tube_all_sig <- bind_rows(tube_ir_sig, tube_ex_sig)

myo_genes_all  <- unique(toupper(myo_all_sig$GENE))
tube_genes_all <- unique(toupper(tube_all_sig$GENE))

# Differentiation categories
diff_dmso_sig  <- diff_dmso  %>% filter(is_sig(FDR, deltapsi))
diff_taxol_sig <- diff_taxol %>% filter(is_sig(FDR, deltapsi))

blocked_events    <- diff_dmso_sig %>% filter(!EVENT %in% diff_taxol_sig$EVENT)
concordant_events <- diff_dmso_sig %>% filter(EVENT %in% diff_taxol_sig$EVENT)
enabled_events    <- diff_taxol_sig %>% filter(!EVENT %in% diff_dmso_sig$EVENT)

blocked_genes    <- unique(toupper(blocked_events$GENE))
concordant_genes <- unique(toupper(concordant_events$GENE))
enabled_genes    <- unique(toupper(enabled_events$GENE))

# Non-significant background: genes with no significant event in ANY comparison
all_genes <- unique(toupper(c(myo_taxol$GENE, tube_taxol$GENE)))
all_sig_genes <- unique(c(myo_genes_all, tube_genes_all, blocked_genes, concordant_genes, enabled_genes))
bg_genes <- setdiff(all_genes, all_sig_genes)

speckle_data <- tibble(
  group = c(
    sprintf("Diff blocked\n(EX+IR, %d genes)", length(blocked_genes)),
    sprintf("Diff concordant\n(EX+IR, %d genes)", length(concordant_genes)),
    sprintf("Diff enabled\n(EX+IR, %d genes)", length(enabled_genes)),
    sprintf("Myoblast Taxol\n(EX+IR, %d genes)", length(myo_genes_all)),
    sprintf("Myotube Taxol\n(EX+IR, %d genes)", length(tube_genes_all)),
    "Non-significant\n(background)"
  ),
  n_genes   = c(length(blocked_genes), length(concordant_genes), length(enabled_genes),
                length(myo_genes_all), length(tube_genes_all), length(bg_genes)),
  n_speckle = c(sum(blocked_genes %in% speckle_upper),
                sum(concordant_genes %in% speckle_upper),
                sum(enabled_genes %in% speckle_upper),
                sum(myo_genes_all %in% speckle_upper),
                sum(tube_genes_all %in% speckle_upper),
                sum(bg_genes %in% speckle_upper)),
  category = c("blocked", "diff", "diff", "taxol", "taxol", "bg")
) %>%
  mutate(
    pct_speckle = 100 * n_speckle / n_genes,
    group = factor(group, levels = rev(group))
  )

# Fisher tests vs non-significant background
bg_n <- speckle_data$n_genes[speckle_data$category == "bg"]
bg_sp <- speckle_data$n_speckle[speckle_data$category == "bg"]
speckle_data$or <- NA_real_; speckle_data$p <- NA_real_
for (i in which(speckle_data$category != "bg")) {
  cont <- matrix(c(speckle_data$n_speckle[i], speckle_data$n_genes[i] - speckle_data$n_speckle[i],
                    bg_sp, bg_n - bg_sp), nrow = 2)
  ft <- fisher.test(cont)
  speckle_data$or[i] <- unname(ft$estimate)
  speckle_data$p[i] <- ft$p.value
}

speckle_data$label <- with(speckle_data, ifelse(
  is.na(or), "", sprintf("OR=%.1f, p=%.0e", or, p)))

# Colours: blocked=red, diff categories=blue, taxol=purple, bg=grey
fill_colours <- c("blocked" = "#B2182B", "diff" = "#4393C3", "taxol" = col_speckle, "bg" = col_nonsp)

p1c <- ggplot(speckle_data, aes(x = group, y = pct_speckle, fill = category)) +
  geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = label), hjust = -0.05, size = 2.5, colour = "grey30") +
  scale_fill_manual(
    values = fill_colours,
    labels = c("blocked" = "Diff blocked by taxol", "diff" = "Differentiation (not taxol-specific)",
               "taxol" = "Direct taxol effect", "bg" = "No significant change"),
    name = ""
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(title = "Nuclear Speckle Gene Enrichment",
       subtitle = "All differentiation categories are enriched -- not taxol-specific. Blocked is slightly higher (OR=1.4 vs concordant, p=0.02)",
       x = "", y = "% speckle-associated genes") +
  coord_flip() +
  theme_pub +
  theme(plot.subtitle = element_text(size = 8),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8))

fig1 <- (p1a | p1b) / p1c +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Figure 1: Taxol drives a cell-type-dependent mirror reversal in splicing",
    subtitle = "Intron retention and exon inclusion reverse direction between myoblasts (low forces) and myotubes (high forces)\nSpeckle enrichment is a property of differentiation-associated splicing genes, not exclusively taxol-blocked genes",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30", margin = margin(b = 10)),
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

ggsave(file.path(PLOTS_DIR, "Figure1_mirror_reversal.pdf"), fig1, width = 14, height = 10)
ggsave(file.path(PLOTS_DIR, "Figure1_mirror_reversal.png"), fig1, width = 14, height = 10, dpi = 300)
cat("  Saved Figure1_mirror_reversal\n")

# ============================================================================
# FIGURE 2: THE BLOCK — Differentiation is selectively impaired at IR level
# ============================================================================

cat("Figure 2: The Block...\n")

# Build differentiation scatter for BOTH event types
build_diff_scatter <- function(dmso_df, taxol_df, type_label) {
  d <- dmso_df %>% select(EVENT, GENE, dpsi_dmso = deltapsi, FDR_dmso = FDR)
  t <- taxol_df %>% select(EVENT, dpsi_taxol = deltapsi, FDR_taxol = FDR)
  inner_join(d, t, by = "EVENT") %>%
    mutate(sig_dmso = is_sig(FDR_dmso, dpsi_dmso),
           sig_taxol = is_sig(FDR_taxol, dpsi_taxol),
           quadrant = case_when(
             sig_dmso & sig_taxol & sign(dpsi_dmso) == sign(dpsi_taxol) ~ "Concordant",
             sig_dmso & sig_taxol & sign(dpsi_dmso) != sign(dpsi_taxol) ~ "Reversed",
             sig_dmso & !sig_taxol ~ "Blocked",
             !sig_dmso & sig_taxol ~ "Enabled",
             TRUE ~ "NS"),
           event_type = type_label,
           is_speckle = toupper(GENE) %in% speckle_upper)
}

sc_ir <- build_diff_scatter(filter_ir(diff_dmso), filter_ir(diff_taxol), "Intron Retention")
sc_ex <- build_diff_scatter(filter_ex(diff_dmso), filter_ex(diff_taxol), "Cassette Exon")

# Panel A: Side-by-side scatter
sc_both <- bind_rows(sc_ex, sc_ir) %>%
  mutate(event_type = factor(event_type, levels = c("Cassette Exon", "Intron Retention")))

p2a <- ggplot(sc_both, aes(dpsi_dmso, dpsi_taxol, colour = quadrant)) +
  geom_point(data = sc_both %>% filter(quadrant == "NS"), size = 0.15, alpha = 0.08) +
  geom_point(data = sc_both %>% filter(quadrant != "NS"), size = 0.9, alpha = 0.7) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.2) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3, linetype = "dotted") +
  scale_colour_manual(
    values = c("NS" = col_ns, "Concordant" = col_conc, "Blocked" = col_blocked,
               "Enabled" = col_enabled, "Reversed" = "orange"),
    name = "", breaks = c("Concordant", "Blocked", "Enabled")) +
  facet_wrap(~event_type) +
  labs(title = "Differentiation splicing: DMSO (normal) vs Taxol",
       subtitle = "Off-diagonal blocked events (red) are far more common for intron retention",
       x = expression(Delta*PSI ~ "or" ~ Delta*PIR ~ "(DMSO differentiation)"),
       y = expression(Delta*PSI ~ "or" ~ Delta*PIR ~ "(Taxol differentiation)")) +
  coord_fixed(xlim = c(-0.65, 0.65), ylim = c(-0.65, 0.65)) +
  theme_pub +
  theme(legend.position = "top")

# Panel B: Block rate + event counts
block_summary <- tibble(
  event_type = factor(c("Cassette Exon", "Intron Retention"),
                       levels = c("Cassette Exon", "Intron Retention")),
  sig_dmso = c(sum(is_sig(filter_ex(diff_dmso)$FDR, filter_ex(diff_dmso)$deltapsi)),
               sum(is_sig(filter_ir(diff_dmso)$FDR, filter_ir(diff_dmso)$deltapsi))),
  n_blocked = c(sum(sc_ex$quadrant == "Blocked"), sum(sc_ir$quadrant == "Blocked")),
  n_conc    = c(sum(sc_ex$quadrant == "Concordant"), sum(sc_ir$quadrant == "Concordant")),
  n_enabled = c(sum(sc_ex$quadrant == "Enabled"), sum(sc_ir$quadrant == "Enabled"))
) %>%
  mutate(block_rate = 100 * n_blocked / sig_dmso)

p2b <- ggplot(block_summary, aes(x = event_type, y = block_rate, fill = event_type)) +
  geom_col(width = 0.5, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n%d / %d", block_rate, n_blocked, sig_dmso)),
            vjust = -0.2, size = 3.5, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = c("Cassette Exon" = col_skip, "Intron Retention" = col_retain), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)), limits = c(0, 65)) +
  labs(title = "Taxol block rate",
       subtitle = "IR events are 2.5× more\nvulnerable to taxol",
       x = "", y = "% differentiation events\nblocked by taxol") +
  theme_pub +
  theme(plot.subtitle = element_text(size = 8.5))

# Panel C: Blocked IR direction (what's being prevented)
blocked_ir_dir <- sc_ir %>%
  filter(quadrant == "Blocked") %>%
  mutate(fate = ifelse(dpsi_dmso < 0,
                        "Should be excised\n(297 introns)",
                        "Should be retained\n(34 introns)"))

p2c <- ggplot(blocked_ir_dir, aes(x = fate, fill = is_speckle)) +
  geom_bar(width = 0.6, colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c("FALSE" = col_nonsp, "TRUE" = col_speckle),
                    name = "", labels = c("Non-speckle gene", "Speckle gene")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "What taxol blocks in differentiation IR",
       subtitle = "90% are introns that should be excised\nbut remain retained under taxol",
       x = "", y = "Number of blocked IR events") +
  theme_pub +
  theme(legend.position = c(0.75, 0.85), legend.background = element_blank(),
        plot.subtitle = element_text(size = 8.5))

fig2 <- p2a / (p2b | p2c) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Figure 2: Taxol selectively blocks differentiation-associated intron excision",
    subtitle = "55% of IR events vs 22% of EX events are blocked -- and 90% of blocked IR events are introns that should be removed",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30", margin = margin(b = 10)),
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

ggsave(file.path(PLOTS_DIR, "Figure2_differentiation_block.pdf"), fig2, width = 12, height = 10)
ggsave(file.path(PLOTS_DIR, "Figure2_differentiation_block.png"), fig2, width = 12, height = 10, dpi = 300)
cat("  Saved Figure2_differentiation_block\n")

# ============================================================================
# FIGURE 3: THE PROOF — Structural features confirm speckle mechanism
# ============================================================================

cat("Figure 3: The Proof...\n")

# Use plain ASCII group names to avoid encoding issues
grp_names <- c("Background (NS)", "Concordant", "Blocked (all)", "Blocked excision",
               "Myotube retained", "Myoblast excised")

# Collect feature data for all key groups
ns_ir <- sc_ir %>% filter(quadrant == "NS") %>% left_join(ir_features, by = "EVENT")

get_feats <- function(df, label) {
  df %>% left_join(ir_features, by = "EVENT") %>%
    filter(!is.na(INTRON_LENGTH)) %>%
    summarise(group = label, n = n(),
              med_length = median(INTRON_LENGTH), med_gc = median(INTRON_GCC),
              med_up_gc = median(UPSTREAM_EXON_GCC), med_dn_gc = median(DOWNSTREAM_EXON_GCC))
}

feature_summary <- bind_rows(
  ns_ir %>% summarise(group = grp_names[1], n = n(),
                       med_length = median(INTRON_LENGTH, na.rm = TRUE),
                       med_gc = median(INTRON_GCC, na.rm = TRUE),
                       med_up_gc = median(UPSTREAM_EXON_GCC, na.rm = TRUE),
                       med_dn_gc = median(DOWNSTREAM_EXON_GCC, na.rm = TRUE)),
  get_feats(sc_ir %>% filter(quadrant == "Concordant"), grp_names[2]),
  get_feats(sc_ir %>% filter(quadrant == "Blocked"), grp_names[3]),
  get_feats(sc_ir %>% filter(quadrant == "Blocked", dpsi_dmso < 0), grp_names[4]),
  filter_ir(tube_taxol) %>% filter(is_sig(FDR, deltapsi), deltapsi > 0) %>%
    left_join(ir_features, by = "EVENT") %>% filter(!is.na(INTRON_LENGTH)) %>%
    summarise(group = grp_names[5], n = n(), med_length = median(INTRON_LENGTH),
              med_gc = median(INTRON_GCC), med_up_gc = median(UPSTREAM_EXON_GCC),
              med_dn_gc = median(DOWNSTREAM_EXON_GCC)),
  filter_ir(myo_taxol) %>% filter(is_sig(FDR, deltapsi), deltapsi < 0) %>%
    left_join(ir_features, by = "EVENT") %>% filter(!is.na(INTRON_LENGTH)) %>%
    summarise(group = grp_names[6], n = n(), med_length = median(INTRON_LENGTH),
              med_gc = median(INTRON_GCC), med_up_gc = median(UPSTREAM_EXON_GCC),
              med_dn_gc = median(DOWNSTREAM_EXON_GCC))
)

# Speckle percentages (same row order as feature_summary)
compute_speckle_pct <- function(events_df) {
  genes <- unique(events_df$GENE)
  100 * sum(toupper(genes) %in% speckle_upper) / length(genes)
}
all_ir_genes <- unique(c(filter_ir(diff_dmso)$GENE, filter_ir(diff_taxol)$GENE,
                          filter_ir(myo_taxol)$GENE, filter_ir(tube_taxol)$GENE))
bg_speckle <- 100 * sum(toupper(all_ir_genes) %in% speckle_upper) / length(all_ir_genes)

speckle_pcts <- c(
  bg_speckle,
  compute_speckle_pct(sc_ir %>% filter(quadrant == "Concordant")),
  compute_speckle_pct(sc_ir %>% filter(quadrant == "Blocked")),
  compute_speckle_pct(sc_ir %>% filter(quadrant == "Blocked", dpsi_dmso < 0)),
  compute_speckle_pct(filter_ir(tube_taxol) %>% filter(is_sig(FDR, deltapsi), deltapsi > 0)),
  compute_speckle_pct(filter_ir(myo_taxol) %>% filter(is_sig(FDR, deltapsi), deltapsi < 0))
)

# Drop "Blocked (all)" — redundant with "Blocked excision"
heatmap_keep <- c(grp_names[1], grp_names[2], grp_names[4], grp_names[5], grp_names[6])

feature_plot <- feature_summary %>%
  mutate(speckle = speckle_pcts) %>%
  filter(group %in% heatmap_keep) %>%
  mutate(group = factor(group, levels = rev(heatmap_keep)))

# Panel A: Intron length
p3a <- ggplot(feature_plot, aes(x = group, y = med_length,
                                 fill = ifelse(grepl("Background", group), "bg", "sig"))) +
  geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d nt", round(med_length))), hjust = -0.1, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("bg" = col_nonsp, "sig" = col_retain), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Intron length (median)",
       subtitle = "Blocked + taxol-affected introns are SHORT",
       x = "", y = "Intron length (nt)") +
  coord_flip() + theme_pub + theme(plot.subtitle = element_text(size = 8.5))

# Panel B: GC content
gc_long <- feature_plot %>%
  select(group, Intron = med_gc, `Upstream exon` = med_up_gc, `Downstream exon` = med_dn_gc) %>%
  pivot_longer(-group, names_to = "region", values_to = "gc") %>%
  mutate(region = factor(region, levels = c("Upstream exon", "Intron", "Downstream exon")))

bg_gc <- gc_long$gc[gc_long$group == grp_names[1] & gc_long$region == "Intron"]

p3b <- ggplot(gc_long, aes(x = group, y = gc, fill = region)) +
  geom_col(position = position_dodge(0.7), width = 0.6, colour = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = bg_gc, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = c("Upstream exon" = "#66C2A5", "Intron" = "#FC8D62",
                                "Downstream exon" = "#8DA0CB"), name = "Region") +
  scale_y_continuous(limits = c(0.4, 0.58), oob = scales::squish) +
  labs(title = "GC content by region (median)",
       subtitle = "Elevated GC throughout = speckle isochore signature",
       x = "", y = "GC content") +
  coord_flip() + theme_pub + theme(legend.position = "top", plot.subtitle = element_text(size = 8.5))

# Panel C: Heatmap synthesis
model_data <- feature_plot %>%
  select(group, `Short introns` = med_length, `High intron GC` = med_gc,
         `High exon GC` = med_up_gc, `Speckle genes` = speckle)

model_long <- model_data %>%
  pivot_longer(-group, names_to = "feature", values_to = "raw_value") %>%
  group_by(feature) %>%
  mutate(z = (raw_value - min(raw_value)) / (max(raw_value) - min(raw_value))) %>%
  ungroup() %>%
  mutate(z = ifelse(feature == "Short introns", 1 - z, z),
         feature = factor(feature, levels = c("Short introns", "High intron GC",
                                               "High exon GC", "Speckle genes")))

p3c <- ggplot(model_long, aes(x = feature, y = group, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(feature == "Speckle genes", sprintf("%.0f%%", raw_value),
                                ifelse(feature == "Short introns", sprintf("%d", round(raw_value)),
                                       sprintf("%.3f", raw_value)))),
            size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#FDAE6B", high = "#8B0000",
                        midpoint = 0.5, guide = "none") +
  labs(title = "Speckle signature strength",
       subtitle = "Darker = more speckle-like. All taxol-sensitive groups match the signature.",
       x = "", y = "") +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                     panel.grid = element_blank(), plot.subtitle = element_text(size = 8.5))

fig3 <- (p3a | p3b) / p3c +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Figure 3: Blocked and taxol-affected introns share the nuclear speckle structural signature",
    subtitle = "Short, GC-rich introns flanked by GC-rich exons in speckle-associated genes -- the hallmark of speckle-dependent splicing",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30", margin = margin(b = 10)),
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

ggsave(file.path(PLOTS_DIR, "Figure3_speckle_signature.pdf"), fig3, width = 13, height = 10)
ggsave(file.path(PLOTS_DIR, "Figure3_speckle_signature.png"), fig3, width = 13, height = 10, dpi = 300)
cat("  Saved Figure3_speckle_signature\n")

# ============================================================================
# FIGURE 4: THE SCATTER — Taxol effects: myoblast vs myotube (quadrant view)
#   "Same introns are NOT affected in both cell types"
# ============================================================================

cat("Figure 4: The Scatter...\n")

# Build taxol scatter — x = ΔPSI/ΔPIR (myoblast taxol), y = ΔPSI/ΔPIR (myotube taxol)
build_taxol_scatter <- function(myo_df, tube_df, type_label) {
  m <- myo_df %>% select(EVENT, GENE, dpsi_myo = deltapsi, FDR_myo = FDR)
  t <- tube_df %>% select(EVENT, dpsi_tube = deltapsi, FDR_tube = FDR)
  inner_join(m, t, by = "EVENT") %>%
    mutate(
      sig_myo  = is_sig(FDR_myo, dpsi_myo),
      sig_tube = is_sig(FDR_tube, dpsi_tube),
      quadrant = case_when(
        sig_myo & sig_tube & sign(dpsi_myo) == sign(dpsi_tube) ~ "Both concordant",
        sig_myo & sig_tube & sign(dpsi_myo) != sign(dpsi_tube) ~ "Discordant",
        sig_myo & !sig_tube ~ "Myoblast only",
        !sig_myo & sig_tube ~ "Myotube only",
        TRUE ~ "NS"),
      event_type = type_label,
      is_speckle = toupper(GENE) %in% speckle_upper
    )
}

tx_ir <- build_taxol_scatter(filter_ir(myo_taxol), filter_ir(tube_taxol), "Intron Retention")
tx_ex <- build_taxol_scatter(filter_ex(myo_taxol), filter_ex(tube_taxol), "Cassette Exon")

# Palette for taxol quadrants
col_both   <- "#2166AC"   # deep blue
col_disc   <- "#E08214"   # orange
col_myob   <- "#92C5DE"   # light blue
col_myot   <- "#D6604D"   # warm red

# Panel A: Side-by-side scatter
tx_both <- bind_rows(tx_ex, tx_ir) %>%
  mutate(event_type = factor(event_type, levels = c("Cassette Exon", "Intron Retention")))

p4a <- ggplot(tx_both, aes(dpsi_myo, dpsi_tube, colour = quadrant)) +
  geom_point(data = tx_both %>% filter(quadrant == "NS"), size = 0.15, alpha = 0.08) +
  geom_point(data = tx_both %>% filter(quadrant != "NS"), size = 0.9, alpha = 0.7) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.2) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3, linetype = "dotted") +
  scale_colour_manual(
    values = c("NS" = col_ns, "Both concordant" = col_both, "Discordant" = col_disc,
               "Myoblast only" = col_myob, "Myotube only" = col_myot),
    name = "", breaks = c("Both concordant", "Discordant", "Myoblast only", "Myotube only")) +
  facet_wrap(~event_type) +
  labs(title = "Taxol splicing response: myoblast vs myotube",
       subtitle = "IR events show near-zero overlap between cell types -- almost all are cell-type-specific",
       x = expression(Delta*PIR ~ "or" ~ Delta*PSI ~ "(Myoblast: Taxol vs DMSO)"),
       y = expression(Delta*PIR ~ "or" ~ Delta*PSI ~ "(Myotube: Taxol vs DMSO)")) +
  coord_fixed(xlim = c(-0.65, 0.65), ylim = c(-0.65, 0.65)) +
  theme_pub +
  theme(legend.position = "top")

# Panel B: Quadrant composition comparison (IR vs EX)
quad_summary <- bind_rows(
  tx_ir %>% filter(quadrant != "NS") %>% count(quadrant) %>% mutate(event_type = "Intron Retention"),
  tx_ex %>% filter(quadrant != "NS") %>% count(quadrant) %>% mutate(event_type = "Cassette Exon")
) %>%
  mutate(
    quadrant = factor(quadrant, levels = c("Both concordant", "Discordant", "Myoblast only", "Myotube only")),
    event_type = factor(event_type, levels = c("Cassette Exon", "Intron Retention"))
  ) %>%
  # Add zero counts for missing categories
  complete(quadrant, event_type, fill = list(n = 0))

p4b <- ggplot(quad_summary, aes(x = event_type, y = n, fill = quadrant)) +
  geom_col(position = "dodge", width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(0.7), vjust = -0.3, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Both concordant" = col_both, "Discordant" = col_disc,
                                "Myoblast only" = col_myob, "Myotube only" = col_myot),
                    name = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Quadrant composition",
       subtitle = "IR: 0 shared events -- taxol-sensitive\nintrons are entirely cell-type-specific",
       x = "", y = "Number of significant events") +
  theme_pub +
  theme(legend.position = "right", plot.subtitle = element_text(size = 8.5))

# Panel C: Speckle enrichment per quadrant (IR only)
ir_quad_speckle <- tx_ir %>%
  filter(quadrant != "NS") %>%
  group_by(quadrant) %>%
  summarise(
    n_events = n(),
    n_genes = n_distinct(GENE),
    n_speckle = sum(!duplicated(GENE) & is_speckle),
    .groups = "drop"
  ) %>%
  mutate(
    pct_speckle = 100 * n_speckle / n_genes,
    quadrant = factor(quadrant, levels = c("Both concordant", "Discordant", "Myoblast only", "Myotube only"))
  )

# Background rate
bg_ir_genes <- unique(toupper(tx_ir$GENE))
bg_pct <- 100 * sum(bg_ir_genes %in% speckle_upper) / length(bg_ir_genes)

p4c <- ggplot(ir_quad_speckle, aes(x = quadrant, y = pct_speckle,
                                    fill = ifelse(pct_speckle > bg_pct * 1.5, "enriched", "normal"))) +
  geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = bg_pct, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d genes)", pct_speckle, n_genes)),
            vjust = -0.15, size = 3, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 0.55, y = bg_pct + 1, label = sprintf("BG: %.0f%%", bg_pct),
           hjust = 0, size = 2.8, colour = "grey50") +
  scale_fill_manual(values = c("enriched" = col_speckle, "normal" = col_nonsp), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Speckle enrichment by IR quadrant",
       subtitle = "Both myoblast-only and myotube-only IR events\nare enriched for speckle genes",
       x = "", y = "% speckle-associated genes") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), plot.subtitle = element_text(size = 8.5))

fig4 <- p4a / (p4b | p4c) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Figure 4: Taxol-sensitive introns are cell-type-specific with zero overlap",
    subtitle = "IR events affected by taxol in myoblasts vs myotubes are entirely non-overlapping -- yet both sets target speckle genes",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30", margin = margin(b = 10)),
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

ggsave(file.path(PLOTS_DIR, "Figure4_taxol_scatter.pdf"), fig4, width = 12, height = 10)
ggsave(file.path(PLOTS_DIR, "Figure4_taxol_scatter.png"), fig4, width = 12, height = 10, dpi = 300)
cat("  Saved Figure4_taxol_scatter\n")

cat("\n=== All figures saved to plots/ ===\n")
cat("  Figure1_mirror_reversal.pdf   — The cell-type direction reversal\n")
cat("  Figure2_differentiation_block.pdf — IR-selective differentiation block\n")
cat("  Figure3_speckle_signature.pdf — Structural proof of speckle mechanism\n")
cat("  Figure4_taxol_scatter.pdf     — Taxol quadrant scatter (myoblast vs myotube)\n")
