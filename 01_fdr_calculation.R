# ============================================================================
# 01_fdr_calculation.R
# Taxol Effect on mRNA Splicing in C2C12 Myoblasts and Myotubes
# FDR Calculation Using betAS (1000 Simulations)
#
# PURPOSE: Compute FDR-corrected differential splicing statistics for all
#          pairwise comparisons using betAS::prepareTableVolcanoFDR.
#          This script is designed to run on a cluster (no plotting).
#
# EXPERIMENTAL DESIGN:
#   12 samples total:
#   - 3 Myoblast DMSO (control)
#   - 3 Myoblast Taxol
#   - 3 Myotube DMSO (control)
#   - 3 Myotube Taxol
#
# COMPARISONS:
#   1. Myoblast:  Taxol vs DMSO (Taxol effect in undifferentiated cells)
#   2. Myotube:   Taxol vs DMSO (Taxol effect in differentiated cells)
#   3. DMSO:      Myotube vs Myoblast (differentiation effect, no drug)
#   4. Taxol:     Myotube vs Myoblast (differentiation effect under Taxol)
#
# EVENT TYPES: All types processed together (C1, C2, C3, S, MIC, IR, ANN, ALTD, ALTA)
#
# USAGE:
#   Rscript 01_fdr_calculation.R
#
# Author: Andrés Gordo Ortiz
# ============================================================================

cat(strrep("=", 70), "\n")
cat("01_fdr_calculation.R — betAS FDR Calculation (nsim = 1000)\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(betAS)
  library(dplyr)
  library(readr)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

NSIM        <- 1000     # Number of beta distribution simulations
NPOINTS     <- 500      # Resolution of the density estimation
MIN_READS   <- 10       # Minimum reads to keep an event (N parameter)
SEED        <- 42       # Reproducibility
BASECOL     <- "#89C0AE"
INTCOL      <- "#E69A9C"

set.seed(SEED)

# ============================================================================
# PATHS — Adjust these to your cluster paths
# ============================================================================

INCLUSION_TABLE <- file.path(getwd(), "inclusion_tables",
                             "INCLUSION_LEVELS_FULL-mm10.tab")
METADATA_FILE   <- file.path(getwd(), "metadata", "metadata.csv")
RESULTS_DIR     <- file.path(getwd(), "results")

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

cat("Inclusion table: ", INCLUSION_TABLE, "\n")
cat("Metadata:        ", METADATA_FILE, "\n")
cat("Results dir:     ", RESULTS_DIR, "\n\n")

# ============================================================================
# LOAD DATA
# ============================================================================

cat("--- Loading inclusion table ---\n")
data <- getDataset(pathTables = INCLUSION_TABLE, tool = "vast-tools")

cat("--- Extracting events (min reads = ", MIN_READS, ") ---\n")
all_events <- filterEvents(getEvents(data, tool = "vast-tools"), N = MIN_READS)

cat(sprintf("  Total events: %d\n", nrow(all_events$PSI)))

# ============================================================================
# METADATA & GROUP SETUP
# ============================================================================

cat("\n--- Loading metadata ---\n")
metadata <- read.csv(METADATA_FILE, stringsAsFactors = FALSE) %>%
  arrange(condition)

cat("  Samples:    ", nrow(metadata), "\n")
cat("  Conditions: ", paste(unique(metadata$condition), collapse = ", "), "\n\n")

# Build group list (same pattern as myoblast_updated.Rmd)
groups <- unique(metadata$condition)
samples_all <- metadata$sample_id

groupList <- lapply(seq_along(groups), function(i) {
  list(
    name    = groups[i],
    samples = metadata$sample_id[metadata$condition == groups[i]],
    color   = c("#4BA3C3", "#D62839", "#2166AC", "#B2182B")[i]
  )
})
names(groupList) <- groups

# Print sample assignment
for (g in names(groupList)) {
  cat(sprintf("  %s: %s\n", g, paste(groupList[[g]]$samples, collapse = ", ")))
}

# ============================================================================
# DEFINE COMPARISONS
# ============================================================================

comparisons <- list(
  # Taxol effect in Myoblasts
  list(labA = "Myoblast_DMSO",  labB = "Myoblast_Taxol",
       name = "Myoblast_Taxol_vs_DMSO"),
  # Taxol effect in Myotubes
  list(labA = "Myotube_DMSO",   labB = "Myotube_Taxol",
       name = "Myotube_Taxol_vs_DMSO"),
  # Differentiation effect (no drug)
  list(labA = "Myoblast_DMSO",  labB = "Myotube_DMSO",
       name = "Differentiation_DMSO"),
  # Differentiation effect under Taxol
  list(labA = "Myoblast_Taxol", labB = "Myotube_Taxol",
       name = "Differentiation_Taxol")
)

# ============================================================================
# RUN FDR CALCULATIONS
# ============================================================================

run_fdr <- function(event_data, comp, groupList) {
  samplesA <- groupList[[comp$labA]]$samples
  samplesB <- groupList[[comp$labB]]$samples
  colsA    <- convertCols(event_data$PSI, samplesA)
  colsB    <- convertCols(event_data$PSI, samplesB)

  outname <- comp$name
  cat(sprintf("  Running: %s (%s vs %s) ...\n", outname, comp$labA, comp$labB))

  t_start <- Sys.time()

  result <- tryCatch({
    prepareTableVolcanoFDR(
      psitable      = event_data$PSI,
      qualtable     = event_data$Qual,
      npoints       = NPOINTS,
      colsA         = colsA,
      colsB         = colsB,
      labA          = comp$labA,
      labB          = comp$labB,
      basalColor    = BASECOL,
      interestColor = INTCOL,
      maxDevTable   = maxDevSimulationN100,
      nsim          = NSIM,
      seed          = TRUE,
      CoverageWeight = FALSE
    )
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message))
    return(NULL)
  })

  elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 2)
  cat(sprintf("    Done in %s min\n", elapsed))

  if (!is.null(result)) {
    outfile <- file.path(RESULTS_DIR, paste0(outname, ".csv"))
    write.csv(result, outfile, row.names = FALSE)
    cat(sprintf("    Saved: %s (%d events)\n", outfile, nrow(result)))

    # Quick summary
    sig <- result[!is.na(result$FDR) & result$FDR <= 0.05 &
                  abs(result$deltapsi) >= 0.1, ]
    n_inc  <- sum(sig$deltapsi > 0, na.rm = TRUE)
    n_skip <- sum(sig$deltapsi < 0, na.rm = TRUE)
    cat(sprintf("    Significant (FDR<=0.05, |dPSI|>=0.1): %d total (↑%d ↓%d)\n",
                nrow(sig), n_inc, n_skip))
  }

  return(result)
}

cat(strrep("=", 70), "\n")
cat("STARTING FDR CALCULATIONS (nsim = ", NSIM, ")\n")
cat(strrep("=", 70), "\n\n")

all_results <- list()

for (comp in comparisons) {
  cat(strrep("-", 50), "\n")
  cat("Comparison: ", comp$name, "\n")
  cat(strrep("-", 50), "\n")

  all_results[[comp$name]] <-
    run_fdr(all_events, comp, groupList)

  cat("\n")
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat(strrep("=", 70), "\n")
cat("SUMMARY OF ALL COMPARISONS\n")
cat(strrep("=", 70), "\n\n")

summary_rows <- list()
for (name in names(all_results)) {
  res <- all_results[[name]]
  if (!is.null(res)) {
    sig <- res[!is.na(res$FDR) & res$FDR <= 0.05 & abs(res$deltapsi) >= 0.1, ]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      comparison   = name,
      total_events = nrow(res),
      significant  = nrow(sig),
      included     = sum(sig$deltapsi > 0, na.rm = TRUE),
      skipped      = sum(sig$deltapsi < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
print(summary_df)
write.csv(summary_df, file.path(RESULTS_DIR, "fdr_summary.csv"), row.names = FALSE)

cat("\n")
cat(strrep("=", 70), "\n")
cat("ALL FDR CALCULATIONS COMPLETE\n")
cat("Results saved to: ", RESULTS_DIR, "\n")
cat(strrep("=", 70), "\n")
cat("\nSession info:\n")
sessionInfo()
