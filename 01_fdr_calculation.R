# ============================================================================
# 01_fdr_calculation.R
# Taxol Effect on mRNA Splicing in C2C12 Myoblasts and Myotubes
# FDR Calculation Using betAS (1000 Simulations)
#
# PURPOSE: Compute FDR-corrected differential splicing statistics for a
#          single pairwise comparison using betAS::prepareTableVolcanoFDR.
#          Designed to run as a SLURM job array (one task per comparison).
#
# EXPERIMENTAL DESIGN:
#   12 samples total:
#   - 3 Myoblast DMSO (control)
#   - 3 Myoblast Taxol
#   - 3 Myotube DMSO (control)
#   - 3 Myotube Taxol
#
# COMPARISONS (mapped to SLURM_ARRAY_TASK_ID 1-4):
#   1. Myoblast:  Taxol vs DMSO (Taxol effect in undifferentiated cells)
#   2. Myotube:   Taxol vs DMSO (Taxol effect in differentiated cells)
#   3. DMSO:      Myotube vs Myoblast (differentiation effect, no drug)
#   4. Taxol:     Myotube vs Myoblast (differentiation effect under Taxol)
#
# EVENT TYPES: All types processed together (C1, C2, C3, S, MIC, IR, ANN, ALTD, ALTA)
#
# USAGE:
#   Rscript 01_fdr_calculation.R <comparison_index>  # 1-4
#   Or set SLURM_ARRAY_TASK_ID environment variable
#
# Author: Andrés Gordo Ortiz
# ============================================================================

# ============================================================================
# PARSE COMPARISON INDEX (from CLI argument or SLURM_ARRAY_TASK_ID)
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  COMP_INDEX <- as.integer(args[1])
} else if (Sys.getenv("SLURM_ARRAY_TASK_ID") != "") {
  COMP_INDEX <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
} else {
  stop("No comparison index provided. Pass as argument or set SLURM_ARRAY_TASK_ID.")
}

if (is.na(COMP_INDEX) || COMP_INDEX < 1 || COMP_INDEX > 4) {
  stop("Comparison index must be between 1 and 4. Got: ", COMP_INDEX)
}

# ============================================================================
# RESOLVE SCRIPT DIRECTORY (robust for SLURM / Singularity)
# getwd() may not be the project dir on the cluster. We derive it from
# the script path itself, or fall back to SLURM_SUBMIT_DIR, then getwd().
# ============================================================================

SCRIPT_DIR <- tryCatch({
  # Works when called as: Rscript /path/to/01_fdr_calculation.R
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    normalizePath(dirname(sub("^--file=", "", file_arg[1])))
  } else {
    stop("no --file arg")
  }
}, error = function(e) {
  # Fallback: SLURM_SUBMIT_DIR -> getwd()
  d <- Sys.getenv("SLURM_SUBMIT_DIR", unset = getwd())
  normalizePath(d)
})

cat(sprintf("Resolved project directory: %s\n", SCRIPT_DIR))
setwd(SCRIPT_DIR)

cat(strrep("=", 70), "\n")
cat(sprintf("01_fdr_calculation.R -- betAS FDR (nsim=1000) -- Comparison %d/4\n",
            COMP_INDEX))
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

INCLUSION_TABLE <- file.path(SCRIPT_DIR, "INCLUSION_LEVELS_FULL-mm10-12.tab")
METADATA_FILE   <- file.path(SCRIPT_DIR, "metadata", "metadata.csv")
RESULTS_DIR     <- file.path(SCRIPT_DIR, "results")

# showWarnings=FALSE: safe when multiple array tasks create this simultaneously
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE,
                                         showWarnings = FALSE)

cat("Inclusion table: ", INCLUSION_TABLE, "\n")
cat("Metadata:        ", METADATA_FILE, "\n")
cat("Results dir:     ", RESULTS_DIR, "\n\n")

# --- Validate input files exist before doing anything expensive ---
if (!file.exists(INCLUSION_TABLE)) {
  stop("Inclusion table not found: ", INCLUSION_TABLE)
}
if (!file.exists(METADATA_FILE)) {
  stop("Metadata file not found: ", METADATA_FILE)
}

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
    # --- Save RDS first (fast, binary, no precision loss) as safety backup ---
    rds_file <- file.path(RESULTS_DIR, paste0(outname, ".rds"))
    tryCatch({
      saveRDS(result, rds_file)
      cat(sprintf("    RDS backup saved: %s\n", rds_file))
    }, error = function(e) {
      cat(sprintf("    WARNING: RDS save failed: %s\n", e$message))
    })

    # --- Write CSV atomically: write to temp file, then rename ---
    outfile  <- file.path(RESULTS_DIR, paste0(outname, ".csv"))
    tmpfile  <- file.path(RESULTS_DIR, paste0(".", outname, ".csv.tmp"))
    tryCatch({
      write.csv(result, tmpfile, row.names = FALSE)
      file.rename(tmpfile, outfile)
      cat(sprintf("    Saved: %s (%d events)\n", outfile, nrow(result)))
    }, error = function(e) {
      cat(sprintf("    ERROR writing CSV: %s\n", e$message))
      cat(sprintf("    Result is preserved in RDS: %s\n", rds_file))
    })

    # Quick summary
    sig <- result[!is.na(result$FDR) & result$FDR <= 0.05 &
                  abs(result$deltapsi) >= 0.1, ]
    n_inc  <- sum(sig$deltapsi > 0, na.rm = TRUE)
    n_skip <- sum(sig$deltapsi < 0, na.rm = TRUE)
    cat(sprintf("    Significant (FDR<=0.05, |dPSI|>=0.1): %d total (+%d -%d)\n",
                nrow(sig), n_inc, n_skip))
  }

  return(result)
}

cat(strrep("=", 70), "\n")
cat(sprintf("RUNNING FDR CALCULATION (nsim = %d) — Comparison %d\n", NSIM, COMP_INDEX))
cat(strrep("=", 70), "\n\n")

comp <- comparisons[[COMP_INDEX]]

cat(strrep("-", 50), "\n")
cat("Comparison: ", comp$name, "\n")
cat(strrep("-", 50), "\n")

result <- run_fdr(all_events, comp, groupList)

# ============================================================================
# DONE — exit with appropriate status code for SLURM
# ============================================================================

cat("\n")
if (is.null(result)) {
  cat(strrep("=", 70), "\n")
  cat(sprintf("FDR CALCULATION FAILED -- %s\n", comp$name))
  cat(strrep("=", 70), "\n")
  cat("\nSession info:\n")
  print(sessionInfo())
  quit(status = 1, save = "no")
} else {
  cat(strrep("=", 70), "\n")
  cat(sprintf("FDR CALCULATION COMPLETE -- %s\n", comp$name))
  cat("Results saved to: ", RESULTS_DIR, "\n")
  cat(strrep("=", 70), "\n")
  cat("\nSession info:\n")
  print(sessionInfo())
  quit(status = 0, save = "no")
}
