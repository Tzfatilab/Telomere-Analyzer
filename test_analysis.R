#!/usr/bin/env Rscript
# =============================================================================
# test_analysis.R
# Standalone test for the --analysis post-processing step of NanoTel.R
#
# Usage:
#   Rscript test_analysis.R --summary_csv "path/to/BCxx_summary.csv" \
#                           --save_path   "path/to/output_dir"        \
#                           [--bias_prediction_model "path/to/poly_regression_model.rds"]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggprism)
  library(survival)
})

BUFFER <- 50L

.script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("--file=", args, value = TRUE)
  if (length(file_flag) > 0) dirname(normalizePath(sub("--file=", "", file_flag[1]))) else getwd()
}, error = function(e) getwd())

option_list <- list(
  make_option("--summary_csv", type = "character", default = NULL,
              help = "Path to the existing BCxx_summary.csv (or filtered_sorted) file."),
  make_option("--save_path", type = "character", default = NULL,
              help = "Output directory for the generated files."),
  make_option("--bias_prediction_model", type = "character", default = NULL,
              help = "Path to poly_regression_model.rds. Defaults to models/ next to this script.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$summary_csv) || is.null(opt$save_path)) {
  stop("Both --summary_csv and --save_path are required.")
}
if (!file.exists(opt$summary_csv)) {
  stop(paste("summary CSV not found:", opt$summary_csv))
}

dir.create(opt$save_path, recursive = TRUE, showWarnings = FALSE)

# Extract barcode name from filename -> "barcodeNN"
raw_name     <- sub("_(summary|filtered_sorted|filtered_sorted_summary)\\.csv$", "",
                    basename(opt$summary_csv))
barcode_num  <- gsub("[^0-9]", "", raw_name)
barcode_name <- paste0("barcode", barcode_num)
cat("Barcode name detected:", barcode_name, "\n")

# ---- Load CSV ---------------------------------------------------------------
df_summary <- read_csv(opt$summary_csv, show_col_types = FALSE)
cat("Rows in input CSV:", nrow(df_summary), "\n")

# ---- Step 1: Filter ---------------------------------------------------------
df_step1_filtered <- df_summary %>%
  dplyr::filter(telo_density_mismatch >= 0.75,
                Telomere_start_mismatch <= 134)
cat("Rows after density/start filter:", nrow(df_step1_filtered), "\n")

# KM median is computed after the density/start filter and before the
# running-median filtration steps.
df_km <- df_step1_filtered %>%
  dplyr::mutate(
    margin = sequence_length - Telomere_end_mismatch,
    event  = as.integer(margin >= BUFFER)
  )
km_fit    <- survfit(Surv(Telomere_length_mismatch, event) ~ 1, data = df_km)
km_median <- summary(km_fit)$table[["median"]]

# ---- Step 2: Sort by sequence_length descending -----------------------------
df_filtered <- df_step1_filtered %>%
  dplyr::arrange(desc(sequence_length))

# ---- Step 3: Running median + difference column -----------------------------
df_filtered <- df_filtered %>%
  dplyr::mutate(
    TelLenMM_RunningMed = sapply(seq_len(dplyr::n()),
                                  function(i) median(Telomere_length_mismatch[1:i])),
    SeqLen_minus_RunMed = sequence_length - TelLenMM_RunningMed
  )

df_for_plot <- df_filtered %>%
  dplyr::mutate(read_index = seq_len(dplyr::n()))

# ---- Step 4: Remove rows where SeqLen_minus_RunMed < 134 -------------------
df_filtered <- df_filtered %>%
  dplyr::filter(SeqLen_minus_RunMed >= 134)
cat("Rows after SeqLen_minus_RunMed filter:", nrow(df_filtered), "\n")

# ---- Write filtered_sorted_summary.csv --------------------------------------
out_csv <- file.path(opt$save_path, paste0(barcode_name, "_filtered_sorted_summary.csv"))
write_csv(df_filtered, out_csv)
cat("Written:", out_csv, "\n")

# ---- Stats .txt -------------------------------------------------------------

# Event: 1 = complete (telomere fully captured), 0 = censored
df_filtered <- df_filtered %>%
  dplyr::mutate(
    margin = sequence_length - Telomere_end_mismatch,
    event  = as.integer(margin >= BUFFER)
  )

n_reads        <- nrow(df_filtered)
n_complete     <- sum(df_filtered$event == 1)
n_censored     <- sum(df_filtered$event == 0)
censoring_rate <- n_censored / n_reads

med_telo  <- median(df_filtered$Telomere_length_mismatch)
pct_short <- round(100 * sum(df_filtered$Telomere_length_mismatch < 2000) / n_reads, 1)
min_len   <- min(df_filtered$sequence_length)
max_len   <- max(df_filtered$sequence_length)

# WORK IN PROGRESS: the model-based expected KM bias calculation below is a
# provisional working feature, not the final calibrated implementation.
# Keep it functional for now, but treat this block as subject to revision.
# Expected KM bias from polynomial calibration model
model_path <- if (!is.null(opt$bias_prediction_model)) {
  opt$bias_prediction_model
} else {
  file.path(.script_dir, "models", "poly_regression_model.rds")
}
calibration_obj <- readRDS(model_path)
expected_bias   <- predict(
  calibration_obj$calibration_models$fit_km_bias,
  newdata = data.frame(
    censoring_for_model = censoring_rate,
    log10_n_reads       = log10(n_reads)
  )
)
bias_label <- ifelse(expected_bias >= 0,
                     paste0("+", round(expected_bias), " bp"),
                     paste0(round(expected_bias), " bp"))

fmt <- function(x) format(round(x), big.mark = ",", scientific = FALSE)

results_lines <- c(
  paste0("Results for ", barcode_name),
  "==========================================",
  paste0("Total Reads             : ", fmt(n_reads)),
  paste0("Complete Reads          : ", fmt(n_complete)),
  paste0("Censored Reads          : ", fmt(n_censored)),
  paste0("Censoring Rate          : ", round(100 * censoring_rate, 1), "%"),
  "",
  paste0("Regular Median (post-filtration) : ", fmt(med_telo), " bp"),
  paste0("KM Median               : ", fmt(km_median), " bp"),
  "",
  paste0("Expected KM Median Bias : ", bias_label),
  paste0("Min Read Length         : ", fmt(min_len), " bp"),
  paste0("Max Read Length         : ", fmt(max_len), " bp"),
  paste0("% < 2kb                 : ", pct_short, "%")
)

out_txt <- file.path(opt$save_path, paste0(barcode_name, "_results.txt"))
write_lines(results_lines, out_txt)
cat("Written:", out_txt, "\n")
cat("\n--- Stats preview ---\n")
cat(paste(results_lines, collapse = "\n"), "\n")

# ---- Telomere plot ----------------------------------------------------------
p_telo <- ggplot(df_for_plot, aes(x = read_index)) +
  geom_line(aes(y = sequence_length,          color = "Read Length")) +
  geom_line(aes(y = Telomere_length_mismatch, color = "Telomere Length (mismatch)")) +
  geom_line(aes(y = TelLenMM_RunningMed,      color = "Running Median Telomere Length")) +
  scale_color_manual(
    values = c("Read Length"                    = "#E8735A",
               "Telomere Length (mismatch)"     = "#228B22",
               "Running Median Telomere Length" = "#4169E1")
  ) +
  labs(title = "Telomere Analysis",
       x     = "Read (sorted by length, longest to shortest)",
       y     = "Length (bp)",
       color = NULL) +
  theme_prism() +
  theme(legend.position = "bottom")

out_png <- file.path(opt$save_path, paste0(barcode_name, "_telomere_plot.png"))
ggsave(filename = out_png, plot = p_telo, width = 12, height = 6, dpi = 150)
cat("Written:", out_png, "\n")

cat("\nDone. All files written to:", opt$save_path, "\n")
