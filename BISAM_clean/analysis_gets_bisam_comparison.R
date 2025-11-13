################################################################################
# BISAM vs GETS Comparison Analysis
# 
# Description: Comprehensive comparison of BISAM and GETS methods across 
#              different threshold levels (SD breakpoints)
# 
# Output: Performance metrics, plots, and summary statistics
################################################################################

# ==============================================================================
# 1. SETUP AND DATA LOADING
# ==============================================================================

# Clear workspace
rm(list = ls())

# Load required libraries (only for data manipulation)
library(dplyr)

# Set up paths
data_path <- "./Simulations/2025-11-06/gets_bisam_comparison_gets-0.05_bisam_prior-imom/"

# Get list of all RDS files in the folder
file_list <- list.files(
  path = data_path, 
  pattern = "\\.RDS$", 
  full.names = TRUE
)

# ==============================================================================
# 2. DATA PROCESSING
# ==============================================================================

# Read and combine all simulation results
combined_data <- lapply(file_list, function(file) {
  # Read the matrix
  mat <- readRDS(file)
  
  # Extract metadata from filename
  filename <- basename(file)
  breaksize <- as.numeric(sub(".*breaksize-([0-9.]+)SD_rep.*", "\\1", filename))
  replicate <- as.integer(sub(".*_rep([0-9]+)\\.RDS", "\\1", filename))
  
  # Convert to data frame with identifiers
  df <- as.data.frame(mat)
  df$row_name <- rownames(mat)
  df$breaksize <- breaksize
  df$replicate <- replicate
  
  return(df)
}) |> 
  bind_rows()

# ==============================================================================
# 3. SUMMARY STATISTICS
# ==============================================================================

# Get unique breaksizes and sort them
breaksize_vals <- combined_data$breaksize |> unique() |> sort()

# Calculate summary table for each breaksize
summary_table <- sapply(breaksize_vals, function(bs) {
  data_subset <- combined_data |>
    filter(breaksize == bs)
  
  # Sum across all replicates
  summary_stats <- data_subset |>
    dplyr::select(-row_name, -breaksize, -replicate) |> 
    colSums()
  
  return(summary_stats)
})

# Set column names
colnames(summary_table) <- sprintf("%.1fSD", breaksize_vals)

# Display summary table
cat("\n\nSummary Table:\n")
print(summary_table)

# ==============================================================================
# 4. PERFORMANCE METRICS CALCULATION
# ==============================================================================

# Define threshold levels
breaksize_labels <- colnames(summary_table)
sd_numeric <- as.numeric(stringr::str_extract(breaksize_labels, "\\d+\\.\\d(?=SD)"))

# Initialize metrics data frame
metrics <- data.frame(
  SD = rep(breaksize_labels, 2),
  SD_num = rep(sd_numeric, 2),
  Method = rep(c("BISAM", "GETS"), each = 6),
  TP = c(summary_table["tr.ssvs", ], summary_table["tr.gets", ]),
  FP = c(summary_table["fp.ssvs", ], summary_table["fp.gets", ]),
  FN = c(summary_table["fn.ssvs", ], summary_table["fn.gets", ]),
  Total_Found = c(summary_table["ssvs", ], summary_table["gets", ]),
  True_Total = rep(summary_table["true", ], 2)
)

# Calculate derived metrics
metrics <- within(metrics, {
  Precision <- TP / (TP + FP)
  Recall <- TP / (TP + FN)
  F1 <- 2 * (Precision * Recall) / (Precision + Recall)
  FDR <- FP / (TP + FP)  # False Discovery Rate
  Specificity <- 1 - FDR
})

# ==============================================================================
# 1. SOPHISTICATED COLOR PALETTE
# ==============================================================================

# Nature/Science-inspired color palette with excellent print/screen reproduction
colors <- list(
  ssvs = "#0072B2",      # Deep blue
  gets = "#D55E00",      # Vermillion/burnt orange  
  window = "#009E73",    # Bluish green
  standard = "#CC79A7",  # Reddish purple
  gray = "#808080",      # Neutral gray
  lightgray = "#E5E5E5", # Light gray for grids
  darkgray = "#404040"   # Dark gray for text
)

# ==============================================================================
# 2. ENHANCED PLOT SETUP FUNCTION
# ==============================================================================

setup_publication_plot <- function(title = "", no_top_margin = FALSE) {
  par(
    mar = if(no_top_margin) c(4.5, 4.8, 2, 1) else c(4.5, 4.8, 3.5, 1),
    mgp = c(3.2, 0.7, 0),
    las = 1,
    cex.lab = 1.15,
    cex.axis = 1.0,
    cex.main = 1.25,
    col.lab = colors$darkgray,
    col.axis = colors$darkgray,
    col.main = colors$darkgray,
    family = "sans",
    tcl = -0.3,
    bty = "n"  # No box
  )
}

# Enhanced axis function without box
add_clean_axes <- function(side = 1:2, at_x = NULL, labels_x = NULL, 
                           at_y = NULL, labels_y = NULL) {
  if (1 %in% side) {
    if (!is.null(at_x)) {
      axis(1, at = at_x, labels = labels_x, col = colors$gray, 
           col.axis = colors$darkgray, lwd = 1.5)
    } else {
      axis(1, col = colors$gray, col.axis = colors$darkgray, lwd = 1.5)
    }
  }
  if (2 %in% side) {
    if (!is.null(at_y)) {
      axis(2, at = at_y, labels = labels_y, col = colors$gray, 
           col.axis = colors$darkgray, lwd = 1.5)
    } else {
      axis(2, col = colors$gray, col.axis = colors$darkgray, lwd = 1.5)
    }
  }
}

# ==============================================================================
# 3. CREATE PUBLICATION-QUALITY PLOTS
# ==============================================================================

# Open PDF with better dimensions for journal publication
pdf("Simulations/ssvs_gets_publication_plots.pdf", width = 14, height = 10)

# Use 2x3 layout with better spacing
par(mfrow = c(2, 3), oma = c(0, 0, 0, 0))

# ------------------------------------------------------------------------------
# Plot 1: Precision and Recall
# ------------------------------------------------------------------------------
setup_publication_plot()

ssvs_idx <- metrics$Method == "BISAM"
gets_idx <- metrics$Method == "GETS"

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Performance Score",
  main = "Precision and Recall",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

# Subtle grid
abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = 0.8)

# Plot Precision (solid lines)
lines(log10(metrics$SD_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
      col = colors$ssvs, lwd = 2.5, lty = 1)
points(log10(metrics$SD_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = 1.3)

lines(log10(metrics$SD_num[gets_idx]), metrics$Precision[gets_idx], 
      col = colors$gets, lwd = 2.5, lty = 1)
points(log10(metrics$SD_num[gets_idx]), metrics$Precision[gets_idx], 
       col = colors$gets, pch = 16, cex = 1.3)

# Plot Recall (dashed lines)
lines(log10(metrics$SD_num[ssvs_idx]), metrics$Recall[ssvs_idx], 
      col = colors$ssvs, lwd = 2.5, lty = 2)
points(log10(metrics$SD_num[ssvs_idx]), metrics$Recall[ssvs_idx], 
       col = colors$ssvs, pch = 17, cex = 1.3)

lines(log10(metrics$SD_num[gets_idx]), metrics$Recall[gets_idx], 
      col = colors$gets, lwd = 2.5, lty = 2)
points(log10(metrics$SD_num[gets_idx]), metrics$Recall[gets_idx], 
       col = colors$gets, pch = 17, cex = 1.3)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM Precision", "GETS Precision", 
                  "BISAM Recall", "GETS Recall"),
       col = c(colors$ssvs, colors$gets, colors$ssvs, colors$gets),
       lty = c(1, 1, 2, 2),
       pch = c(16, 16, 17, 17),
       lwd = 2.5,
       bty = "n",
       cex = 0.95,
       pt.cex = 1.2)

# Panel label
mtext("A", side = 3, line = 1.5, at = par("usr")[1], cex = 1.4, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 2: F1 Score Comparison
# ------------------------------------------------------------------------------
setup_publication_plot()

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "F1 Score",
  main = "F1 Score Comparison",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = 0.8)

lines(log10(metrics$SD_num[ssvs_idx]), metrics$F1[ssvs_idx], 
      col = colors$ssvs, lwd = 2.5)
points(log10(metrics$SD_num[ssvs_idx]), metrics$F1[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = 1.3)

lines(log10(metrics$SD_num[gets_idx]), metrics$F1[gets_idx], 
      col = colors$gets, lwd = 2.5)
points(log10(metrics$SD_num[gets_idx]), metrics$F1[gets_idx], 
       col = colors$gets, pch = 16, cex = 1.3)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM", "GETS"),
       col = c(colors$ssvs, colors$gets),
       lty = 1,
       lwd = 2.5,
       pch = 16,
       bty = "n",
       cex = 0.95,
       pt.cex = 1.2)

mtext("B", side = 3, line = 1.5, at = par("usr")[1], cex = 1.4, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 3: Detection Rates
# ------------------------------------------------------------------------------
setup_publication_plot()

detection_rate_ssvs <- summary_table["ssvs", ] / summary_table["true", ]
detection_rate_gets <- summary_table["gets", ] / summary_table["true", ]
true_detection_rate_ssvs <- summary_table["tr.ssvs", ] / summary_table["true", ]
true_detection_rate_gets <- summary_table["tr.gets", ] / summary_table["true", ]

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1.5),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Detection Rate",
  main = "Detection Rates",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1.5, 0.25), col = colors$lightgray, lty = 1, lwd = 0.8)
abline(h = 1, col = colors$gray, lty = 2, lwd = 1.2)

# True detection rate (solid)
lines(log10(sd_numeric), true_detection_rate_ssvs, 
      col = colors$ssvs, lwd = 2.5, lty = 1)
points(log10(sd_numeric), true_detection_rate_ssvs, 
       col = colors$ssvs, pch = 16, cex = 1.3)

lines(log10(sd_numeric), true_detection_rate_gets, 
      col = colors$gets, lwd = 2.5, lty = 1)
points(log10(sd_numeric), true_detection_rate_gets, 
       col = colors$gets, pch = 16, cex = 1.3)

# Total detection rate (dashed)
lines(log10(sd_numeric), detection_rate_ssvs, 
      col = colors$ssvs, lwd = 2.5, lty = 3)
points(log10(sd_numeric), detection_rate_ssvs, 
       col = colors$ssvs, pch = 1, cex = 1.3, lwd = 1.5)

lines(log10(sd_numeric), detection_rate_gets, 
      col = colors$gets, lwd = 2.5, lty = 3)
points(log10(sd_numeric), detection_rate_gets, 
       col = colors$gets, pch = 1, cex = 1.3, lwd = 1.5)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1.5, 0.25)
)

legend("topright", 
       legend = c("BISAM True", "GETS True", "BISAM Total", "GETS Total"),
       col = c(colors$ssvs, colors$gets, colors$ssvs, colors$gets),
       lty = c(1, 1, 3, 3),
       pch = c(16, 16, 1, 1),
       lwd = 2.5,
       bty = "n",
       cex = 0.85,
       pt.cex = 1.2)

mtext("C", side = 3, line = 1.5, at = par("usr")[1], cex = 1.4, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 4: True Positives
# ------------------------------------------------------------------------------
setup_publication_plot()

y_vals_tp <- c(summary_table["tr.ssvs", ], summary_table["tr.gets", ])
max_y_tp <- max(y_vals_tp)

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)) + c(-0.15, 0.15),
  ylim = c(0, max_y_tp * 1.15),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Number of Breaks",
  main = "True Positives",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

bar_width <- 0.03
offset <- 0.001

for (j in 1:6) {
  rect(log10(sd_numeric[j]) - bar_width - offset, 0,
       log10(sd_numeric[j]) - offset, summary_table["tr.ssvs", j],
       col = colors$ssvs, border = NA)
  
  rect(log10(sd_numeric[j]) + offset, 0,
       log10(sd_numeric[j]) + bar_width + offset, summary_table["tr.gets", j],
       col = colors$gets, border = NA)
}

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels
)

legend("topright", 
       legend = c("BISAM", "GETS"),
       fill = c(colors$ssvs, colors$gets),
       border = NA,
       bty = "n",
       cex = 0.95)

mtext("D", side = 3, line = 1.5, at = par("usr")[1], cex = 1.4, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 5: False Positives
# ------------------------------------------------------------------------------
setup_publication_plot()

y_vals_fp <- c(summary_table["fp.ssvs", ], summary_table["fp.gets", ])
max_y_fp <- max(y_vals_fp)

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)) + c(-0.15, 0.15),
  ylim = c(0, max_y_fp * 1.15),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Number of Breaks",
  main = "False Positives",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

for (j in 1:6) {
  rect(log10(sd_numeric[j]) - bar_width - offset, 0,
       log10(sd_numeric[j]) - offset, summary_table["fp.ssvs", j],
       col = colors$ssvs, border = NA)
  
  rect(log10(sd_numeric[j]) + offset, 0,
       log10(sd_numeric[j]) + bar_width + offset, summary_table["fp.gets", j],
       col = colors$gets, border = NA)
}

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels
)

legend("topright", 
       legend = c("BISAM", "GETS"),
       fill = c(colors$ssvs, colors$gets),
       border = NA,
       bty = "n",
       cex = 0.95)

mtext("E", side = 3, line = 1.5, at = par("usr")[1], cex = 1.4, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 6: BISAM Window Analysis
# ------------------------------------------------------------------------------
setup_publication_plot()

standard_rate <- summary_table["tr.ssvs", ] / summary_table["true", ]
window_rate <- summary_table["ssvs_t3_tr", ] / summary_table["true", ]

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)) + c(-0.15, 0.15),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "True Positive Rate",
  main = "BISAM: Standard vs Window (±2)",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = 0.8)

for (j in 1:6) {
  rect(log10(sd_numeric[j]) - bar_width - offset, 0,
       log10(sd_numeric[j]) - offset, standard_rate[j],
       col = colors$standard, border = NA)
  
  rect(log10(sd_numeric[j]) + offset, 0,
       log10(sd_numeric[j]) + bar_width + offset, window_rate[j],
       col = colors$window, border = NA)
}

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("Standard (exact)", "Window (±2)"),
       fill = c(colors$standard, colors$window),
       border = NA,
       bty = "n",
       cex = 0.95)

mtext("F", side = 3, line = 1.5, at = par("usr")[1], cex = 1.4, font = 2, adj = 0)

dev.off()

# ==============================================================================
# 4. CREATE INDIVIDUAL HIGH-RESOLUTION PLOTS
# ==============================================================================

# Individual plot for Precision and Recall (larger format)
pdf("Simulations/precision_recall_publication.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
setup_publication_plot()

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Performance Score",
  main = "Precision and Recall",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = 0.8)

lines(log10(metrics$SD_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
      col = colors$ssvs, lwd = 3, lty = 1)
points(log10(metrics$SD_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = 1.5)

lines(log10(metrics$SD_num[gets_idx]), metrics$Precision[gets_idx], 
      col = colors$gets, lwd = 3, lty = 1)
points(log10(metrics$SD_num[gets_idx]), metrics$Precision[gets_idx], 
       col = colors$gets, pch = 16, cex = 1.5)

lines(log10(metrics$SD_num[ssvs_idx]), metrics$Recall[ssvs_idx], 
      col = colors$ssvs, lwd = 3, lty = 2)
points(log10(metrics$SD_num[ssvs_idx]), metrics$Recall[ssvs_idx], 
       col = colors$ssvs, pch = 17, cex = 1.5)

lines(log10(metrics$SD_num[gets_idx]), metrics$Recall[gets_idx], 
      col = colors$gets, lwd = 3, lty = 2)
points(log10(metrics$SD_num[gets_idx]), metrics$Recall[gets_idx], 
       col = colors$gets, pch = 17, cex = 1.5)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM Precision", "GETS Precision", 
                  "BISAM Recall", "GETS Recall"),
       col = c(colors$ssvs, colors$gets, colors$ssvs, colors$gets),
       lty = c(1, 1, 2, 2),
       pch = c(16, 16, 17, 17),
       lwd = 3,
       bty = "n",
       cex = 1.1,
       pt.cex = 1.4)

dev.off()

# Individual plot for F1 Score
pdf("Simulations/f1_score_publication.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
setup_publication_plot()

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "F1 Score",
  main = "F1 Score Comparison",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = 0.8)

lines(log10(metrics$SD_num[ssvs_idx]), metrics$F1[ssvs_idx], 
      col = colors$ssvs, lwd = 3)
points(log10(metrics$SD_num[ssvs_idx]), metrics$F1[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = 1.5)

lines(log10(metrics$SD_num[gets_idx]), metrics$F1[gets_idx], 
      col = colors$gets, lwd = 3)
points(log10(metrics$SD_num[gets_idx]), metrics$F1[gets_idx], 
       col = colors$gets, pch = 16, cex = 1.5)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM", "GETS"),
       col = c(colors$ssvs, colors$gets),
       lty = 1,
       lwd = 3,
       pch = 16,
       bty = "n",
       cex = 1.1,
       pt.cex = 1.4)

dev.off()

# ==============================================================================
# 5. FINAL OUTPUT MESSAGE
# ==============================================================================

cat("\n\n=====================================\n")
cat("PUBLICATION-QUALITY PLOTS CREATED!\n")
cat("=====================================\n")
cat("Files saved:\n")
cat("- Multi-panel plot: ssvs_gets_publication_plots.pdf (14×10 in)\n")
cat("- Individual plots:\n")
cat("  * precision_recall_publication.pdf (8×6 in)\n")
cat("  * f1_score_publication.pdf (8×6 in)\n")
cat("=====================================\n")
cat("Design features:\n")
cat("- No boxes (clean open design)\n")
cat("- Sophisticated color palette\n")
cat("- Better aspect ratios\n")
cat("- Subtle gridlines\n")
cat("- Panel labels (A-F)\n")
cat("- Enhanced typography\n")
cat("=====================================\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================

