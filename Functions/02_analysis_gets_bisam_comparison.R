################################################################################
# BISAM vs GETS vs ALASSO Comparison Analysis
# 
# Description: Comprehensive comparison of BISAM, GETS, and Adaptive Lasso 
#              methods across different threshold levels (SD breakpoints)
# 
# Output: Performance metrics, plots, and summary statistics
#         Now with SLIDE-FRIENDLY options!
################################################################################

# Clear workspace
rm(list = ls())

# ==============================================================================
# PLOTTING MODE CONFIGURATION
# ==============================================================================

# SET THIS TO EITHER "publication" OR "slide"
PLOT_MODE <- "slide"  # Change to "publication" for journal-style plots

# ==============================================================================
# 1. SETUP AND DATA LOADING
# ==============================================================================

# Load required libraries (only for data manipulation)
library(dplyr)

date <- "2025-12-12_sparse"
gets_lvl <- "0.01"
bisam_prior <- "imom"
tau <- "1.92072941034706"

# Set up paths
if(tau == "") {
  data_path <- sprintf("./Simulations/%s/gets_bisam_alasso_comparison_gets-%s_bisam_prior-%s/", 
                       date, 
                       gets_lvl,
                       bisam_prior)
} else {
  data_path <- sprintf("./Simulations/%s/gets_bisam_alasso_comparison_gets-%s_bisam_prior-%s_tau-%s/", 
                       date, 
                       gets_lvl,
                       bisam_prior,
                       tau)
}


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
  # Read the results
  result_obj <- readRDS(file)
  
  # Check if it's the new format (list) or old format (matrix)
  if (is.list(result_obj) && "break_comparison" %in% names(result_obj)) {
    mat <- result_obj$break_comparison
  } else {
    mat <- result_obj
  }
  
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

# Initialize metrics data frame for all three methods
metrics <- data.frame(
  SD = rep(breaksize_labels, 3),
  SD_num = rep(sd_numeric, 3),
  Method = rep(c("BISAM", "GETS", "ALASSO"), each = length(breaksize_vals)),
  TP = c(summary_table["tr.ssvs", ], 
         summary_table["tr.gets", ],
         summary_table["tr.alasso", ]),
  FP = c(summary_table["fp.ssvs", ], 
         summary_table["fp.gets", ],
         summary_table["fp.alasso", ]),
  FN = c(summary_table["fn.ssvs", ], 
         summary_table["fn.gets", ],
         summary_table["fn.alasso", ]),
  Total_Found = c(summary_table["ssvs", ], 
                  summary_table["gets", ],
                  summary_table["alasso", ]),
  True_Total = rep(summary_table["true", ], 3)
)

# Calculate derived metrics
metrics <- within(metrics, {
  Precision <- TP / (TP + FP)
  Recall <- TP / (TP + FN)
  F1 <- 2 * (Precision * Recall) / (Precision + Recall)
  FDR <- FP / (TP + FP)  # False Discovery Rate
  Specificity <- 1 - FDR
})

# Handle NaN values (when TP + FP = 0 or TP + FN = 0)
metrics$Precision[is.nan(metrics$Precision)] <- 0
metrics$Recall[is.nan(metrics$Recall)] <- 0
metrics$F1[is.nan(metrics$F1)] <- 0

# Display metrics summary
cat("\n\nPerformance Metrics by Method:\n")
print(metrics)

# ==============================================================================
# 5. COLOR PALETTE
# ==============================================================================

# Nature/Science-inspired color palette with excellent print/screen reproduction
colors <- list(
  ssvs = "#0072B2",      # Deep blue
  gets = "#D55E00",      # Vermillion/burnt orange  
  alasso = "#009E73",    # Bluish green (was window color)
  window = "#CC79A7",    # Reddish purple (was standard color)
  standard = "#E69F00",  # Orange (new)
  gray = "#808080",      # Neutral gray
  lightgray = "#E5E5E5", # Light gray for grids
  darkgray = "#404040"   # Dark gray for text
)

# ==============================================================================
# 6. MODE-SPECIFIC SETTINGS
# ==============================================================================

if (PLOT_MODE == "slide") {
  settings <- list(
    # Font sizes
    cex.lab = 1.8,
    cex.axis = 1.6,
    cex.main = 2.0,
    cex.legend = 1.4,
    
    # Line and point sizes
    lwd.line = 4.5,
    lwd.axis = 2.5,
    cex.point = 2.2,
    
    # Grid
    lwd.grid = 1.5,
    
    # Margins
    mar = c(5.5, 6.0, 4.0, 2.0),
    mgp = c(4.0, 1.0, 0),
    
    # PDF dimensions
    pdf.width = 16,
    pdf.height = 9,
    pdf.width.single = 12,
    pdf.height.single = 9,
    
    # File suffix
    suffix = "slides"
  )
} else {  # publication mode
  settings <- list(
    # Font sizes
    cex.lab = 1.15,
    cex.axis = 1.0,
    cex.main = 1.25,
    cex.legend = 0.95,
    
    # Line and point sizes
    lwd.line = 2.5,
    lwd.axis = 1.5,
    cex.point = 1.3,
    
    # Grid
    lwd.grid = 0.8,
    
    # Margins
    mar = c(4.5, 4.8, 3.5, 1),
    mgp = c(3.2, 0.7, 0),
    
    # PDF dimensions
    pdf.width = 14,
    pdf.height = 10,
    pdf.width.single = 8,
    pdf.height.single = 6,
    
    # File suffix
    suffix = "publication"
  )
}

# ==============================================================================
# 7. ENHANCED PLOT SETUP FUNCTION
# ==============================================================================

setup_plot <- function(title = "", no_top_margin = FALSE) {
  par(
    mar = if(no_top_margin) c(settings$mar[1], settings$mar[2], 2, settings$mar[4]) else settings$mar,
    mgp = settings$mgp,
    las = 1,
    cex.lab = settings$cex.lab,
    cex.axis = settings$cex.axis,
    cex.main = settings$cex.main,
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
           col.axis = colors$darkgray, lwd = settings$lwd.axis)
    } else {
      axis(1, col = colors$gray, col.axis = colors$darkgray, lwd = settings$lwd.axis)
    }
  }
  if (2 %in% side) {
    if (!is.null(at_y)) {
      axis(2, at = at_y, labels = labels_y, col = colors$gray, 
           col.axis = colors$darkgray, lwd = settings$lwd.axis)
    } else {
      axis(2, col = colors$gray, col.axis = colors$darkgray, lwd = settings$lwd.axis)
    }
  }
}

# ==============================================================================
# 8. CREATE PLOTS
# ==============================================================================

# Output file names
if(tau == "") {
  multi_panel_file <- sprintf("./Simulations/%s/%s_multi_gets-%s_bisam-%s_tau-%s_%s.pdf", 
                              date, 
                              settings$suffix,
                              gets_lvl,
                              bisam_prior, 
                              "auto")
  f1_score_file <- sprintf("./Simulations/%s/%s_f1_gets-%s_bisam-%s_tau-%s_%s.pdf", 
                           date, 
                           settings$suffix,
                           gets_lvl,
                           bisam_prior, 
                           "auto")
  
} else {
  multi_panel_file <- sprintf("./Simulations/%s/%s_multi_gets-%s_bisam-%s_tau-%s.pdf", 
                              date, 
                              settings$suffix,
                              gets_lvl,
                              bisam_prior, 
                              tau)
  f1_score_file <- sprintf("./Simulations/%s/%s_f1_gets-%s_bisam-%s_tau-%s.pdf", 
                           date, 
                           settings$suffix,
                           gets_lvl,
                           bisam_prior, 
                           tau)
}

# Open PDF for multi-panel plot
pdf(multi_panel_file, width = settings$pdf.width, height = settings$pdf.height)

# Use 2x3 layout
par(mfrow = c(2, 3), oma = c(0, 0, 0, 0))

# Extract indices for each method
ssvs_idx <- metrics$Method == "BISAM"
gets_idx <- metrics$Method == "GETS"
alasso_idx <- metrics$Method == "ALASSO"

# ------------------------------------------------------------------------------
# Plot 1: Precision
# ------------------------------------------------------------------------------
setup_plot()

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Precision",
  main = "Precision",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

# Subtle grid
abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

# Plot Precision
lines(log10(metrics$SD_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line)
points(log10(metrics$SD_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point)

lines(log10(metrics$SD_num[gets_idx]), metrics$Precision[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line)
points(log10(metrics$SD_num[gets_idx]), metrics$Precision[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point)

lines(log10(metrics$SD_num[alasso_idx]), metrics$Precision[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line)
points(log10(metrics$SD_num[alasso_idx]), metrics$Precision[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       pch = 16,
       lwd = settings$lwd.line,
       bty = "n",
       cex = settings$cex.legend,
       pt.cex = settings$cex.point * 0.9)

# Panel label
mtext("A", side = 3, line = 1.5, at = par("usr")[1], cex = settings$cex.main, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 2: Recall
# ------------------------------------------------------------------------------
setup_plot()

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Recall (Sensitivity)",
  main = "Recall",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

# Subtle grid
abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

# Plot Recall
lines(log10(metrics$SD_num[ssvs_idx]), metrics$Recall[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line)
points(log10(metrics$SD_num[ssvs_idx]), metrics$Recall[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point)

lines(log10(metrics$SD_num[gets_idx]), metrics$Recall[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line)
points(log10(metrics$SD_num[gets_idx]), metrics$Recall[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point)

lines(log10(metrics$SD_num[alasso_idx]), metrics$Recall[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line)
points(log10(metrics$SD_num[alasso_idx]), metrics$Recall[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       pch = 16,
       lwd = settings$lwd.line,
       bty = "n",
       cex = settings$cex.legend,
       pt.cex = settings$cex.point * 0.9)

# Panel label
mtext("B", side = 3, line = 1.5, at = par("usr")[1], cex = settings$cex.main, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 3: F1 Score Comparison
# ------------------------------------------------------------------------------
setup_plot()

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

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

lines(log10(metrics$SD_num[ssvs_idx]), metrics$F1[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line)
points(log10(metrics$SD_num[ssvs_idx]), metrics$F1[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point)

lines(log10(metrics$SD_num[gets_idx]), metrics$F1[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line)
points(log10(metrics$SD_num[gets_idx]), metrics$F1[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point)

lines(log10(metrics$SD_num[alasso_idx]), metrics$F1[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line)
points(log10(metrics$SD_num[alasso_idx]), metrics$F1[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       lwd = settings$lwd.line,
       pch = 16,
       bty = "n",
       cex = settings$cex.legend,
       pt.cex = settings$cex.point * 0.9)

mtext("C", side = 3, line = 1.5, at = par("usr")[1], cex = settings$cex.main, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 4: True Positives
# ------------------------------------------------------------------------------
setup_plot()

y_vals_tp <- c(summary_table["tr.ssvs", ], summary_table["tr.gets", ], 
               summary_table["tr.alasso", ])
max_y_tp <- max(y_vals_tp)

bar_width <- if (PLOT_MODE == "slide") 0.03 else 0.02
offset <- if (PLOT_MODE == "slide") 0.005 else 0.003

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

for (j in 1:length(sd_numeric)) {
  # BISAM (left)
  rect(log10(sd_numeric[j]) - 1.5*bar_width - offset, 0,
       log10(sd_numeric[j]) - 0.5*bar_width - offset, summary_table["tr.ssvs", j],
       col = colors$ssvs, border = NA)
  
  # GETS (middle)
  rect(log10(sd_numeric[j]) - bar_width/2, 0,
       log10(sd_numeric[j]) + bar_width/2, summary_table["tr.gets", j],
       col = colors$gets, border = NA)
  
  # ALASSO (right)
  rect(log10(sd_numeric[j]) + 0.5*bar_width + offset, 0,
       log10(sd_numeric[j]) + 1.5*bar_width + offset, summary_table["tr.alasso", j],
       col = colors$alasso, border = NA)
}

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels
)

legend("topleft", 
       legend = c("BISAM", "GETS", "ALASSO"),
       fill = c(colors$ssvs, colors$gets, colors$alasso),
       border = NA,
       bty = "n",
       cex = settings$cex.legend)

mtext("D", side = 3, line = 1.5, at = par("usr")[1], cex = settings$cex.main, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 5: False Positives
# ------------------------------------------------------------------------------
setup_plot()

y_vals_fp <- c(summary_table["fp.ssvs", ], summary_table["fp.gets", ])#, summary_table["fp.alasso", ])
max_y_fp <- max(y_vals_fp)

y1_min <- 0
y1_max <- max_y_fp * 1.15          # left axis (BISAM/GETS)
y2_min <- 0
y2_max <- max(summary_table["fp.alasso", ])

# map ALASSO values to left scale
scale_to_left <- function(z) y1_min + (z - y2_min) * (y1_max - y1_min) / (y2_max - y2_min)
alasso_scaled <- scale_to_left(summary_table["fp.alasso", ])

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)) + c(-0.15, 0.15),
  ylim = c(y1_min, y1_max),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Number of Breaks (BISAM/GETS)",
  main = "False Positives",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

for (j in 1:length(sd_numeric)) {
  # BISAM (left)
  rect(log10(sd_numeric[j]) - 1.5*bar_width - offset, 0,
       log10(sd_numeric[j]) - 0.5*bar_width - offset, summary_table["fp.ssvs", j],
       col = colors$ssvs, border = NA)
  
  # GETS (middle)
  rect(log10(sd_numeric[j]) - bar_width/2, 0,
       log10(sd_numeric[j]) + bar_width/2, summary_table["fp.gets", j],
       col = colors$gets, border = NA)
  
  # ALASSO (right) – scaled to left axis
  rect(log10(sd_numeric[j]) + 0.5*bar_width + offset, 0,
       log10(sd_numeric[j]) + 1.5*bar_width + offset, alasso_scaled[j],
       col = colors$alasso, border = NA)
}

# left x and y axes as before
add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels
)

# add right y-axis for ALASSO
axis(
  side = 4,
  at = pretty(c(y1_min, y1_max)),
  labels = round(y2_min + (pretty(c(y1_min, y1_max)) - y1_min) * (y2_max - y2_min)/(y1_max - y1_min)),
  line = -2   # try 0, then -0.5, -1, ...
)
# mtext("Number of Breaks (ALASSO)", side = 4, line = 1.5)

legend("topright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       fill = c(colors$ssvs, colors$gets, colors$alasso),
       border = NA,
       bty = "n",
       cex = settings$cex.legend)

mtext("E", side = 3, line = 1.5, at = par("usr")[1],
      cex = settings$cex.main, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 6: Near Misses (1-Period Neighbor False Positives)
# ------------------------------------------------------------------------------
setup_plot()

# Calculate near-miss rates (FP within 1 period of true breaks)
near_miss_ssvs <- summary_table["ssvs_1nn_fp", ]
near_miss_gets <- summary_table["gets_1nn_fp", ]
near_miss_alasso <- summary_table["alasso_1nn_fp", ]

# Calculate as proportion of all false positives
near_miss_prop_ssvs <- near_miss_ssvs / pmax(summary_table["fp.ssvs", ], 1)
near_miss_prop_gets <- near_miss_gets / pmax(summary_table["fp.gets", ], 1)
near_miss_prop_alasso <- near_miss_alasso / pmax(summary_table["fp.alasso", ], 1)

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Proportion of False Positives",
  main = "Near Misses (±1 Period)",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

# Plot near-miss proportions
lines(log10(sd_numeric), near_miss_prop_ssvs, 
      col = colors$ssvs, lwd = settings$lwd.line)
points(log10(sd_numeric), near_miss_prop_ssvs, 
       col = colors$ssvs, pch = 16, cex = settings$cex.point)

lines(log10(sd_numeric), near_miss_prop_gets, 
      col = colors$gets, lwd = settings$lwd.line)
points(log10(sd_numeric), near_miss_prop_gets, 
       col = colors$gets, pch = 16, cex = settings$cex.point)

lines(log10(sd_numeric), near_miss_prop_alasso, 
      col = colors$alasso, lwd = settings$lwd.line)
points(log10(sd_numeric), near_miss_prop_alasso, 
       col = colors$alasso, pch = 16, cex = settings$cex.point)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("topright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       pch = 16,
       lwd = settings$lwd.line,
       bty = "n",
       cex = settings$cex.legend,
       pt.cex = settings$cex.point * 0.9)

mtext("F", side = 3, line = 1.5, at = par("usr")[1], cex = settings$cex.main, font = 2, adj = 0)

dev.off()

# ==============================================================================
# 9. CREATE INDIVIDUAL HIGH-RESOLUTION PLOTS
# ==============================================================================

# Individual plot for Precision
precision_file_individual <- sprintf("./Simulations/%s/%s_precision_gets-%s_bisam-%s_tau-%s.pdf", 
                                     date, 
                                     settings$suffix,
                                     gets_lvl,
                                     bisam_prior, 
                                     ifelse(tau == "", "auto", tau))

pdf(precision_file_individual, width = settings$pdf.width.single, height = settings$pdf.height.single)
par(mfrow = c(1, 1))
setup_plot()

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Precision",
  main = "Precision",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

lines(log10(metrics$SD_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point * 1.15)

lines(log10(metrics$SD_num[gets_idx]), metrics$Precision[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[gets_idx]), metrics$Precision[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point * 1.15)

lines(log10(metrics$SD_num[alasso_idx]), metrics$Precision[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[alasso_idx]), metrics$Precision[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point * 1.15)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       pch = 16,
       lwd = settings$lwd.line * 1.2,
       bty = "n",
       cex = settings$cex.legend * 1.15,
       pt.cex = settings$cex.point)

dev.off()

# Individual plot for Recall
recall_file_individual <- sprintf("./Simulations/%s/%s_recall_gets-%s_bisam-%s_tau-%s.pdf", 
                                  date, 
                                  settings$suffix,
                                  gets_lvl,
                                  bisam_prior, 
                                  ifelse(tau == "", "auto", tau))

pdf(recall_file_individual, width = settings$pdf.width.single, height = settings$pdf.height.single)
par(mfrow = c(1, 1))
setup_plot()

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Recall (Sensitivity)",
  main = "Recall",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

lines(log10(metrics$SD_num[ssvs_idx]), metrics$Recall[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[ssvs_idx]), metrics$Recall[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point * 1.15)

lines(log10(metrics$SD_num[gets_idx]), metrics$Recall[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[gets_idx]), metrics$Recall[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point * 1.15)

lines(log10(metrics$SD_num[alasso_idx]), metrics$Recall[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[alasso_idx]), metrics$Recall[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point * 1.15)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       pch = 16,
       lwd = settings$lwd.line * 1.2,
       bty = "n",
       cex = settings$cex.legend * 1.15,
       pt.cex = settings$cex.point)

dev.off()

# Individual plot for Near Misses
near_miss_file_individual <- sprintf("./Simulations/%s/%s_nearmiss_gets-%s_bisam-%s_tau-%s.pdf", 
                                     date, 
                                     settings$suffix,
                                     gets_lvl,
                                     bisam_prior, 
                                     ifelse(tau == "", "auto", tau))

pdf(near_miss_file_individual, width = settings$pdf.width.single, height = settings$pdf.height.single)
par(mfrow = c(1, 1))
setup_plot()

# Calculate near-miss rates
near_miss_ssvs <- summary_table["ssvs_1nn_fp", ]
near_miss_gets <- summary_table["gets_1nn_fp", ]
near_miss_alasso <- summary_table["alasso_1nn_fp", ]

near_miss_prop_ssvs <- near_miss_ssvs / pmax(summary_table["fp.ssvs", ], 1)
near_miss_prop_gets <- near_miss_gets / pmax(summary_table["fp.gets", ], 1)
near_miss_prop_alasso <- near_miss_alasso / pmax(summary_table["fp.alasso", ], 1)

plot(
  x = log10(sd_numeric), 
  y = NULL,
  xlim = range(log10(sd_numeric)),
  ylim = c(0, 1),
  xlab = "Threshold Level (SD, log scale)",
  ylab = "Proportion of False Positives",
  main = "Near Misses (±1 Period from True Break)",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

lines(log10(sd_numeric), near_miss_prop_ssvs, 
      col = colors$ssvs, lwd = settings$lwd.line * 1.2)
points(log10(sd_numeric), near_miss_prop_ssvs, 
       col = colors$ssvs, pch = 16, cex = settings$cex.point * 1.15)

lines(log10(sd_numeric), near_miss_prop_gets, 
      col = colors$gets, lwd = settings$lwd.line * 1.2)
points(log10(sd_numeric), near_miss_prop_gets, 
       col = colors$gets, pch = 16, cex = settings$cex.point * 1.15)

lines(log10(sd_numeric), near_miss_prop_alasso, 
      col = colors$alasso, lwd = settings$lwd.line * 1.2)
points(log10(sd_numeric), near_miss_prop_alasso, 
       col = colors$alasso, pch = 16, cex = settings$cex.point * 1.15)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("topright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       pch = 16,
       lwd = settings$lwd.line * 1.2,
       bty = "n",
       cex = settings$cex.legend * 1.15,
       pt.cex = settings$cex.point)

dev.off()

# Individual plot for F1 Score
pdf(f1_score_file, width = settings$pdf.width.single, height = settings$pdf.height.single)
par(mfrow = c(1, 1))
setup_plot()

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

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

lines(log10(metrics$SD_num[ssvs_idx]), metrics$F1[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[ssvs_idx]), metrics$F1[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point * 1.15)

lines(log10(metrics$SD_num[gets_idx]), metrics$F1[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[gets_idx]), metrics$F1[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point * 1.15)

lines(log10(metrics$SD_num[alasso_idx]), metrics$F1[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line * 1.2)
points(log10(metrics$SD_num[alasso_idx]), metrics$F1[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point * 1.15)

add_clean_axes(
  at_x = log10(sd_numeric), 
  labels_x = breaksize_labels,
  at_y = seq(0, 1, 0.2)
)

legend("bottomright", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       lwd = settings$lwd.line * 1.2,
       pch = 16,
       bty = "n",
       cex = settings$cex.legend * 1.15,
       pt.cex = settings$cex.point)

dev.off()

# ==============================================================================
# 10. CREATE ADDITIONAL COMPARISON PLOTS
# ==============================================================================

# Create a detailed comparison plot (Precision-Recall tradeoff)
pr_tradeoff_file <- sprintf("./Simulations/%s/%s_pr_tradeoff_gets-%s_bisam-%s_tau-%s.pdf", 
                            date, 
                            settings$suffix,
                            gets_lvl,
                            bisam_prior, 
                            ifelse(tau == "", "auto", tau))

pdf(pr_tradeoff_file, width = settings$pdf.width.single, height = settings$pdf.height.single)
par(mfrow = c(1, 1))
setup_plot()

plot(
  metrics$Recall[ssvs_idx], 
  metrics$Precision[ssvs_idx],
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "Recall (Sensitivity)",
  ylab = "Precision",
  main = "Precision-Recall Trade-off",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)
abline(v = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)
abline(0, 1, col = colors$gray, lty = 2, lwd = settings$lwd.axis)

# Plot trajectories
lines(metrics$Recall[ssvs_idx], metrics$Precision[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line * 1.2)
points(metrics$Recall[ssvs_idx], metrics$Precision[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point * 1.15)

lines(metrics$Recall[gets_idx], metrics$Precision[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line * 1.2)
points(metrics$Recall[gets_idx], metrics$Precision[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point * 1.15)

lines(metrics$Recall[alasso_idx], metrics$Precision[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line * 1.2)
points(metrics$Recall[alasso_idx], metrics$Precision[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point * 1.15)

# Add labels for breaksize levels
text(metrics$Recall[ssvs_idx], metrics$Precision[ssvs_idx], 
     labels = breaksize_labels, pos = 3, cex = 0.7 * settings$cex.legend, 
     col = colors$ssvs)

add_clean_axes(
  at_x = seq(0, 1, 0.2),
  at_y = seq(0, 1, 0.2)
)

legend("bottomleft", 
       legend = c("BISAM", "GETS", "ALASSO"),
       col = c(colors$ssvs, colors$gets, colors$alasso),
       lty = 1,
       lwd = settings$lwd.line * 1.2,
       pch = 16,
       bty = "n",
       cex = settings$cex.legend * 1.15,
       pt.cex = settings$cex.point)

dev.off()

cat("\n\n=================================================================\n")
cat("All plots have been saved successfully!\n")
cat("=================================================================\n\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================