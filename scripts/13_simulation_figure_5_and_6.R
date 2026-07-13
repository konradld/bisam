################################################################################
# BISAM vs GETS vs ALASSO Comparison Analysis
# 
# Description: Comprehensive comparison of BISAM, GETS, and Adaptive Lasso 
#              methods across different N and T
# 
# Output: Performance metrics, plots, and summary statistics
################################################################################

# Clear workspace
rm(list = ls())

#===============================================================================
#           Which Figure to replicate (see paper for numbering)
# FIGURE <- 5
DATE <- "2026-07-07"
#===============================================================================

# ==============================================================================
# 1. SETUP AND DATA LOADING
# ==============================================================================

# Load required libraries (only for data manipulation)
library(dplyr)

date <- paste0(DATE, "_no-FEs_relative-breaks_NT")
gets_lvl <- "0.01"
bisam_prior <- "imom"
tau <- "3.31744830051061"

# Set up paths
data_path <- sprintf("./output/simulation/%s/gets_bisam_comparison_gets-%s_bisam_prior-%s_tau-%s/", 
                     date, 
                     gets_lvl,
                     bisam_prior,
                     tau)

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
  breaksize <- as.numeric(sub(".*breaksize-([0-9.]+)SD_breaknumber.*", "\\1", filename))
  Nt <- as.numeric(sub(".*Nt-([0-9.]+)_rep.*", "\\1", filename))
  Ni <- as.numeric(sub(".*Ni-([0-9.]+)_Nt.*", "\\1", filename))
  breaknumber <- sub(".*breaknumber-([a-z]+)_Ni.*", "\\1", filename)
  replicate <- as.integer(sub(".*_rep([0-9]+)\\.RDS", "\\1", filename))
  
  
  # Convert to data frame with identifiers
  df <- as.data.frame(mat)
  df$row_name <- rownames(mat)
  df$breaksize <- breaksize
  df$breaknumber <- breaknumber
  df$Nt <- Nt
  df$Ni <- Ni
  df$replicate <- replicate
  
  return(df)
}) |> 
  bind_rows()

### subsetting here

FIGURE <- 6

setting <- if(FIGURE == 5) "sparse" else if(FIGURE == 6) "dense"

# break SD to consider
breaksize_slct <- 3
Nt_slct <- seq(from = 20, to = 100, by = 5)
Ni_slct <- 10 # seq(from = 5, to = 30, by = 5)



filtered_data <- combined_data |> 
  filter(Nt %in% Nt_slct,
         Ni %in% Ni_slct, 
         breaknumber == setting, 
         breaksize == breaksize_slct)

# ==============================================================================
# 3. SUMMARY STATISTICS
# ==============================================================================

# Get unique breaknumbers and sort them
Nt_vals <- filtered_data$Nt |> unique() |> sort()
Ni_vals <- filtered_data$Ni |> unique() |> sort()

# Calculate summary table for each breaknumber
summary_table <- sapply(Nt_vals, function(nn) {
  data_subset <- filtered_data |>
    filter(Nt == nn)
  
  # Sum across all replicates
  summary_stats <- data_subset |>
    dplyr::select(-row_name, -breaksize, -breaknumber, -Nt, -Ni, -replicate) |> 
    colSums()
  
  return(summary_stats)
})

# Set column names
colnames(summary_table) <- sprintf("%.0f", Nt_vals)

# Display summary table
cat("\n\nSummary Table:\n")
print(summary_table)

# ==============================================================================
# 4. PERFORMANCE METRICS CALCULATION
# ==============================================================================

# Define threshold levels
Nt_labels <- colnames(summary_table)
Nt_numeric <- as.numeric(Nt_labels)

# Initialize metrics data frame for all three methods
metrics <- data.frame(
  NT = rep(Nt_labels, 3),
  Nt_num = rep(Nt_numeric, 3),
  Method = rep(c("BISAM", "GETS", "ALASSO"), each = length(Nt_vals)),
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
  suffix = "paper"
)

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
multi_panel_file <- sprintf("./output/simulation/figure_%d_no-FEs-relative-breaks_Ni-%d.pdf", FIGURE, Ni_slct)

# Open PDF for multi-panel plot
pdf(multi_panel_file, width = settings$pdf.width, height = settings$pdf.height)

# Use 2x2 layout
par(mfrow = c(2, 2), oma = c(0, 0, 0, 0))

# Extract indices for each method
ssvs_idx <- metrics$Method == "BISAM"
gets_idx <- metrics$Method == "GETS"
alasso_idx <- metrics$Method == "ALASSO"

# ------------------------------------------------------------------------------
# Plot 1: Precision
# ------------------------------------------------------------------------------
setup_plot()

plot(
  x = (Nt_numeric), 
  y = NULL,
  xlim = range((Nt_numeric)),
  ylim = c(0, 1),
  xlab = "Number of time periods",
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
lines((metrics$Nt_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line)
points((metrics$Nt_num[ssvs_idx]), metrics$Precision[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point)

lines((metrics$Nt_num[gets_idx]), metrics$Precision[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line)
points((metrics$Nt_num[gets_idx]), metrics$Precision[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point)

lines((metrics$Nt_num[alasso_idx]), metrics$Precision[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line)
points((metrics$Nt_num[alasso_idx]), metrics$Precision[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point)

add_clean_axes(
  at_x = (Nt_numeric), 
  labels_x = Nt_labels,
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
# Plot 2: F1 Score
# ------------------------------------------------------------------------------
setup_plot()

plot(
  x = (Nt_numeric), 
  y = NULL,
  xlim = range((Nt_numeric)),
  ylim = c(0, 1),
  xlab = "Number of time periods",
  ylab = "F1 Score",
  main = "F1 Score",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

lines((metrics$Nt_num[ssvs_idx]), metrics$F1[ssvs_idx], 
      col = colors$ssvs, lwd = settings$lwd.line)
points((metrics$Nt_num[ssvs_idx]), metrics$F1[ssvs_idx], 
       col = colors$ssvs, pch = 16, cex = settings$cex.point)

lines((metrics$Nt_num[gets_idx]), metrics$F1[gets_idx], 
      col = colors$gets, lwd = settings$lwd.line)
points((metrics$Nt_num[gets_idx]), metrics$F1[gets_idx], 
       col = colors$gets, pch = 16, cex = settings$cex.point)

lines((metrics$Nt_num[alasso_idx]), metrics$F1[alasso_idx], 
      col = colors$alasso, lwd = settings$lwd.line)
points((metrics$Nt_num[alasso_idx]), metrics$F1[alasso_idx], 
       col = colors$alasso, pch = 16, cex = settings$cex.point)

add_clean_axes(
  at_x = (Nt_numeric), 
  labels_x = Nt_labels,
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

mtext("B", side = 3, line = 1.5, at = par("usr")[1], cex = settings$cex.main, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 3: True Positive and False Positive Rates
# ------------------------------------------------------------------------------
setup_plot()

# Calculate rates
tp_rate_ssvs <- summary_table["tr.ssvs", ] / summary_table["true", ]
tp_rate_gets <- summary_table["tr.gets", ] / summary_table["true", ]
tp_rate_alasso <- summary_table["tr.alasso", ] / summary_table["true", ]

fp_rate_ssvs <- summary_table["fp.ssvs", ] / summary_table["true", ]
fp_rate_gets <- summary_table["fp.gets", ] / summary_table["true", ]
fp_rate_alasso <- summary_table["fp.alasso", ] / summary_table["true", ]

# Determine y-axis range
all_rates <- c(tp_rate_ssvs, tp_rate_gets, #tp_rate_alasso,
               fp_rate_ssvs, fp_rate_gets) #, fp_rate_alasso)
ylim_max <- max(all_rates) * 1.1

plot(
  x = (Nt_numeric), 
  y = NULL,
  xlim = range((Nt_numeric)),
  ylim = c(0, ylim_max),
  xlab = "Number of time periods",
  ylab = "Rate (relative to true breaks)",
  main = "True and False Positive Rates",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

# Grid
y_ticks <- pretty(c(0, ylim_max), n = 5)
abline(h = y_ticks, col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)
abline(h = 1, col = colors$gray, lty = 2, lwd = settings$lwd.axis)

# Plot True Positive Rates (solid lines)
lines((Nt_numeric), tp_rate_ssvs, 
      col = colors$ssvs, lwd = settings$lwd.line, lty = 1)
points((Nt_numeric), tp_rate_ssvs, 
       col = colors$ssvs, pch = 16, cex = settings$cex.point)

lines((Nt_numeric), tp_rate_gets, 
      col = colors$gets, lwd = settings$lwd.line, lty = 1)
points((Nt_numeric), tp_rate_gets, 
       col = colors$gets, pch = 16, cex = settings$cex.point)

lines((Nt_numeric), tp_rate_alasso, 
      col = colors$alasso, lwd = settings$lwd.line, lty = 1)
points((Nt_numeric), tp_rate_alasso, 
       col = colors$alasso, pch = 16, cex = settings$cex.point)

# Plot False Positive Rates (dashed lines)
lines((Nt_numeric), fp_rate_ssvs, 
      col = colors$ssvs, lwd = settings$lwd.line, lty = 2)
points((Nt_numeric), fp_rate_ssvs, 
       col = colors$ssvs, pch = 1, cex = settings$cex.point, lwd = settings$lwd.axis * 0.8)

lines((Nt_numeric), fp_rate_gets, 
      col = colors$gets, lwd = settings$lwd.line, lty = 2)
points((Nt_numeric), fp_rate_gets, 
       col = colors$gets, pch = 1, cex = settings$cex.point, lwd = settings$lwd.axis * 0.8)

lines((Nt_numeric), fp_rate_alasso, 
      col = colors$alasso, lwd = settings$lwd.line, lty = 2)
points((Nt_numeric), fp_rate_alasso, 
       col = colors$alasso, pch = 1, cex = settings$cex.point, lwd = settings$lwd.axis * 0.8)

add_clean_axes(
  at_x = (Nt_numeric), 
  labels_x = Nt_labels,
  at_y = y_ticks
)

legend("topright", 
       legend = c("BISAM (TP)", "GETS (TP)", "ALASSO (TP)",
                  "BISAM (FP)", "GETS (FP)", "ALASSO (FP)"),
       col = c(colors$ssvs, colors$gets, colors$alasso,
               colors$ssvs, colors$gets, colors$alasso),
       lty = c(1, 1, 1, 2, 2, 2),
       pch = c(16, 16, 16, 1, 1, 1),
       lwd = settings$lwd.line,
       bty = "n",
       cex = settings$cex.legend * 0.85,
       pt.cex = settings$cex.point * 0.85,
       ncol = 1)

mtext("C", side = 3, line = 1.5, at = par("usr")[1], cex = settings$cex.main, font = 2, adj = 0)

# ------------------------------------------------------------------------------
# Plot 4: Near Misses (1-Period Neighbor False Positives)
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
  x = (Nt_numeric), 
  y = NULL,
  xlim = range((Nt_numeric)),
  ylim = c(0, 1),
  xlab = "Number of time periods",
  ylab = "Proportion of False Positives",
  main = "Near Misses (±1 Period)",
  type = "n",
  axes = FALSE,
  xaxs = "i",
  yaxs = "i"
)

abline(h = seq(0, 1, 0.2), col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)

# Plot near-miss proportions
lines(Nt_numeric, near_miss_prop_ssvs, 
      col = colors$ssvs, lwd = settings$lwd.line)
points(Nt_numeric, near_miss_prop_ssvs, 
       col = colors$ssvs, pch = 16, cex = settings$cex.point)

lines(Nt_numeric, near_miss_prop_gets, 
      col = colors$gets, lwd = settings$lwd.line)
points(Nt_numeric, near_miss_prop_gets, 
       col = colors$gets, pch = 16, cex = settings$cex.point)

lines(Nt_numeric, near_miss_prop_alasso, 
      col = colors$alasso, lwd = settings$lwd.line)
points(Nt_numeric, near_miss_prop_alasso, 
       col = colors$alasso, pch = 16, cex = settings$cex.point)

add_clean_axes(
  at_x = Nt_numeric, 
  labels_x = Nt_labels,
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

mtext("D", side = 3, line = 1.5, at = par("usr")[1], cex = settings$cex.main, font = 2, adj = 0)

dev.off()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
