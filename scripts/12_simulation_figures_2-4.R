################################################################################
# BISAM vs GETS vs ALASSO Comparison Analysis
#
# Description: Comprehensive comparison of BISAM, GETS, and Adaptive Lasso
#              methods across different threshold levels (SD) or break numbers (BN)
#
# Figures 2 & 3: vary sparse vs dense setting (SD breakpoints), 
# Figure  4    : vary number of breaks (BN breakpoints) at fixed break size
#
# Output: Performance metrics, plots, and summary statistics
################################################################################

rm(list = ls())

# ==============================================================================
#                         USER CONFIGURATION
# ==============================================================================

FIGURE <- 3          # 2 (sparse SD), 3 (dense SD), or 4 (BN)
DATE   <- "2026-01-23"
GETS_LVL <- 0.01

# Figure 4 only: which break size (in SD) to condition on
BREAKSIZE_SLCT <- 3

# ==============================================================================
# 1. SETUP
# ==============================================================================

library(dplyr)

gets_lvl    <- as.character(GETS_LVL)
bisam_prior <- "imom"
tau         <- "3.31744830051061"

# Mode derived from FIGURE choice
mode <- if (FIGURE %in% c(2, 3)) "SD" else if (FIGURE == 4) "BN" else stop("FIGURE must be 2, 3, or 4")

date      <- paste0(DATE, "_", mode)
setting   <- switch(as.character(FIGURE), "2" = "sparse", "3" = "dense", "4" = NULL)

data_path <- sprintf(
  "./output/simulation/%s/gets_bisam_comparison_gets-%s_bisam_prior-%s_tau-%s/",
  date, gets_lvl, bisam_prior, tau
)

file_list <- list.files(path = data_path, pattern = "\\.RDS$", full.names = TRUE)

# ==============================================================================
# 2. DATA PROCESSING
# ==============================================================================

combined_data <- lapply(file_list, function(file) {
  result_obj <- readRDS(file)
  mat <- if (is.list(result_obj) && "break_comparison" %in% names(result_obj)) {
    result_obj$break_comparison
  } else {
    result_obj
  }
  
  filename   <- basename(file)
  breaksize  <- as.numeric(sub(".*breaksize-([0-9.]+)SD_breaknumber.*",  "\\1", filename))
  replicate  <- as.integer(sub(".*_rep([0-9]+)\\.RDS", "\\1", filename))
  
  # breaknumber is character in SD mode, numeric in BN mode
  breaknumber <- if (mode == "SD")
    sub(".*breaknumber-([a-z]+)_rep.*", "\\1", filename)
  else
    as.numeric(sub(".*breaknumber-([0-9.]+)_rep.*", "\\1", filename))
  
  df             <- as.data.frame(mat)
  df$row_name    <- rownames(mat)
  df$breaksize   <- breaksize
  df$breaknumber <- breaknumber
  df$replicate   <- replicate
  df
}) |> bind_rows()

# Filter to the relevant subset
combined_data <- if (mode == "SD") {
  combined_data |> filter(breaknumber == setting)
} else {
  combined_data |> filter(breaksize == BREAKSIZE_SLCT)
}

# ==============================================================================
# 3. SUMMARY STATISTICS
# ==============================================================================

# The grouping variable differs by mode
group_var  <- if (mode == "SD") "breaksize" else "breaknumber"
group_vals <- combined_data[[group_var]] |> unique() |> sort()

summary_table <- sapply(group_vals, function(gv) {
  combined_data |>
    filter(.data[[group_var]] == gv) |>
    dplyr::select(-row_name, -breaksize, -breaknumber, -replicate) |>
    colSums()
})

colnames(summary_table) <- if (mode == "SD") {
  sprintf("%.1fSD", group_vals)
} else {
  sprintf("%.0f", group_vals)
}
cat("\n\nSummary Table:\n"); print(summary_table)

# ==============================================================================
# 4. PERFORMANCE METRICS
# ==============================================================================

group_labels <- colnames(summary_table)
group_numeric <- if (mode == "SD") {
  as.numeric(stringr::str_extract(group_labels, "\\d+\\.\\d(?=SD)"))
} else {
  as.numeric(group_labels)
}

metrics <- data.frame(
  GRP     = rep(group_labels,  3),
  GRP_num = rep(group_numeric, 3),
  Method  = rep(c("BISAM", "GETS", "ALASSO"), each = length(group_vals)),
  TP      = c(summary_table["tr.ssvs",  ], summary_table["tr.gets",  ], summary_table["tr.alasso",  ]),
  FP      = c(summary_table["fp.ssvs",  ], summary_table["fp.gets",  ], summary_table["fp.alasso",  ]),
  FN      = c(summary_table["fn.ssvs",  ], summary_table["fn.gets",  ], summary_table["fn.alasso",  ]),
  Total_Found = c(summary_table["ssvs", ], summary_table["gets",     ], summary_table["alasso",     ]),
  True_Total  = rep(summary_table["true", ], 3)
)

metrics <- within(metrics, {
  Precision <- TP / (TP + FP)
  Recall    <- TP / (TP + FN)
  F1        <- 2 * (Precision * Recall) / (Precision + Recall)
  FDR       <- FP / (TP + FP)
  Specificity <- 1 - FDR
})

metrics$Precision[is.nan(metrics$Precision)] <- 0
metrics$Recall[is.nan(metrics$Recall)]       <- 0
metrics$F1[is.nan(metrics$F1)]               <- 0

cat("\n\nPerformance Metrics by Method:\n"); print(metrics)

# ==============================================================================
# 5. COLORS & PLOT SETTINGS
# ==============================================================================

colors <- list(
  ssvs      = "#0072B2",
  gets      = "#D55E00",
  alasso    = "#009E73",
  gray      = "#808080",
  lightgray = "#E5E5E5",
  darkgray  = "#404040"
)

settings <- list(
  cex.lab    = 1.8,  cex.axis  = 1.6,  cex.main  = 2.0,  cex.legend = 1.4,
  lwd.line   = 4.5,  lwd.axis  = 2.5,  cex.point = 2.2,  lwd.grid   = 1.5,
  mar        = c(5.5, 6.0, 4.0, 2.0),
  mgp        = c(4.0, 1.0, 0),
  pdf.width  = 16,   pdf.height = 9
)

# ==============================================================================
# 6. HELPER FUNCTIONS
# ==============================================================================

setup_plot <- function(no_top_margin = FALSE) {
  par(
    mar     = if (no_top_margin) c(settings$mar[1], settings$mar[2], 2, settings$mar[4]) else settings$mar,
    mgp     = settings$mgp,
    las     = 1,
    cex.lab = settings$cex.lab,  cex.axis = settings$cex.axis,  cex.main = settings$cex.main,
    col.lab = colors$darkgray,   col.axis = colors$darkgray,    col.main = colors$darkgray,
    family  = "sans",  tcl = -0.3,  bty = "n"
  )
}

add_clean_axes <- function(at_x = NULL, labels_x = NULL, at_y = NULL) {
  axis(1, at = at_x, labels = labels_x, col = colors$gray, col.axis = colors$darkgray, lwd = settings$lwd.axis)
  axis(2, at = at_y,                    col = colors$gray, col.axis = colors$darkgray, lwd = settings$lwd.axis)
}

# Convenience: draw lines + points for one method
draw_series <- function(x, y, col, lty = 1, pch = 16, lwd_pt_factor = 1) {
  lines(x, y, col = col, lwd = settings$lwd.line, lty = lty)
  points(x, y, col = col, pch = pch, cex = settings$cex.point,
         lwd = settings$lwd.axis * lwd_pt_factor)
}

# ==============================================================================
# 7. AXIS TRANSFORMATION (log scale for SD mode, linear for BN mode)
# ==============================================================================

x_transform <- if (mode == "SD") log10 else identity
x_vals      <- x_transform(group_numeric)
x_label     <- if (mode == "SD") "Threshold Level (SD, log scale)" else "Number of breaks"

# ==============================================================================
# 8. CREATE PLOTS
# ==============================================================================

out_file <- sprintf("./output/simulation/figure_%d.pdf", FIGURE)
pdf(out_file, width = settings$pdf.width, height = settings$pdf.height)
par(mfrow = c(2, 2), oma = c(0, 0, 0, 0))

ssvs_idx   <- metrics$Method == "BISAM"
gets_idx   <- metrics$Method == "GETS"
alasso_idx <- metrics$Method == "ALASSO"

x_ssvs   <- x_transform(metrics$GRP_num[ssvs_idx])
x_gets   <- x_transform(metrics$GRP_num[gets_idx])
x_alasso <- x_transform(metrics$GRP_num[alasso_idx])

# --- helper to open a blank panel ---
open_panel <- function(ylab, main, ylim = c(0, 1)) {
  setup_plot()
  plot(x_vals, NULL,
       xlim = range(x_vals), ylim = ylim,
       xlab = x_label, ylab = ylab, main = main,
       type = "n", axes = FALSE, xaxs = "i", yaxs = "i")
  abline(h = seq(ylim[1], ylim[2], length.out = 6),
         col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)
}

std_legend <- function(position = "bottomright") {
  legend(position,
         legend = c("BISAM", "GETS", "ALASSO"),
         col    = c(colors$ssvs, colors$gets, colors$alasso),
         lty = 1, pch = 16, lwd = settings$lwd.line, bty = "n",
         cex = settings$cex.legend, pt.cex = settings$cex.point * 0.9)
}

panel_label <- function(label) {
  mtext(label, side = 3, line = 1.5, at = par("usr")[1],
        cex = settings$cex.main, font = 2, adj = 0)
}

# ------------------------------------------------------------------------------
# Panel A: Precision
# ------------------------------------------------------------------------------
open_panel("Precision", "Precision")

draw_series(x_ssvs,   metrics$Precision[ssvs_idx],   colors$ssvs)
draw_series(x_gets,   metrics$Precision[gets_idx],   colors$gets)
draw_series(x_alasso, metrics$Precision[alasso_idx], colors$alasso)

add_clean_axes(at_x = x_vals, labels_x = group_labels, at_y = seq(0, 1, 0.2))
std_legend("bottomright")
panel_label("A")

# ------------------------------------------------------------------------------
# Panel B: F1 Score
# ------------------------------------------------------------------------------
open_panel("F1 Score", "F1 Score")

draw_series(x_ssvs,   metrics$F1[ssvs_idx],   colors$ssvs)
draw_series(x_gets,   metrics$F1[gets_idx],   colors$gets)
draw_series(x_alasso, metrics$F1[alasso_idx], colors$alasso)

add_clean_axes(at_x = x_vals, labels_x = group_labels, at_y = seq(0, 1, 0.2))
std_legend("bottomright")
panel_label("B")

# ------------------------------------------------------------------------------
# Panel C: True and False Positive Rates
# ------------------------------------------------------------------------------
tp_ssvs   <- summary_table["tr.ssvs",  ] / summary_table["true", ]
tp_gets   <- summary_table["tr.gets",  ] / summary_table["true", ]
tp_alasso <- summary_table["tr.alasso",] / summary_table["true", ]
fp_ssvs   <- summary_table["fp.ssvs",  ] / summary_table["true", ]
fp_gets   <- summary_table["fp.gets",  ] / summary_table["true", ]
fp_alasso <- summary_table["fp.alasso",] / summary_table["true", ]

ylim_max  <- max(c(tp_ssvs, tp_gets, fp_ssvs, fp_gets)) * 1.1
y_ticks   <- pretty(c(0, ylim_max), n = 5)

open_panel("Rate (relative to true breaks)", "True and False Positive Rates",
           ylim = c(0, ylim_max))
# abline(h = y_ticks, col = colors$lightgray, lty = 1, lwd = settings$lwd.grid)
abline(h = 1,       col = colors$gray,      lty = 2, lwd = settings$lwd.axis)

draw_series(x_vals, tp_ssvs,   colors$ssvs,   lty = 1, pch = 16)
draw_series(x_vals, tp_gets,   colors$gets,   lty = 1, pch = 16)
draw_series(x_vals, tp_alasso, colors$alasso, lty = 1, pch = 16)
draw_series(x_vals, fp_ssvs,   colors$ssvs,   lty = 2, pch = 1, lwd_pt_factor = 0.8)
draw_series(x_vals, fp_gets,   colors$gets,   lty = 2, pch = 1, lwd_pt_factor = 0.8)
draw_series(x_vals, fp_alasso, colors$alasso, lty = 2, pch = 1, lwd_pt_factor = 0.8)

add_clean_axes(at_x = x_vals, labels_x = group_labels, at_y = y_ticks)

legend("topright",
       legend = c("BISAM (TP)", "GETS (TP)", "ALASSO (TP)",
                  "BISAM (FP)", "GETS (FP)", "ALASSO (FP)"),
       col    = rep(c(colors$ssvs, colors$gets, colors$alasso), 2),
       lty    = c(1, 1, 1, 2, 2, 2),
       pch    = c(16, 16, 16, 1, 1, 1),
       lwd    = settings$lwd.line, bty = "n",
       cex    = settings$cex.legend * 0.85, pt.cex = settings$cex.point * 0.85, ncol = 1)

panel_label("C")

# ------------------------------------------------------------------------------
# Panel D: Near Misses (±1 Period)
# ------------------------------------------------------------------------------
nm_ssvs   <- summary_table["ssvs_1nn_fp",  ] / pmax(summary_table["fp.ssvs",  ], 1)
nm_gets   <- summary_table["gets_1nn_fp",  ] / pmax(summary_table["fp.gets",  ], 1)
nm_alasso <- summary_table["alasso_1nn_fp",] / pmax(summary_table["fp.alasso",], 1)

open_panel("Proportion of False Positives", "Near Misses (\u00b11 Period)")

draw_series(x_vals, nm_ssvs,   colors$ssvs)
draw_series(x_vals, nm_gets,   colors$gets)
draw_series(x_vals, nm_alasso, colors$alasso)

add_clean_axes(at_x = x_vals, labels_x = group_labels, at_y = seq(0, 1, 0.2))
std_legend("topright")
panel_label("D")

dev.off()
cat("\nSaved:", out_file, "\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================