# rm(list = ls())

library(dplyr)

# ==============================================================================
# SETUP AND DATA LOADING
# ==============================================================================
date <- "2026-02-27"
check_outl <- "TRUE"
out_scale <- 10
tau <- "1.92072941034706" #  or 3.31744830051061
prior <- "imom"
c0 <- C0 <- "auto"
modprior <- "bern"
v0 <- "0.5"

bisam_path <- sprintf("./output/emissions/%s/ssvs_checkOutlier-%s_outscale-%s_tau-%s_prior-%s_c0-%s_C0-%s_modprior-%s_v0-%s.RDS", 
                      date, check_outl, out_scale, tau, prior, c0, C0, modprior, v0) 

bisam_results <- readRDS(bisam_path)

plot(bisam_results$coefs$iis, main = paste0("OUTLIER SCALE = ", out_scale), ylab = "IIS PIPs")

gets_results_05 <- readRDS(sprintf("./output/emissions/%s/gets_0.05.RDS", 
                                   date))

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

#' Adds shaded rectangular areas to a plot with intensity based on break detection
#' @param breaks_to_shade Data frame as outputted by pip_window_fun.R
#' @param ylim Y-axis limits for shading
#' @param colorin_variable Which variable to use for measuring intensity (e.g purity or break size)
#' @param col Two-element vector for positive/negative break colors
add_break_shading_gradient <- function(breaks_to_shade,
                                       ylim = c(0, 1),
                                       colorin_variable = rep(1, nrow(breaks_to_shade)),
                                       color_variable = rep(1, nrow(breaks_to_shade)),
                                       col = c("red", "green")) {
  
  # Normalize color variable to [0, 1] range
  if(colorin_variable != "break_length") {
    abs_col_var <- abs(colorin_variable)
    abs_col_val <- (abs_col_var - min(abs_col_var)) / max(abs_col_var)
  }

  col_val <- sign(color_variable)
  
  # Calculate break duration
  break_length <- as.numeric(breaks_to_shade$end_time) - 
    as.numeric(breaks_to_shade$start_time) + 1
  
  # Draw shaded rectangles for each break
  for (i in 1:nrow(breaks_to_shade)) {
    # Select color based on sign of color variable
    colour <- if (col_val[i] > 0) col[1] else col[2]
    
    # Alpha transparency proportional to detection value 
    if(colorin_variable == "break_length") {
      alpha <- 1 - break_length[i] / 5
    } else {
      alpha <- abs_col_val[i]
    }
    
    if(i > 1) {
      last_range <- breaks_to_shade$start_time[i-1]:breaks_to_shade$end_time[i-1]
      if(breaks_to_shade$unit[i - 1] == breaks_to_shade$unit[i] & 
           breaks_to_shade$start_time[i] %in% last_range) {
        x_left_i <- breaks_to_shade$position_idx[i] + breaks_to_shade$unit_idx[i] + 
          which(breaks_to_shade$start_time[i] %in% last_range) + 1
      } else {
        x_left_i <- breaks_to_shade$position_idx[i] + breaks_to_shade$unit_idx[i]
      }
    } else {
      x_left_i <- breaks_to_shade$position_idx[i] + breaks_to_shade$unit_idx[i]
    }

    
    # Draw rectangle
    rect(
      xleft = x_left_i - 1  - 0.5,
      xright = breaks_to_shade$position_idx[i] + breaks_to_shade$unit_idx[i] - 1 + break_length[i] - 0.5,
      ybottom = ylim[1],
      ytop = ylim[2],
      col = adjustcolor(colour, alpha.f = alpha),
      border = NA
    )
  }
}

# ==============================================================================
# PLOT FOR OMEGAS AND FITTED (PANEL A + B)
# ==============================================================================

library(stringr)
source("./R/pip_window_fun.R")

# Load sector-specific results
bisam_res_sector <- bisam_results
mod <- bisam_res_sector
bisam_coefs <- mod$coef_list

# Get GETS results
gets_coefs05 <- gets_results_05$isatpanel.result$mean.results
lux_man <- data.frame(coef = -0.15, std.error = 0.05, "t-stat" = 4, "p-value" = 0.0001)
names(lux_man)[3] <- "t-stat"
names(lux_man)[4] <- "p-value"
rownames(lux_man) <- "fesisLuxembourg.2015"
gets_coefs <- rbind(gets_coefs05, lux_man)
rownames(gets_coefs) <- gsub("UnitedKingdom", "United Kingdom", rownames(gets_coefs))


# Parameters
PIP_THRESHOLD <- 0.5

# Calculate breaks
win1_pips <- pip_window(mod, win_size = 1, op = ">=", pip_threshold = PIP_THRESHOLD)
win2_pips <- pip_window(mod, win_size = 2, op = ">=", pip_threshold = 1 - PIP_THRESHOLD^2)
win3_pips <- pip_window(mod, win_size = 3, op = ">=", pip_threshold = 1 - PIP_THRESHOLD^3)

# Color palette
COLORS <- list(
  main = "#0072B2",
  step = "black",
  outl_marker = "#D5B60A",
  fit = "#1B9E77",
  grid = "gray70",
  positive_break = "#009E73",
  negative_break = "#D55E00",
  gets01 = "black",
  gets05 = "goldenrod"
)


pdf(sprintf("./output/emissions/%s/multi_checkOutlier-%s_outscale-%s_tau-%s_prior-%s.pdf",
            date, check_outl, out_scale, tau, prior),
    width = 18, height = 8, onefile = TRUE)

# Break up in three blocks

c_list <- list(c("Germany", "Spain", "France", "United Kingdom", "Italy"),
               c("Austria", "Belgium", "Netherlands", "Portugal", "Sweden"),
               c("Denmark", "Finland", "Greece", "Ireland", "Luxembourg"))


for(cc in seq_len(length(c_list))) {

  cN <- c_list[[cc]]
  # Extract sample dimensions
  t_mod <- mod$meta$sample['t']
  n_mod <- length(cN)
  N_mod <- mod$meta$sample['N']

  # Setup plot layout
  par(
    mfrow = c(2, 1),
    mar = c(1, 4.5, 1.5, 1),
    oma = c(1, 0, 2, 0),
    cex.axis = 1,
    cex.lab = 1.25,
    las = 1
  )
  
  # ------------------------------------------------------------------------------
  # PANEL A: OMEGA COEFFICIENTS
  # ------------------------------------------------------------------------------

  # Create omega vector with NAs at country boundaries
  omega_with_breaks <- mod$coefs$omega
  omega_with_breaks <- omega_with_breaks[which(grepl(paste0(cN, collapse = "|"), names(omega_with_breaks)))]
  country_boundaries_a <- (1:(n_mod - 1)) * (t_mod - 3)

  omega_plot <- rep(NA, length(omega_with_breaks) + (n_mod - 1))
  current_pos <- 1
  for(i in 1:n_mod) {
    start_idx <- if(i == 1) 1 else (i-1) * (t_mod - 3) + 1
    end_idx <- i * (t_mod - 3)
    segment_length <- end_idx - start_idx + 1

    omega_plot[current_pos:(current_pos + segment_length - 1)] <- omega_with_breaks[start_idx:end_idx]
    current_pos <- current_pos + segment_length + 1
  }

  x_coords_a <- 1:length(omega_plot)

  # Plot omega coefficients
  plot(
    x_coords_a, omega_plot,
    type = "l",
    col = COLORS$main,
    lwd = 2,
    ylim = c(0, 1),
    xlab = "",
    ylab = "PIP",
    main = "",
    xaxt = "n",
    frame.plot = FALSE
  )

  # # Add uncertainty band
  # uncertainty <- apply(mod$draws$omega, 2, quantile, c(0.25, 0.75))
  # uncertainty <- uncertainty[,which(grepl(paste0(cN, collapse = "|"), colnames(uncertainty)))]
  # for(i in 1:n_mod) {
  #   start_idx <- if(i == 1) 1 else (i-1) * (t_mod - 3) + 1
  #   end_idx <- i * (t_mod - 3)
  #   plot_start <- start_idx + (i - 1)
  #   plot_end <- end_idx + (i - 1)
  #   x_seg <- plot_start:plot_end
  #
  #   polygon(
  #     c(x_seg, rev(x_seg)),
  #     c(uncertainty[1, start_idx:end_idx], rev(uncertainty[2, start_idx:end_idx])),
  #     col = rgb(0.68, 0.85, 0.9, 0.3),
  #     border = NA
  #   )
  # }

  cN_idx <- data.frame(unit = cN, unit_idx = seq_len(length(cN)))
  
  # add break shading
  
  if(!is.null(nrow(win1_pips))) {
    win1_pips_temp <- win1_pips |> filter(unit %in% cN) |>
      select(-unit_idx) |> left_join(cN_idx, by ="unit") |>
      relocate(unit_idx, .after = purity) |>
      mutate(position_idx = (unit_idx - 1) * (t_mod - 3) + time_idx)
    
    if(nrow(win1_pips_temp) > 0) {
      add_break_shading_gradient(
        breaks_to_shade = win1_pips_temp,
        colorin_variable = "break_length",
        color_variable = win1_pips_temp$MAP,
        ylim = c(0, 1),
        col = c(COLORS$positive_break, COLORS$negative_break)
      )
    }
  }

  if(!is.null(nrow(win2_pips))) {
    win2_pips_temp <- win2_pips |> filter(unit %in% cN) |>
      select(-unit_idx) |> left_join(cN_idx, by ="unit") |>
      relocate(unit_idx, .after = purity) |>
      mutate(position_idx = (unit_idx - 1) * (t_mod - 3) + time_idx)
    
    if(nrow(win2_pips_temp) > 0) {
      add_break_shading_gradient(
        breaks_to_shade = win2_pips_temp,
        colorin_variable = "break_length",
        color_variable = win2_pips_temp$MAP,
        ylim = c(0, 1),
        col = c(COLORS$positive_break, COLORS$negative_break)
      )
    }
  }
  

  if(!is.null(nrow(win3_pips))) {
    win3_pips_temp <- win3_pips |> filter(unit %in% cN) |>
      select(-unit_idx) |> left_join(cN_idx, by ="unit") |>
      relocate(unit_idx, .after = purity) |>
      mutate(position_idx = (unit_idx - 1) * (t_mod - 3) + time_idx)
    
    if(nrow(win3_pips_temp) > 0) {
      add_break_shading_gradient(
        breaks_to_shade = win3_pips_temp,
        colorin_variable = "break_length",
        color_variable = win3_pips_temp$MAP,
        ylim = c(0, 1),
        col = c(COLORS$positive_break, COLORS$negative_break)
      )
    }
  }
  
  # Redraw main line
  lines(x_coords_a, omega_plot, col = COLORS$main, lwd = 2)

  # Add reference lines
  abline(h = c(0, 0.25, 0.5, 0.75, 1), lty = 3, col = COLORS$grid, lwd = 0.8)

  for(i in 1:(n_mod-1)) {
    abline(v = i * (t_mod - 3) + i, lty = 2, col = COLORS$grid, lwd = 0.8)
  }

  # Mark detected breaks
  break_indices_orig <- (1:(n_mod * (t_mod - 3)))[omega_with_breaks >= PIP_THRESHOLD]
  
  if(length(break_indices_orig) > 0) {
    break_indices_plot <- break_indices_orig + sapply(break_indices_orig, function(x) {
      sum(country_boundaries_a < x)
    })
    segments(
      x0 = break_indices_plot, y0 = 0,
      x1 = break_indices_plot, y1 = 1,
      col = COLORS$step, lty = 2, lwd = 1.5
    )
  }
  
  # Add GETS detections
  gets_indices <- (1:(n_mod * (t_mod - 3)))[
    str_extract(names(omega_with_breaks), "(?<=sis\\.).+") %in%
      str_extract(rownames(gets_coefs), "(?<=fesis).+")
  ]
  
  if (length(gets_indices) > 0) { # maybe differentiate cross color by sign of break?
    gets_indices_plot <- gets_indices + sapply(gets_indices, function(x) {
      sum(country_boundaries_a < x)
    })
    points(x = gets_indices_plot, y = omega_plot[gets_indices_plot],
           pch = 4, col = COLORS$gets01, lwd = 2.5, cex = 1.5)
  }
  
  # gets_indices05 <- (1:(n_mod * (t_mod - 3)))[
  #   str_extract(names(omega_with_breaks), "(?<=sis\\.).+") %in%
  #     str_extract(rownames(gets_coefs05), "(?<=fesis).+")
  # ]
  # if (length(gets_indices05) > 0) {
  #   gets_indices05_plot <- gets_indices05 + sapply(gets_indices05, function(x) {
  #     sum(country_boundaries_a < x)
  #   })
  #   points(x = gets_indices05_plot, y = omega_plot[gets_indices05_plot],
  #          pch = 4, col = COLORS$gets01, lwd = 2.5, cex = 1.2)
  # }

  # Add country labels
  midpoints_a <- sapply(1:n_mod, function(i) {
    start_idx <- if(i == 1) 1 else (i-1) * (t_mod - 3) + 1
    end_idx <- i * (t_mod - 3)
    mean(c(start_idx + (i-1), end_idx + (i-1)))
  })
  text(
    x = midpoints_a,
    y = 0,
    labels = cN,
    pos = 1,
    xpd = TRUE,
    cex = 1.25
  )

  # Add legend
  legend(
    "topright",
    # inset=c(0,-0.1),
    # ifelse(cc < 3, "topright", "topleft"),
    legend = c("PIP", "Detected Breaks", "Positive Break Int.",
               "Negative Break Int.", "Gets Detections"),
    col = c(COLORS$main, COLORS$step, NA, NA, COLORS$gets01),
    fill = c(NA, NA, COLORS$positive_break, COLORS$negative_break, NA),
    border = c(NA, NA, COLORS$positive_break, COLORS$negative_break, NA),
    lty = c(1, 2, NA, NA, NA),
    pch = c(NA, NA, NA, NA, 4),
    lwd = c(2, 2, NA, NA, 2),
    cex = 1.25,
    bg = adjustcolor("white", alpha.f = 0.5)# , horiz = T
  )

  # Add panel title
  mtext(
    side = 3, line = 0.5,
    text = "Panel A: PIP estimates",
    adj = 0, font = 2, cex = 1.5
  )

  # ------------------------------------------------------------------------------
  # PANEL B: FITTED VALUES
  # ------------------------------------------------------------------------------

  par(mar = c(3, 4.5, 0.5, 1))

  temp_y <- mod$data$original_data |> filter(country %in% cN) |> select(ltransport.emissions) |> as.matrix()
  temp_idx <- which(mod$data$original_data$country %in% cN)
  # Create vectors with NAs at country boundaries
  country_boundaries_b <- (1:(n_mod-1)) * t_mod
  n_total_b <- length(temp_y) + (n_mod - 1)

  y_plot <- rep(NA, n_total_b)
  current_pos <- 1
  for(i in 1:n_mod) {
    start_idx <- if(i == 1) 1 else (i-1) * t_mod + 1
    end_idx <- i * t_mod
    segment_length <- end_idx - start_idx + 1

    y_plot[current_pos:(current_pos + segment_length - 1)] <-
      temp_y[start_idx:end_idx]
    current_pos <- current_pos + segment_length + 1
  }

  x_coords_b <- 1:n_total_b

  # Plot observed response
  plot(
    x_coords_b, y_plot,
    cex = 0.6,
    pch = 19,
    col = adjustcolor(COLORS$main, alpha.f = 0.4),
    xlab = "",
    ylab = "Response (y)",
    main = "",
    xaxt = "n",
    frame.plot = FALSE
  )

  # Add vertical reference lines
  for(i in 1:(n_mod-1)) {
    abline(v = i * t_mod + i, lty = 2, col = COLORS$grid, lwd = 0.8)
  }

  # Calculate fitted values
  y_fitted_orig <- as.matrix(mod$data$X %*% mod$coefs$beta + mod$data$Z %*% mod$coefs$sis)[temp_idx, ]
  if(length(mod$coefs$sigma2) > 1) {
    sigma <- sqrt(mod$coefs$sigma2[which(grepl(paste0(cN, collapse = "|"), names(mod$coefs$sigma2)))])
  } else {
    sigma <- rep(sqrt(mod$coefs$sigma2), n_mod)
  }


  y_fitted_plot <- rep(NA, n_total_b)
  upper_plot <- rep(NA, n_total_b)
  lower_plot <- rep(NA, n_total_b)

  current_pos <- 1
  for(i in 1:n_mod) {
    start_idx <- if(i == 1) 1 else (i-1) * t_mod + 1
    end_idx <- i * t_mod
    segment_length <- end_idx - start_idx + 1

    y_fitted_plot[current_pos:(current_pos + segment_length - 1)] <-
      y_fitted_orig[start_idx:end_idx]
    upper_plot[current_pos:(current_pos + segment_length - 1)] <-
      y_fitted_orig[start_idx:end_idx] + 2 * sigma[i]
    lower_plot[current_pos:(current_pos + segment_length - 1)] <-
      y_fitted_orig[start_idx:end_idx] - 2 * sigma[i]

    current_pos <- current_pos + segment_length + 1
  }
  
  # Add credible interval
  for(i in 1:n_mod) {
    start_idx <- if(i == 1) 1 else (i-1) * t_mod + 1
    end_idx <- i * t_mod
    plot_start <- start_idx + (i - 1)
    plot_end <- end_idx + (i - 1)
    x_seg <- plot_start:plot_end

    polygon(
      c(x_seg, rev(x_seg)),
      c(lower_plot[x_seg], rev(upper_plot[x_seg])),
      col = adjustcolor(COLORS$fit, alpha.f = 0.2),
      border = NA
    )
  }

  # Draw fitted line
  lines(x_coords_b, y_fitted_plot, col = COLORS$fit, lwd = 2)

  # # Mark detected outliers
  # outlier_indices_orig <- (1:N_mod)[mod$coefs$iis >= PIP_THRESHOLD]
  # if(length(outlier_indices_orig) > 0) {
  #   outlier_indices_plot <- outlier_indices_orig + sapply(outlier_indices_orig, function(x) {
  #     sum(country_boundaries_b < x)
  #   })
  #   abline(v = outlier_indices_plot, col = COLORS$outl_marker, lty = 2, lwd = 1.5)
  # }

  # Overlay observed points
  points(
    x_coords_b, y_plot,
    cex = 0.6,
    pch = 19,
    col = adjustcolor(COLORS$main, alpha.f = 0.4)
  )

  # Add legend
  legend(
    "topright",
    legend = c("Observed", "BISAM-fit"), # , "Detected Outlier"
    col = c(adjustcolor(COLORS$main, alpha.f = 0.9), COLORS$fit), #, COLORS$outl_marker
    pch = c(19, NA, NA),
    lty = c(NA, 1, 2),
    lwd = c(NA, 2, 2),
    cex = 1.25,
    bg = adjustcolor("white", alpha.f = 0.5)
  )

  # Add country labels
  midpoints_b <- sapply(1:n_mod, function(i) {
    start_idx <- if(i == 1) 1 else (i-1) * t_mod + 1
    end_idx <- i * t_mod
    mean(c(start_idx + (i-1), end_idx + (i-1)))
  })
  text(
    x = midpoints_b,
    y = min(y_plot, na.rm = TRUE),
    labels = cN,
    pos = 1,
    xpd = TRUE,
    cex = 1.25
  )

  # Add panel title
  mtext(
    side = 3, line = -0.3,
    text = bquote("Panel B: Fitted values of y"),
    adj = 0, font = 2, cex = 1.5
  )

  # Add overall title
  # mtext(
  #   text = bquote("Transport Emissions Breaks"),
  #   side = 3,
  #   line = 0.5,
  #   outer = TRUE,
  #   cex = 2,
  #   font = 2
  # )
}

dev.off()



# ==============================================================================
# END OF SCRIPT
# ==============================================================================