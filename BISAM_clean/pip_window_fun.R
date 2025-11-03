pip_window <- function(mod, win_size, op=c(">=","=="), pip_threshold = 0.50, cred_int = 0.95) {
  
  data_pip <- mod$draws$omega
  data_steps <- mod$draws$sis
  # Extract column names
  col_names <- colnames(data_pip)
  
  # Parse unit and time from column names
  parsed <- strsplit(col_names, "\\.")
  units <- sapply(parsed, function(x) x[2])
  times <- sapply(parsed, function(x) x[3])
  
  times_all <- unique(times)[order(unique(times))] # used for gap check below
  
  # Get unique units
  unique_units <- unique(units)
  
  # Results list
  result <- list()
  
  # For each unit
  for (unit in unique_units) {
    # Get indices for this unit
    unit_indices <- which(units == unit)
    unit_times <- times[unit_indices]
    
    # Sort by time
    time_order <- order(as.numeric(unit_times))
    sorted_indices <- unit_indices[time_order]
    sorted_times <- unit_times[time_order]
    
    # Store results for all windows
    unit_results_pip <- list()
    unit_results_step <- list()
    unit_sign_purity <- list()
    
    # If the unit has fewer times than window length, use all available times
    if (length(sorted_times) < win_size) {
      # get all available date for that window
      window_pips <- data_pip[, sorted_indices, drop = FALSE]
      window_steps<- data_steps[, sorted_indices, drop = FALSE]

      # create window statistics per MCMC sample in all available times
      row_counts_pip <- rowSums(window_pips)
      row_sums_step <- rowSums(window_steps)
      
      # get average step purity
      step_purity <-  sum(sign(colMeans(window_steps)) == sign(mean(row_sums_step))) / 
        length(window_indices) * sign(mean(row_sums_step))
      
      # wrap up
      time_range <- paste(min(sorted_times), "-", max(sorted_times))
      unit_results_pip[[time_range]] <- row_counts_pip
      unit_results_step[[time_range]] <- row_sums_step
      unit_sign_purity[[time_range]] <- step_purity
      
      
    } else {
      # Sliding window approach
      for (i in 1:(length(sorted_times) - win_size + 1)) {
        window_indices <- sorted_indices[i:(i + win_size - 1)]
        window_times <- sorted_times[i:(i + win_size - 1)]
        
        pos_times <- which(times_all %in% window_times)
        
        # Check if times are consecutive (no gaps)
        if (all(diff(pos_times) == 1)) {
          # get data for this window
          window_pips <- data_pip[, window_indices, drop = FALSE]
          window_steps<- data_steps[, window_indices, drop = FALSE]
          
          # create window statistics per MCMC sample
          row_counts_pip <- rowSums(window_pips)
          row_sums_step <- rowSums(window_steps)
          
          # get average step purity
          step_purity <-  sum(sign(colMeans(window_steps)) == sign(mean(row_sums_step))) / 
            length(window_indices) * sign(mean(row_sums_step))

          # wrap up
          time_range <- paste(window_times[1], "-", window_times[length(window_times)])
          unit_results_pip[[time_range]] <- row_counts_pip
          unit_results_step[[time_range]] <- row_sums_step
          unit_sign_purity[[time_range]] <- step_purity
        }
      }
    }
    unit_results_pip <- do.call(cbind, unit_results_pip)
    unit_results_step <- do.call(cbind, unit_results_step)
    unit_sign_purity <- do.call(cbind, unit_sign_purity)
    
    if(op == ">=") {
      unit_results_pip <- colMeans(unit_results_pip >= 1)
    } else if(op == "==") {
      unit_results_pip <- colMeans(unit_results_pip == 1)
    } else {
      print("Only op=c('>=','==') is implemented")
    }
    
    unit_results_step_lower <- apply(unit_results_step, 2, quantile, p = (1 - cred_int) / 2) 
    unit_results_step_upper <- apply(unit_results_step, 2, quantile, p = 1 - (1 - cred_int) / 2) 
    unit_results_step <- colMeans(unit_results_step)

    if(any(unit_results_pip>=pip_threshold)) {
      result[[unit]] <- rbind(
        MAP = unit_results_step[unit_results_pip >= pip_threshold],
        CI_Lower = unit_results_step_lower[unit_results_pip >= pip_threshold],
        CI_Upper = unit_results_step_upper[unit_results_pip >= pip_threshold],
        PIP = unit_results_pip[unit_results_pip >= pip_threshold],
        purity = unit_sign_purity[unit_results_pip >= pip_threshold])
    }
  }
  
  results_table <- do.call(rbind, lapply(names(result), function(unit) {
    periods <- colnames(result[[unit]])
    PIP <- as.numeric(result[[unit]]["PIP",])
    MAP <- as.numeric(result[[unit]]["MAP",])
    CI_Lower <- as.numeric(result[[unit]]["CI_Lower",])
    CI_Upper <- as.numeric(result[[unit]]["CI_Upper",])
    purity <- as.numeric(result[[unit]]["purity",])
    start_times <- sapply(strsplit(periods, " - "), function(x) x[1])
    end_times <- sapply(strsplit(periods, " - "), function(x) x[2])
    
    data.frame(
      unit = unit,
      period = periods,
      start_time = start_times,
      end_time = end_times,
      MAP = MAP,
      CI_Lower = CI_Lower, 
      CI_Upper = CI_Upper,
      PIP = PIP,
      purity = purity,
      stringsAsFactors = FALSE
    )
  }))
  
  # Use the position in the FULL unit list, not just units with breaks
  results_table$unit_idx <- match(results_table$unit, unique_units)
  
  # Balanced panel starting in 2003
  base_time <- times_all[1]
  end_time <- times_all[length(times_all)]
  
  # Calculate period index (1-indexed)
  results_table$time_idx <- match(results_table$start_time, times_all) 
  # results_table$start_time - base_time + 1
  
  # Calculate number of possible breaks per unit
  t_admissible <- length(times_all)
  
  # Calculate overall plot index for balanced panel
  if(length(results_table$unit_idx > 0)) {
    results_table <- results_table[order(results_table$unit_idx, results_table$start_time), ]
  }
  results_table$position_idx <- (results_table$unit_idx - 1) * t_admissible + results_table$time_idx
  
  return(results_table)
}
