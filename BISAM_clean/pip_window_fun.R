pip_window <- function(mod, win_size, op=c(">=","=="), pip_threshold = 0.50) {
  
  data_pip <- mod$draws$omega
  data_steps <- mod$draws$sis
  # Extract column names
  col_names <- colnames(data_pip)
  
  # Parse unit and time from column names
  parsed <- strsplit(col_names, "\\.")
  units <- sapply(parsed, function(x) x[2]) # rename to make more general?
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
    time_order <- order(unit_times)
    sorted_indices <- unit_indices[time_order]
    sorted_times <- unit_times[time_order]
    
    # Store results for all windows
    unit_results_pip <- list()
    unit_results_step <- list()
    unit_sign_purity <- list()
    
    # If the unit has fewer times than window length, use all available times
    if (length(sorted_times) < win_size) {
      # get all availabe date for that window
      window_pips <- data_pip[, sorted_indices, drop = FALSE]
      window_steps<- data_steps[, sorted_indices, drop = FALSE]
      window_sign <- sign(window_steps)
      # create window statistics per MCMC sample in all available times
      row_counts_pip <- rowSums(window_pips)
      row_means_step <- rowMeans(window_steps)
      # get relative step purity
      row_rel_sign_cnt <- rowSums(window_sign) / length(window_sign)
      # wrap up
      time_range <- paste(min(sorted_times), "-", max(sorted_times))
      unit_results_pip[[time_range]] <- row_counts_pip
      unit_results_step[[time_range]] <- row_means_step
      unit_sign_purity[[time_range]] <- row_rel_sign_cnt
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
          window_sign <- sign(window_steps)
          # create window statistics per MCMC sample
          row_counts_pip <- rowSums(window_pips)
          row_means_step <- rowMeans(window_steps)
          # get relative step purity
          row_rel_sign_cnt <- rowSums(window_sign) / ncol(window_sign)
          # wrap up
          time_range <- paste(window_times[1], "-", window_times[length(window_times)])
          unit_results_pip[[time_range]] <- row_counts_pip
          unit_results_step[[time_range]] <- row_means_step
          unit_sign_purity[[time_range]] <- row_rel_sign_cnt
        }
      }
    }
    unit_results_pip <- do.call(cbind, unit_results_pip)
    unit_results_step <- do.call(cbind, unit_results_step)
    unit_sign_purity <- do.call(cbind, unit_sign_purity)
    # unit_sign_purity[is.nan(unit_sign_purity)] <- 0
    
    if(op == ">=") {
      unit_results_step_cond_on_incl <- sapply(1:ncol(unit_results_step), \(i) {
        mean(unit_results_step[unit_results_pip[, i] >= 1, i])
      })
      unit_sign_purity_cond_on_incl <- sapply(1:ncol(unit_sign_purity), \(i) {
        mean(unit_sign_purity[unit_results_pip[, i] >= 1, i])
      })
      unit_results_pip <- colMeans(unit_results_pip >= 1)
    } else if(op == "==") {
      unit_results_step_cond_on_incl <- sapply(1:ncol(unit_results_step), \(i) {
        mean(unit_results_step[unit_results_pip[, i] == 1, i])
      })
      unit_sign_purity_cond_on_incl <- sapply(1:ncol(unit_sign_purity), \(i) {
        mean(unit_sign_purity[unit_results_pip[, i] == 1, i])
      })
      unit_results_pip <- colMeans(unit_results_pip == 1)
    } else {
      print("Only op=c('>=','==') is implemented")
    }
    unit_results_step <- colMeans(unit_results_step)
    unit_sign_purity <- colMeans(unit_sign_purity)
    
    
    names(unit_sign_purity_cond_on_incl) <- 
      names(unit_results_step_cond_on_incl) <- 
      names(unit_results_pip)
    
    if(any(unit_results_pip>=pip_threshold)) {
      result[[unit]] <- rbind(
        pip = unit_results_pip[unit_results_pip >= pip_threshold],
        step = unit_results_step[unit_results_pip >= pip_threshold],
        pure = unit_sign_purity[unit_results_pip >= pip_threshold],
        cond_step = unit_results_step_cond_on_incl[unit_results_pip >= pip_threshold],
        cond_pure = unit_sign_purity_cond_on_incl[unit_results_pip >= pip_threshold])
    }
  }
  
  results_table <- do.call(rbind, lapply(names(result), function(unit) {
    periods <- colnames(result[[unit]])
    joint_pip <- as.numeric(result[[unit]]["pip",])
    step <- as.numeric(result[[unit]]["step",])
    purity <- as.numeric(result[[unit]]["pure",])
    cond_step <- as.numeric(result[[unit]]["cond_step",])
    cond_purity <- as.numeric(result[[unit]]["cond_pure",])
    start_times <- sapply(strsplit(periods, " - "), function(x) x[1])
    end_times <- sapply(strsplit(periods, " - "), function(x) x[2])
    
    data.frame(
      unit = unit,
      period = periods,
      start_time = start_times,
      end_time = end_times,
      joint_pip = joint_pip,
      step = step,
      purity = purity,
      cond_step = cond_step,
      cond_purity = cond_purity,
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
  results_table <- results_table[order(results_table$unit_idx, results_table$start_time), ]
  results_table$position_idx <- (results_table$unit_idx - 1) * t_admissible + results_table$time_idx
  
  return(results_table)
}