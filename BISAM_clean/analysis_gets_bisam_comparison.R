# Get list of all RDS files in the folder
file_list <- list.files(path = "./Simulations/gets_bisam_comparison_gets-0.01_bisam_prior-imom/", 
                        pattern = "\\.RDS$", 
                        full.names = TRUE)

library(dplyr)
combined_data <- lapply(file_list, function(file) {
  # Read the matrix
  mat <- readRDS(file)
  
  # Extract identifiers
  filename <- basename(file)
  breaksize <- sub(".*breaksize-([0-9.]+SD)_rep.*", "\\1", filename)
  replicate <- as.integer(sub(".*_rep([0-9]+)\\.RDS", "\\1", filename))
  
  # Convert to data frame and add identifiers
  df <- as.data.frame(mat)
  df$row_name <- rownames(mat)
  df$breaksize <- breaksize
  df$replicate <- replicate
  
  return(df)
}) |> 
  bind_rows()

head(combined_data)



# Filter for breaksize = 1.0SD
data_1SD <- combined_data |> 
  filter(breaksize == "10.0SD")

# Get summary statistics for 1.0SD
summary_1SD <- data_1SD |> 
  select(-row_name, -breaksize, -replicate) |> colSums()
summary_1SD
