################################################################################
# BISAM vs GETS F1 Comparison Analysis
# 
# Description: Comprehensive comparison of BISAM and GETS with different 
#              additional settings 
# 
# Output: Performance metrics, plots, and summary statistics
################################################################################

# Clear workspace
rm(list = ls())

#===============================================================================
#           Date of which the data is coming from
DATE <- "2026-07-24_variations"
#===============================================================================

# ==============================================================================
# 1. SETUP AND DATA LOADING
# ==============================================================================

# Load required libraries (only for data manipulation)
library(dplyr)

date <- DATE
gets_lvl <- "0.01"
bisam_prior <- "imom"
tau <- "3.31744830051061"

# Set up paths
data_path <- sprintf("./output/simulation/%s/", 
                     date)

# /gets_bisam_comparison_gets-%s_bisam_prior-%s_tau-%s, 
# gets_lvl,
# bisam_prior,
# tau

# Get list of all RDS files in the folder
file_list <- list.files(
  path = data_path, 
  pattern = "\\.RDS$", 
  full.names = TRUE
)


# file_list_renamed <- gsub("breaksize", "bs", file_list)
# file_list_renamed <- gsub("_breaknumber", "", file_list_renamed)
# file.rename(file_list[1], file_list_renamed[1])



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
  breaknumber <- sub(".*breaknumber-([a-z]+)_Nt.*", "\\1", filename)
  Nt <- as.numeric(sub(".*Nt-([0-9.]+)_Ni.*", "\\1", filename))
  Ni <- as.numeric(sub(".*Ni-([0-9.]+)_Nx.*", "\\1", filename))
  Nx <- as.numeric(sub(".*Nx-([0-9.]+)_ife.*", "\\1", filename))
  ife <- sub(".*ife-([A-Z]+)_tfe.*", "\\1", filename)
  tfe <- sub(".*tfe-([A-Z]+)_outl.*", "\\1", filename)
  outl <- sub(".*outl-([A-Z]+)_sprior.*", "\\1", filename)
  
  do_sv <- sub(".*do-sv-([A-Z]+)_sprior.*", "\\1", filename)
  if(do_sv != filename) {
    outl <- sub(".*outl-([A-Z]+)_do-sv.*", "\\1", filename)
  } else {
    do_sv <- "FALSE"
  }
  
  sprior <- sub(".*sprior-([a-z]+)_sigma.*", "\\1", filename)
  if(sprior == filename) sprior <- "beta_bern"
  sigma <- sub(".*sigma-([0-9.]+)_rep.*", "\\1", filename)
  if(is.na(as.numeric(sigma))) {
    sigma <- "auto"
  }
  replicate <- as.integer(sub(".*_rep([0-9]+)\\.RDS", "\\1", filename))
  
  
  # Convert to data frame with identifiers
  df <- as.data.frame(mat)
  # df$row_name <- rownames(mat)
  df$breaksize <- breaksize
  df$breaknumber <- breaknumber
  df$Nt <- Nt
  df$Ni <- Ni
  df$Nx <- Nx
  df$ife <- ife
  df$tfe <- tfe
  df$outl <- outl
  df$do_sv <- do_sv
  df$sprior <- sprior
  df$sigma <- sigma
  df$replicate <- replicate
  
  return(df)
}) |> 
  bind_rows()


summarized_data <- combined_data |> 
  group_by(breaksize, breaknumber, Nt, Ni, Nx, ife, tfe, outl, do_sv, sprior, sigma) |>
  summarise(across(true:ssvs_t3_tr, function(x) sum(x))) |> 
  mutate(specification = case_when(
    Ni == 25 ~ "Ni = 25", 
    Ni == 50 ~ "Ni = 50", 
    Nx == 1 ~ "Nx = 1", 
    Nx == 2 ~ "Nx = 2",
    ife == "TRUE" & tfe == "TRUE" ~ "TWFE",
    ife == "TRUE" ~ "IFE",
    tfe == "TRUE" ~ "TFE",
    outl == "TRUE" & do_sv == "TRUE" ~ "Outlier + SV", # probably remove
    outl == "TRUE" ~ "Outlier",
    do_sv == "TRUE" ~ "Stoch. Vol.", # missing - need to decide between internal and stochvol
    sprior == "beta_bern" ~ "Beta-Bernoulli",
    sigma == "0.01" ~ "m0,n0 = 0.01",
    .default = "Baseline"
  ), .after = breaknumber) |> 
  ungroup() |> 
  transmute(breaksize, breaknumber, specification, 
            tr.ssvs, tr.gets, fp.ssvs, fp.gets, fn.ssvs, fn.gets, ssvs, gets, true) |> 
  filter(specification != "Outlier + SV")


metrics <- summarized_data |> 
  transmute(breaksize = paste0(breaksize, "SD"), 
            breaknumber = paste0(toupper(substr(breaknumber, 1, 1)), substr(breaknumber, 2, nchar(breaknumber))), 
            breaknumber = factor(breaknumber, levels = c("Sparse", "Dense")),
            specification = factor(specification, 
                                   levels = c("Baseline", "IFE", "TFE", "TWFE", "Nx = 1", "Nx = 2", 
                                              "Beta-Bernoulli", "Outlier", "Stoch. Vol.", "m0,n0 = 0.01",
                                              "Ni = 25", "Ni = 50")) ,
            prec.ssvs = tr.ssvs / (tr.ssvs + fp.ssvs), 
            rec.ssvs = tr.ssvs / (tr.ssvs + fn.ssvs), 
            f1.ssvs = 2 * (prec.ssvs * rec.ssvs) / (prec.ssvs + rec.ssvs),
            prec.gets = tr.gets / (tr.gets + fp.gets), 
            rec.gets = tr.gets / (tr.gets + fn.gets), 
            f1.gets = 2 * (prec.gets * rec.gets) / (prec.gets + rec.gets), 
            f1.relative = f1.ssvs / f1.gets, 
            f1.absolute = f1.ssvs - f1.gets)


# Create heatmap
p_rel <- ggplot(metrics, aes(x = as.factor(breaksize), y = forcats::fct_rev(specification), fill = f1.relative)) +
  geom_tile(color = "white", size = 0.8) +
  facet_wrap(breaknumber~.) +
  scale_fill_gradient2(
    low = "#D55E00",       # Orange for GETS
    mid = "#FFFFFF",       # White for equal
    high = "#0072B2",      # Blue for BISAM
    midpoint = 1
  ) +
  geom_text(aes(label = sprintf("%.2f", f1.relative)), 
            color = "black", 
            size = 3.5,
            fontface = "bold") +
  labs(
    # title = "Relative F1 Scores (BISAM vs GETS)",
    # subtitle = "Orange = Low Expression | Blue = High Expression",
    x = "Break Size",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 1), 
    strip.text = element_text(size = 16, face = "bold")
  ) # + coord_equal()  # Make tiles square



# Create heatmap
p_abs <- ggplot(metrics, aes(x = as.factor(breaksize), y = forcats::fct_rev(specification), fill = f1.absolute)) +
  geom_tile(color = "white", size = 0.8) +
  facet_wrap(breaknumber~.) +
  scale_fill_gradient2(
    low = "#D55E00",       # Orange for GETS
    mid = "#FFFFFF",       # White for equal
    high = "#0072B2",      # Blue for BISAM
    midpoint = 0
  ) +
  geom_text(aes(label = sprintf("%.2f", f1.absolute)), 
            color = "black", 
            size = 3.5,
            fontface = "bold") +
  labs(
    # title = "Relative F1 Scores (BISAM vs GETS)",
    # subtitle = "Orange = Low Expression | Blue = High Expression",
    x = "Break Size",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 1), 
    strip.text = element_text(size = 16, face = "bold")
  ) # + coord_equal()  # Make tiles square


pdf("./output/simulation/figure_heatmap.pdf", width = 12, height = 9, onefile = TRUE)
p_rel
p_abs
dev.off()

