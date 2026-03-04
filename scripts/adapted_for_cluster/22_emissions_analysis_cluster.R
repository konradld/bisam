# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                       Analysis file for CO2 Break Detection
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# clean environment
rm(list=ls())

# load libraries
library(dplyr)
library(stringr)
library(gets)
library(getspanel)
library(Matrix)
library(ggplot2)
library(mombf)


# Settings ----------------------------------------------------------------

config <- list(
  ssvs_settings = expand.grid(
    sis_prior = c("mom", "imom"),
    tau = 1/seq(1, 10, by = 1)^2 * 4, 
    check_outl = c(TRUE, FALSE),
    beta_prior = c("f"),
    incl_prior = c("bern"),
    stringsAsFactors = FALSE
  ),
  gets_lvl = c(0.05, 0.01), # just use the setting as done by Koch et al (2022)
  date = "2026-02-27"
)

dir_res <- sprintf("./output/emissions/%s/", 
                   config$date)
dir_data <- "./data/CO2DriversEU_dataset_CLEAN.csv"

# To save summarized results in
res <- list()

# Load and prepare data ---------------------------------------------------

data <- read.csv(dir_data)[-1]

# Group specification
EU15   <- c("Austria", "Belgium", "Germany", "Denmark", "Spain", "Finland",
            "France", "United Kingdom", "Ireland", "Italy", "Luxembourg", 
            "Netherlands", "Greece", "Portugal", "Sweden")

data_    <- data[, c('country','year','ltransport.emissions','lgdp','lgdp_sq','lpop')]
dat      <- filter(data_, country %in% EU15, year>=1995)

i_names <- unique(dat$country)
t_names <- unique(dat$year)

n <- length(i_names)
t <- length(t_names)


# GETS --------------------------------------------------------------------

res$gets <- list() # initiate new sub-list

res$ssvs <- list() # initiate new sub-list
ls_gets <- list.files(sprintf("./output/emissions/%s/", 
                              config$date), full.names = TRUE)
ls_gets <- ls_gets[which(grepl("gets", ls_gets))]

for(ff in ls_gets){
  
  p.value <- stringr::str_extract(ff, "(?<=gets_).*?(?=.RDS)")
  
  gets_i <- readRDS(ff)
  
  # fix names
  gets_breaks <- gets_i$isatpanel.result$mean.results %>%
    rownames() %>%
    str_replace_all(c(
      "time" = "tfe.",
      "id" = "ife.",
      "(fesis|sis)" = "sis.",
      "iis" = "iis."
    ))
  
  gets_coefs <- gets_i$isatpanel.result$mean.results
  rownames(gets_coefs) <- gets_breaks
  
  ind_breaks <- gets_coefs |> 
    filter(grepl("sis", rownames(gets_coefs)))
  
  gets_breaks_info <- tibble(
    id = stringr::str_split(rownames(ind_breaks), "\\.", simplify = T)[,2], 
    gets_lvl = p.value,
    period = stringr::str_split(rownames(ind_breaks), "\\.", simplify = T)[,3],
    period_start = period, 
    period_end = period, 
    period_length = 1,
    MAP = round(ind_breaks$coef, 4), 
    PIP = 1-round(ind_breaks$`p-value`, 2), 
    significance = ind_breaks$`p-value` < as.numeric(p.value), 
    purity = sign(MAP)
  )
  
  res$gets[[as.character(p.value)]] <- gets_breaks_info
}


# BISAM -------------------------------------------------------------------

source("./R/pip_window_fun.R")

res$ssvs <- list() # initiate new sub-list
res$iis <- list()
ls_ssvs <- list.files(sprintf("./output/emissions/%s/", 
                              config$date), full.names = TRUE)
ls_ssvs <- ls_ssvs[which(grepl("ssvs", ls_ssvs))]

for(ff in ls_ssvs) {
  
  tau <- stringr::str_extract(ff, "(?<=tau-).*?(?=_)")
  check_outl <- TRUE # stringr::str_extract(ff, "(?<=checkOutlier-).*?(?=_)")
  sis_prior <- stringr::str_extract(ff, "(?<=prior-).*?(?=_)")
  
  ssvs_temp <- readRDS(ff)
  
  ind_breaks <- ssvs_temp$coef_list |> 
    filter(!is.na(PIP), 
           PIP > 0.5)
  if(nrow(ind_breaks) > 0) {
    ind_breaks_info <- tibble(
      id = stringr::str_split(rownames(ind_breaks), "\\.", simplify = T)[,2], 
      period = stringr::str_split(rownames(ind_breaks), "\\.", simplify = T)[,3],
      period_start = period, 
      period_end = period, 
      period_length = 1,
      MAP = round(ind_breaks$MAP, 4), 
      PIP = round(ind_breaks$PIP, 2), 
      significance = ind_breaks$CI_Lower * ind_breaks$CI_Upper > 0, 
      purity = sign(MAP)
    )
  } else {
    ind_breaks_info <- tibble(
      id = NA, 
      period = NA,
      period_start = NA, 
      period_end = NA, 
      period_length = NA,
      MAP = NA,
      PIP = NA, 
      significance = NA,
      purity = NA
    )
  }
  window_breaks <- pip_window(mod = ssvs_temp, 
                              win_size = 2, 
                              op = ">=", 
                              pip_threshold = 0.5)
  
  if(length(window_breaks$unit_idx) > 0) {
    window_breaks_info <- tibble(
      id = window_breaks$unit,
      period  = window_breaks$period, 
      period_start  = window_breaks$start_time, 
      period_end  = window_breaks$end_time,
      period_length = 2,
      MAP = round(window_breaks$MAP, 4),
      PIP = round(window_breaks$PIP, 2), 
      significance = window_breaks$CI_Lower * window_breaks$CI_Upper > 0,
      purity = round(window_breaks$purity, 2))
  } else {
    window_breaks_info <- tibble(
      id = NA, 
      period = NA,
      period_start = NA, 
      period_end = NA, 
      period_length = NA,
      MAP = NA,
      PIP = NA, 
      significance = NA,
      purity = NA
    )
  }
  
  breaks_info <- rbind(ind_breaks_info, window_breaks_info) |> 
    mutate(tau = tau, checkOutlier = check_outl, prior = sis_prior, .after = "id")
  
  res$ssvs[[ff]] <- breaks_info
  
  if(check_outl == TRUE) {
    iis_temp <- colMeans(ssvs_temp$draws$iis)
    iis_df <- tibble(id = stringr::str_split(colnames(ssvs_temp$draws$iis), "\\.", simplify = T)[,2], 
                     tau = tau, checkOutlier = check_outl, prior = sis_prior,
                     period = stringr::str_split(colnames(ssvs_temp$draws$iis), "\\.", simplify = T)[,3],
                     PIP = iis_temp) 
    res$iis[[ff]] <- iis_df
  } 
  
}

df_gets <- do.call(rbind, res$gets)
df_ssvs <- do.call(rbind, res$ssvs)
df_iis <- do.call(rbind, res$iis)



df_ssvs |> filter(prior == "imom", period_length == 1, PIP > 0.5) |>
  mutate(tau = as.numeric(tau)) |>  
  group_by(tau, checkOutlier) # |>
  # count() |>
  # ggplot(aes(x = tau, y = n, group = checkOutlier, color = checkOutlier)) +
  # geom_line() +
  # theme_bw() + labs(y = "Number of breaks", x = "Tau") +
  # theme(legend.position = "bottom")


df_iis |> filter(PIP > 0.5)

