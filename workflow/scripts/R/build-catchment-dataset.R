#!/usr/bin/env Rscript

library(tidyverse)
library(arrow)
library(sf)
library(rnrfa)
library(lubridate)
library(RcppRoll)
library(yaml)

options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  config <- snakemake@config
  stations_file <- snakemake@input[["stations"]]
  aggregation_period <- snakemake@wildcards[["aggr"]]
  season <- snakemake@wildcards[["season"]]
  stn_id <- snakemake@wildcards[["stn"]]
  outputroot <- snakemake@params[["outputroot"]]
  snakemake@source("utils.R")
} else {
  ## TESTING
  config <- read_yaml("config/config_2.yml")
  stations_file <- "resources/station_list_with_grid_coords.parquet"
  aggregation_period <- "yr2to5_lag"
  season <- "DJFM"
  stn_id <- read_parquet(stations_file)$id[1]
  outputroot <- "results/input"
  cwd = "workflow/decadal-prediction-scripts/R"
  source(file.path(cwd, "utils.R"))
}

## Parse aggregation period specification
config[["aggregation_period"]] <- parse_config_aggregation_period(config)
config[["subset"]] <- parse_config_subset(config)
period <- config$aggregation_period[[aggregation_period]]
lead_tm <- period$lead_time
start <- min(lead_tm)
end <- max(lead_tm)
study_period <- period$study_period
save_obs <- period$observed
save_fcst <- period$hindcast
error_label <- period$error_name
error_season <- period$error_season
if (is.na(error_season))
  error_season <- season

stations <-
  read_parquet(stations_file) %>%
  filter(id %in% stn_id)

## TODO put these in config
extended_study_period <- 1960:2015
climate_vars <- c(
  "nao", "ea", "amv", "european_precip", "uk_precip", "uk_temp",
  "precip_field", "precip_field_antecedent", "temp_field"
)
antecedent_season <- get_antecedent_season(season)
antecedent_vars <- c("european_precip", "precip_field")

get_observed_discharge <- function(id, source) {
  if (source == "GRDC") {
    filepath = file.path(outputroot, 'discharge', 'GRDC', season, id)
  } else if (source == "UKBN") {
    filepath = file.path(outputroot, 'discharge', 'NRFA', season, id)
  }
  ds <- open_dataset(filepath) %>% collect()
  ds
}

observed_discharge_data <- get_observed_discharge(
  stations$id,
  stations$source
)

## ################################### ##
## ################################### ##
##
## Load climate data
##
## ################################### ##
## ################################### ##

## Observed
obs <- read_parquet(
  file.path(outputroot, 'meteo', aggregation_period, season, "observed.parquet")
)

obs_antecedent <- read_parquet(
  file.path(
    outputroot, 'meteo', aggregation_period,
    antecedent_season, "observed.parquet"
  )
) %>%
  filter(variable %in% antecedent_vars) %>%
  mutate(variable = paste0(variable, "_antecedent"))

obs <-
  rbind(obs, obs_antecedent) %>%
  arrange(init_year, variable)

## Ensemble [N.B. lagged (n=676)]
ensemble_fcst <- read_parquet(
  file.path(
    outputroot, "nao_matching", aggregation_period,
    season, "matched_ensemble.parquet"
  )
) %>%
  mutate(lag = init_year - init_year_matched + 1) %>%
  dplyr::select(-any_of("error"))

ensemble_fcst_antecedent <- read_parquet(
  file.path(
    outputroot, "nao_matching", aggregation_period,
    antecedent_season, "matched_ensemble.parquet"
  )
) %>%
  mutate(lag = init_year - init_year_matched + 1) %>%
  filter(variable %in% antecedent_vars) %>%
  mutate(variable = paste0(variable, "_antecedent")) %>%
  dplyr::select(-any_of("error"))

ensemble_fcst <-
  rbind(ensemble_fcst, ensemble_fcst_antecedent) %>%
  arrange(init_year, source_id, member, init_year_matched, variable)

## ################################### ##
## ################################### ##
##
## Load local climate
##
## ################################### ##
## ################################### ##

## Observed
obs_local <- open_dataset(
  file.path(
    outputroot, "meteo", aggregation_period,
    season, "observed-field", stations$coord
  )
) %>% collect()

obs_local_antecedent <- open_dataset(
  file.path(
    outputroot, "meteo", aggregation_period,
    antecedent_season, "observed-field", stations$coord
  )
) %>%
  collect() %>%
  filter(variable %in% antecedent_vars) %>%
  mutate(variable = paste0(variable, "_antecedent"))

obs_local <-
  rbind(obs_local, obs_local_antecedent) %>%
  arrange(init_year, variable)

## Ensemble [N.B. unlagged (n=169)]
ensemble_fcst_local <- open_dataset(
  file.path(
    outputroot, "meteo", aggregation_period,
    season, "ensemble-forecast-field", stations$coord
  )
) %>% collect()

ensemble_fcst_local_antecedent <- open_dataset(
  file.path(
    outputroot, "meteo", aggregation_period,
    antecedent_season, "ensemble-forecast-field", stations$coord
  )
) %>%
  collect() %>%
  filter(variable %in% antecedent_vars) %>%
  mutate(variable = paste0(variable, "_antecedent"))

ensemble_fcst_local <-
  rbind(ensemble_fcst_local, ensemble_fcst_local_antecedent) %>%
  arrange(init_year, source_id, member, variable)

## Combine obs & obs_local, ensemble_fcst & ensemble_fcst_local
obs <-
  rbind(obs, obs_local) %>%
  arrange(init_year, variable) %>%
  rename(value = obs)

ensemble_fcst <-
  ensemble_fcst %>%
  pivot_wider(-all_of(c("std")), names_from="variable", values_from="value")

ensemble_fcst_local <-
  ensemble_fcst_local %>%
  pivot_wider(names_from="variable", values_from="value") %>%
  rename(init_year_matched = init_year)

## Joining in this way (i.e. against init_year_matched)
## effectively lags the ensemble-field data
ensemble_fcst <-
  ensemble_fcst %>%
  left_join(
    ensemble_fcst_local,
    by=c("project", "source_id", "mip", "member", "init_year_matched")
  ) %>%
  pivot_longer(-(project:lag), names_to="variable", values_to="value")

## Load ensemble error data (from NAO-matching)
## Note that we can use a different aggregation period
## for the error [`error_label`]
ensemble_fcst_error <- read_parquet(
  file.path(
    outputroot, "nao_matching", error_label,
    error_season, "matched_ensemble_error.parquet"
  )
) %>%
  mutate(across(contains("init_year"), as.integer)) %>%
  arrange(source_id, member, init_year, init_year_matched)

## ################################### ##
## ################################### ##
##
## Load discharge data
##
## ################################### ##
## ################################### ##

## NB season_year is the year of the first month
dis_season <-
  observed_discharge_data %>%
  filter(clim_season %in% season) %>%
  arrange(season_year)

## Compute summary statistics for DJFM
dis_season_aggregated <- rolling_fun(
  dis_season$season_year,
  dis_season,
  cols = c(
    "missing_pct",
    "Q_max", "Q_mean", "Q_05", "Q_50", "Q_90", "Q_95"
  ),
  funs = list(
    missing_pct = mean,
    Q_max = mean, Q_mean = mean, Q_05 = mean,
    Q_50 = mean, Q_90 = mean, Q_95 = mean
  ),
  start = start, end = end
)

## Join with complete timeseries and update missing_pct
complete_annual_ts <- tibble(init_year = study_period)
dis_season_aggregated <-
  complete_annual_ts %>%
  left_join(dis_season_aggregated, by=c("init_year")) %>%
  mutate(missing_pct = ifelse(is.na(missing_pct), 100, missing_pct))

## ################################### ##
## ################################### ##
##
## Save observed data for catchment
##
## ################################### ##
## ################################### ##

## Save observed data (observed discharge + observed climate indices)
standardize <- function(x) {
  return((x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
}

obs_aggregated <-
  obs %>%
  group_by(variable) %>%
  filter(init_year %in% study_period) %>%
  mutate(value=standardize(value)) %>%
  pivot_wider(init_year, names_from="variable", values_from="value")

## Join with observed discharge data
dis_season_obs <-
  dis_season_aggregated %>%
  left_join(obs_aggregated, by="init_year")

## Write dataset to file
dis_season_obs %>%
  rename(year = init_year) %>%
  mutate(
    lead_time = min(lead_tm),
    period = aggregation_period,
    subset = "observed",
    id = stn_id
  ) %>%
  arrange(year) %>%
  group_by(id, subset) %>%
  write_dataset(
    file.path(outputroot, "combined", aggregation_period, season),
    format = "parquet",
    hive_style = FALSE
  )

## ################################### ##
## ################################### ##
##
## Save hindcast data for catchment
##
## ################################### ##
## ################################### ##

## Hindcast dataset consists of observed discharge
## + hindcast climate indices
for (k in 1:length(config$subset)) {
  subset = config$subset[[k]]
  ensemble_fcst_subset <- create_ensemble_forecast(
    ensemble_fcst_error,
    ensemble_fcst,
    vars = climate_vars,
    model_select = subset$models,
    project_select = subset$projects,
    full = subset$full,
    best_n = subset$best_n,
    worst_n = subset$worst_n,
    n_select = 20
  )
  fcst <-
    dis_season_aggregated %>%
    left_join(ensemble_fcst_subset, by = "init_year")
  fcst %>%
    rename(year = init_year) %>%
    mutate(
      lead_time = min(lead_tm),
      period = aggregation_period,
      subset = subset$name,
      id = stn_id
    ) %>%
    arrange(year) %>%
    group_by(id, subset) %>%
    write_dataset(
      file.path(outputroot, "combined", aggregation_period, season),
      format = "parquet",
      hive_style = FALSE
    )
}
