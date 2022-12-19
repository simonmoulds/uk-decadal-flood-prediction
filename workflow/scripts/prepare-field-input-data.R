#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(zoo)
library(yaml)
library(arrow)

options(bitmapType = "cairo")
options(dplyr.summarise.inform = FALSE)

## TESTING
config = read_yaml('config/config_2.yml')
obspath = 'results/intermediate/observed-field'
fcstpath = 'results/intermediate/ensemble-forecast-field'
stations = 'results/exp2/stations.csv'
aggr_period = 'yr2to9_lag'
outputroot = 'results/exp1/analysis'
cwd = 'workflow/scripts'

## ## extract configuration info
## if (sys.nframe() == 0L) {
##   args <- commandArgs(trailingOnly=TRUE)
##   config <- read_yaml(args[1])
##   obspath <- args[2]
##   fcstpath <- args[3]
##   stations <- args[4]
##   aggr_period <- args[5]
##   outputroot <- args[6]
##   args <- commandArgs()
##   m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
##   cwd <- dirname(regmatches(args, m))
## }
source(file.path(cwd, "utils.R"))
config[["aggregation_period"]] = parse_config_aggregation_period(config)

## Parse aggregation period specification
period <- config$aggregation_period[[aggr_period]]
lead_tm <- period$lead_time
start <- min(lead_tm)
end <- max(lead_tm)
study_period <- period$study_period

## Make output directory
outputdir = file.path(outputroot, aggr_period)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)

## TODO put these in config somehow
## vars <- c("nao", "ea", "amv", "european_precip", "uk_precip", "uk_temp")
## antecedent_vars <- c("european_precip")
vars <- c("precip_field", "temp_field")
antecedent_vars <- c("precip_field")
months <- c(12, 1, 2, 3)
antecedent_months <- c(9, 10, 11)

## The antecedent year may need to be offset in some cases,
## e.g. for MAM the antecedent period could be DJF, which
## starts in the previous year
if (antecedent_months[1] >= months[1]) {
  year_offset <- 1
} else {
  year_offset <- 0
}

## Observed field data
metadata <- read_csv(stations, show_col_types = FALSE)

## FOR TESTING ONLY
metadata <- metadata %>% filter(source %in% "UKBN")

station_ids <- metadata$id
n_stations <- length(station_ids)

## ################################### ##
## Load observed data
## ################################### ##

dataset <- open_dataset(obspath)
pb = txtProgressBar(min=0, max=n_stations, initial=0)
for (i in 1:length(station_ids)) {
  id <- station_ids[i]
  subdataset <- dataset %>% filter(ID %in% id) %>% collect()
  obs_field <- get_obs_new(
    subdataset,
    study_period,
    start = start, end = end,
    vars = vars,
    months = months
  )
  obs_field_antecedent <- get_obs_new(
    subdataset,
    study_period,
    start = start,
    end = end,
    vars = antecedent_vars,
    months = antecedent_months
  )
  obs_field_antecedent <-
    obs_field_antecedent %>%
    mutate(variable = paste0(variable, "_antecedent")) %>%
    mutate(init_year = init_year + year_offset)

  ## Join together and save output as dataset
  obs_field <- rbind(obs_field, obs_field_antecedent)
  obs_field <-
    obs_field %>%
    mutate(ID = id, .before = init_year) %>%
    arrange(init_year, variable) %>%
    group_by(ID) %>%
    write_dataset(
      file.path(outputdir, "observed-field"),
      format = "parquet"
    )
  ## Update progress bar
  setTxtProgressBar(pb, i)
}
close(pb)

## ################################### ##
## Load modelled data
## ################################### ##

lead_times <- lead_tm
dataset <- open_dataset(fcstpath)
pb = txtProgressBar(min=0, max=n_stations, initial=0)
for (i in 1:length(station_ids)) {
  id <- station_ids[i]
  ## TODO test this function on ARC - currently taking far too long here
  subdataset <- dataset %>% filter(ID %in% id) # & source_id %in% "CanCM4")
  system.time(ensemble_fcst <- get_hindcast_data_new(
    subdataset, study_period, lead_times,
    vars = vars, months = months
  ))
  ensemble_fcst_antecedent <- get_hindcast_data_new(
    subdataset, study_period, lead_times,
    vars = antecedent_vars, months = antecedent_months
  )
  ensemble_fcst_antecedent <-
    ensemble_fcst_antecedent %>%
    mutate(variable = paste0(variable, "_antecedent")) %>%
    mutate(init_year = init_year + year_offset)

  ## Join together and save output
  ensemble_fcst <- rbind(ensemble_fcst, ensemble_fcst_antecedent)
  ensemble_fcst <-
    ensemble_fcst %>%
    mutate(ID = id) %>%
    group_by(ID) %>%
    write_dataset(
      file.path(outputdir, "ensemble-forecast-field"),
      format = "parquet"
    )
  ## Update progress bar
  setTxtProgressBar(pb, i)
}
close(pb)
