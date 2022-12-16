#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(zoo)
library(cowplot)
library(yaml)
library(arrow)

options(bitmapType = "cairo")
options(dplyr.summarise.inform = FALSE)

## ## TESTING
## config = read_yaml('config/config_1.yml')
## obspath = 'results/intermediate/obs.parquet'
## fcstpath = 'results/intermediate/ensemble-forecast'
## aggr_period = 'yr2to9_lag'
## outputroot = 'results/exp2/analysis'
## cwd = 'workflow/scripts'

## extract configuration info
if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  config = read_yaml(args[1])
  obspath = args[2]
  fcstpath = args[3]
  aggr_period = args[4]
  outputroot = args[5]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "utils.R"))
config[["aggregation_period"]] = parse_config_aggregation_period(config)

## Parse aggregation period specification
period = config$aggregation_period[[aggr_period]]
lead_tm = period$lead_time
start = min(lead_tm)
end = max(lead_tm)
study_period = period$study_period

## Make output directory
outputdir = file.path(outputroot, aggr_period)
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)

## TODO put these in config somehow
vars <- c("nao", "ea", "amv", "european_precip", "uk_precip", "uk_temp")
antecedent_vars <- c("european_precip")
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

## ################################### ##
## Load observed data
## ################################### ##

dataset <- read_parquet(obspath)
obs <- get_obs_new(
  dataset, study_period, start = start, end = end,
  vars = vars, months = months
)
obs_antecedent <- get_obs_new(
  dataset, study_period, start = start, end = end,
  vars = antecedent_vars, months = antecedent_months
)
obs_antecedent <-
  obs_antecedent %>%
  mutate(variable = paste0(variable, "_antecedent")) %>%
  mutate(init_year = init_year + year_offset)

## Join together and save output
obs <- rbind(obs, obs_antecedent)
write_parquet(
  obs,
  file.path(outputdir, "observed.parquet")
)

## ################################### ##
## Load modelled data
## ################################### ##

lead_times <- lead_tm
dataset <- open_dataset(fcstpath)
ensemble_fcst <- get_hindcast_data_new(
  dataset, study_period, lead_times,
  vars = vars, months = months
)
ensemble_fcst_antecedent <- get_hindcast_data_new(
  dataset, study_period, lead_times,
  vars = antecedent_vars, months = antecedent_months
)
ensemble_fcst_antecedent <-
  ensemble_fcst_antecedent %>%
  mutate(variable = paste0(variable, "_antecedent")) %>%
  mutate(init_year = init_year + year_offset)

## Join together and save output
ensemble_fcst <- rbind(ensemble_fcst, ensemble_fcst_antecedent)
write_parquet(
  ensemble_fcst,
  file.path(outputdir, "ensemble_forecast.parquet")
)

