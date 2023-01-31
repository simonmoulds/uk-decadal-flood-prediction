#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(zoo)
library(cowplot)
library(yaml)
library(arrow)

options(bitmapType = "cairo")
options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  config <- snakemake@config
  obspath <- snakemake@input[["obs"]]
  fcstpath <- snakemake@input[["fcst"]]
  aggregation_period <- snakemake@wildcards[["aggr"]]
  season <- snakemake@wildcards[["season"]]
  outputroot <- snakemake@params[["outputroot"]]
  snakemake@source("utils.R")
} else {
  ## TESTING
  config <- read_yaml("config/config_2.yml")
  obspath <- "results/intermediate/observed.parquet"
  fcstpath <- "results/intermediate/ensemble-forecast"
  aggregation_period <- "yr2to5"
  season <- "JJAS"
  outputroot <- "results/meteo"
  cwd = "workflow/decadal-prediction-scripts/R"
  source(file.path(cwd, "utils.R"))
}

## Parse aggregation period specification
config[["aggregation_period"]] = parse_config_aggregation_period(config)
period = config$aggregation_period[[aggregation_period]]
lead_tm = period$lead_time
start = min(lead_tm)
end = max(lead_tm)
study_period = period$study_period

## Make output directory
output_dir = file.path(outputroot, aggregation_period, season)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## TODO put these in config somehow
vars <- c("nao", "ea", "amv", "european_precip", "uk_precip", "uk_temp")
months <- get_month_index(season)

## ################################### ##
## Load observed data
## ################################### ##

dataset <- read_parquet(obspath)
obs <- get_obs(
  dataset, study_period, start = start, end = end,
  vars = vars, months = months
)
write_parquet(
  obs,
  file.path(output_dir, "observed.parquet")
)

## ################################### ##
## Load modelled data
## ################################### ##

lead_times <- lead_tm
dataset <- open_dataset(
  fcstpath,
  partitioning = c("source_id", "member", "init_year", "variable")
)
ensemble_fcst <- get_hindcast_data(
  dataset, study_period, lead_times,
  vars = vars, months = months
)
write_parquet(
  ensemble_fcst,
  file.path(output_dir, "ensemble_forecast.parquet")
)

