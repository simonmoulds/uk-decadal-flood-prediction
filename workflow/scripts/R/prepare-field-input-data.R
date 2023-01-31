#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(zoo)
library(yaml)
library(arrow)

options(bitmapType = "cairo")
options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  config <- snakemake@config
  obspath <- snakemake@input[["obs"]]
  fcstpath <- snakemake@input[["fcst"]]
  grid_coord <- snakemake@wildcards[["grid"]]
  aggregation_period <- snakemake@wildcards[["aggr"]]
  season <- snakemake@wildcards[["season"]]
  outputroot <- snakemake@params[["outputroot"]]
  snakemake@source("utils.R")
} else {
  ## TESTING
  config <- read_yaml("config/config_1.yml")
  obspath <- "results/intermediate/observed-field.parquet"
  fcstpath <- "results/intermediate/ensemble-forecast-field"
  grid_coord <- "s50e115"
  aggregation_period <- "yr2to5"
  season <- "JJAS"
  outputroot <- "results/input/meteo"
  cwd = "workflow/decadal-prediction-scripts/R"
  source(file.path(cwd, "utils.R"))
}

## Parse aggregation period specification
config[["aggregation_period"]] = parse_config_aggregation_period(config)
period <- config$aggregation_period[[aggregation_period]]
lead_tm <- period$lead_time
start <- min(lead_tm)
end <- max(lead_tm)
study_period <- period$study_period

## Make output directory
output_dir = file.path(outputroot, aggregation_period, season)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## TODO put these in config somehow
vars <- c("precip_field", "temp_field")
months <- get_month_index(season)

## ################################### ##
## Load observed data
## ################################### ##

dataset <- read_parquet(obspath)
subdataset <-
  dataset %>%
  filter(coord %in% grid_coord) %>%
  collect()

obs_field <- get_obs(
  subdataset,
  study_period,
  start = start, end = end,
  vars = vars,
  months = months
)

obs_field <-
  obs_field %>%
  mutate(coord = grid_coord, .before = init_year) %>%
  arrange(init_year, variable) %>%
  group_by(coord) %>%
  write_dataset(
    file.path(output_dir, "observed-field"),
    format = "parquet",
    hive_style = FALSE
  )

## ################################### ##
## Load modelled data
## ################################### ##

lead_times <- lead_tm
dataset <- open_dataset(
  file.path(fcstpath, grid_coord),
  partitioning = c("source_id", "member", "variable")
) %>% collect()

ensemble_fcst <- get_hindcast_data(
  dataset, study_period, lead_times,
  vars = vars, months = months
)

ensemble_fcst <-
  ensemble_fcst %>%
  mutate(coord = grid_coord) %>%
  group_by(coord) %>%
  write_dataset(
    file.path(output_dir, "ensemble-forecast-field"),
    format = "parquet",
    hive_style = FALSE
  )
