#!/usr/bin/env Rscript

library(tidyverse)
library(arrow)
library(rnrfa)
library(yaml)

options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  stations_file <- snakemake@input[["stations"]]
  season <- snakemake@wildcards[["season"]]
  outputroot <- snakemake@params[["outputroot"]]
  snakemake@source("utils.R")
} else {
  stations_file <- "resources/station_list_with_grid_coords.parquet"
  season <- "DJFM"
  outputroot <- "results/input/discharge/NRFA"
  cwd = "workflow/decadal-prediction-scripts/R"
  source(file.path(cwd, "utils.R"))
}

stations <- read_parquet(stations_file)
station_ids <- stations %>% filter(source %in% "UKBN") %>% `$`(id)
n_stations <- length(station_ids)

pb = txtProgressBar(min=0, max=n_stations, initial=0)
for (i in 1:n_stations) {
  stn_id <- as.numeric(station_ids[i])
  df <- download_nrfa_data(stn_id)
  df <- summarise_discharge_data(df, season)
  if (any(is.na(df$ID))) stop()
  df %>%
    group_by(ID) %>%
    write_dataset(
      file.path(outputroot, season),
      format = "parquet",
      hive_style = FALSE
    )
  setTxtProgressBar(pb, i)
}
close(pb)

