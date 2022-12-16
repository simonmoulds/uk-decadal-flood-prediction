#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)
library(arrow)
library(sf)
library(rnrfa)
library(RcppRoll)
library(yaml)

options(dplyr.summarise.inform = FALSE)

if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  inputdir <- args[1]
  outputfile <- args[2]
  config <- read_yaml(args[3])
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "utils.R"))
config = parse_config_io(config, inputdir)

metadata = catalogue()
names(metadata) = names(metadata) %>% gsub("-", "_", .)

## First we filter by record length
metadata =
  metadata %>%
  mutate(
    gdf_start_date = as.Date(gdf_start_date),
    gdf_end_date = as.Date(gdf_end_date)
  ) %>%
  mutate(
    gdf_record_length = time_length(
      gdf_end_date - gdf_start_date, "years"
    )
  ) %>%
  filter(
    gdf_record_length >= 40 & gdf_start_date <= as.Date("1975-01-01")
  )

## Next identify stations included in the UKBN2 dataset
ukbn_stations <- read_csv(
  file.path(
    inputdir,
    config$aux_data$ukbn, "UKBN_Station_List_vUKBN2.0_1.csv"
  ),
  show_col_types = FALSE
)
## Allow benchmark scores of 1 (caution) and 2 (suitable)
ukbn_stations <- ukbn_stations[ukbn_stations$High_Score >= 1, "Station"]
ukbn_stations <- unlist(ukbn_stations) %>% unname()

## Now filter UKBN2 stations
metadata = metadata %>% filter(id %in% ukbn_stations)
stations <- metadata$id %>% as.character()

## Now filter UKBN stations
metadata <- metadata %>% filter(id %in% ukbn_stations)
ukbn_stations <-
  metadata %>%
  dplyr::select(id, latitude, longitude) %>%
  mutate(source = "UKBN") %>%
  setNames(c("id", "lat", "lon", "source"))

write_csv(ukbn_stations, outputfile)

## ## Write output
## conn <- file(outputfile)
## writeLines(stations, conn)
## close(conn)
