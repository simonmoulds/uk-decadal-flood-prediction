#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)
library(arrow)
library(sf)
library(rnrfa)
library(RcppRoll)
library(yaml)

options(dplyr.summarise.inform = FALSE)

## Extract configuration info
if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  print(args[1])
  config <- read_yaml(args[1])
  outputfile <- args[2]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "external/R/utils.R"))
config = parse_config_io(config)

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

## Next identify stations included in the CAMELS-GB dataset
camels <- st_read(
  file.path(
    config$aux_data$root,
    config$aux_data$camels$subdirectory,
    "data", "CAMELS_GB_catchment_boundaries.shp"
  )
)
camels_stations <- camels$ID_STRING %>% unique() %>% as.integer()

## ## Now identify stations included in the CAMELS-GB dataset
## ukbn2_stations <- read_csv(
##   file.path(config$aux_data$ukbn, "UKBN_Station_List_vUKBN2.0_1.csv"),
##   show_col_types = FALSE
## )
## ## Allow benchmark scores of 1 (caution) and 2 (suitable)
## ukbn2_stations <- ukbn2_stations[ukbn2_stations$High_Score >= 1, "Station"]
## ukbn2_stations <- unlist(ukbn2_stations) %>% unname()
## camels_stations <- camels_stations[camels_stations %in% ukbn2_stations]

## Now filter UKBN2 stations
metadata <- metadata %>% filter(id %in% camels_stations)
stations <- metadata$id %>% as.character()

## TODO changepoint analysis [see Lopez and Frances (2013)]

## Write output
conn <- file(outputfile)
writeLines(stations, conn)
close(conn)
