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
  inputdir <- args[1]
  outputfile <- args[2]
  config <- read_yaml(args[3])
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "utils.R"))
config = parse_config_io(config, inputdir)

## ################################### ##
## UKBN [via NRFA]
## ################################### ##

ukbn_filename = file.path(inputdir, config$aux_data$ukbn, "UKBN_Station_List_vUKBN2.0_1.csv")

## ## TESTING
## inputdir = '/Users/simonmoulds/projects/decadal-flood-prediction/data-raw'
## ukbn_filename = file.path(inputdir, 'UKBN', 'UKBN_Station_List_vUKBN2.0_1.csv')

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

## Next identify stations included in the UKBN dataset
ukbn_stations <- read_csv(ukbn_filename, show_col_types = FALSE)

## Allow benchmark scores of 1 (caution) and 2 (suitable)
ukbn_stations <- ukbn_stations[ukbn_stations$High_Score >= 1, "Station"]
ukbn_stations <- unlist(ukbn_stations) %>% unname()

## Now filter UKBN stations
metadata <- metadata %>% filter(id %in% ukbn_stations)
ukbn_stations <-
  metadata %>%
  dplyr::select(id, latitude, longitude) %>%
  mutate(source = "UKBN") %>%
  setNames(c("id", "lat", "lon", "source"))

## ################################### ##
## GRDC
## ################################### ##

grdc_filename <- file.path(inputdir, config$aux_data$grdc, "grdc_stations.csv")
grdc_basins_filename <- file.path(inputdir, config$aux_data$grdc, "original_data_2020", "grdc_basins_smoothed.shp")
iso3166_codes_filename <- file.path(inputdir, config$aux_data$grdc, "data_csv.csv")

## ## TESTING
## grdc_filename <- file.path(inputdir, "GRDC", "grdc_stations.csv")
## grdc_basins_filename <- file.path(inputdir, "GRDC", "original_data_2020", "grdc_basins_smoothed.shp")
## iso3166_codes_filename <- file.path(inputdir, "GRDC", "data_csv.csv")

metadata = read_csv(grdc_filename, show_col_types = FALSE)
basins <-
  st_read(grdc_basins_filename) %>%
  as_tibble() %>%
  dplyr::select(GRDC_NO, TYPE, QUALITY) %>%
  setNames(c("grdc_no", "type", "quality"))
country_codes <-
  read_csv(iso3166_codes_filename, show_col_types = FALSE) %>%
  setNames(c("country_name", "country"))

metadata <- basins %>% left_join(metadata, by="grdc_no") %>% left_join(country_codes, by="country")

## First filter by record length
metadata <- metadata %>% filter(t_start <= 1975 & t_yrs >= 40)

## Now filter by quality
metadata <- metadata %>% filter(quality %in% c("High", "Medium"))

## Remove UK stations, since these are already available via NRFA
metadata <- metadata %>% filter(!country %in% "GB") # GB=United Kingdom

## Prepare dataset
grdc_stations <-
  metadata %>%
  dplyr::select(grdc_no, lat, long) %>%
  mutate(source = "GRDC") %>%
  setNames(c("id", "lat", "lon", "source"))

## ################################### ##
## Combine
## ################################### ##

stations <- rbind(ukbn_stations, grdc_stations)
write_csv(stations, outputfile)

## NOT USED:
##
## metadata = catalogue()
## names(metadata) = names(metadata) %>% gsub("-", "_", .)

## ## First we filter by record length
## metadata =
##   metadata %>%
##   mutate(
##     gdf_start_date = as.Date(gdf_start_date),
##     gdf_end_date = as.Date(gdf_end_date)
##   ) %>%
##   mutate(
##     gdf_record_length = time_length(
##       gdf_end_date - gdf_start_date, "years"
##     )
##   ) %>%
##   filter(
##     gdf_record_length >= 40 & gdf_start_date <= as.Date("1975-01-01")
##   )

## ## Next identify stations included in the CAMELS-GB dataset
## camels <- st_read(
##   file.path(
##     inputdir,
##     ## config$input_data_root,
##     config$aux_data$camels$subdirectory,
##     "data", "CAMELS_GB_catchment_boundaries.shp"
##   )
## )
## camels_stations <- camels$ID_STRING %>% unique() %>% as.integer()

## ## Now identify stations included in the CAMELS-GB dataset
## ukbn2_stations <- read_csv(
##   file.path(config$aux_data$ukbn, "UKBN_Station_List_vUKBN2.0_1.csv"),
##   show_col_types = FALSE
## )
## ## Allow benchmark scores of 1 (caution) and 2 (suitable)
## ukbn2_stations <- ukbn2_stations[ukbn2_stations$High_Score >= 1, "Station"]
## ukbn2_stations <- unlist(ukbn2_stations) %>% unname()
## camels_stations <- camels_stations[camels_stations %in% ukbn2_stations]

## ## Now filter UKBN2 stations
## metadata <- metadata %>% filter(id %in% camels_stations)
## stations <- metadata$id %>% as.character()

## ## TODO changepoint analysis [see Lopez and Frances (2013)]

## ## Write output
## conn <- file(outputfile)
## writeLines(stations, conn)
## close(conn)
