#!/usr/bin/env Rscript

library(tidyverse)
library(arrow)
library(rnrfa)
library(yaml)

options(dplyr.summarise.inform = FALSE)

if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  stations_file = args[1]
  outputroot <- args[2]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "utils.R"))

metadata = catalogue()

names(metadata) = names(metadata) %>% gsub("-", "_", .)

## conn <- file(stations_file)
## station_ids <- readLines(conn)
## station_ids <- as.integer(station_ids)
## close(conn)
stations <- read_csv(stations_file)
station_ids <- stations %>% filter(source %in% "UKBN") %>% `$`(id)

metadata <- metadata %>% filter(id %in% station_ids)
## metadata$id <- as.numeric(metadata$id)
station_ids <- metadata$id
n_stations = length(station_ids)

## ################################### ##
## 1 - Time series                     ##
## ################################### ##

pb = txtProgressBar(min=0, max=n_stations, initial=0)
for (i in 1:n_stations) {
  stn_id = as.numeric(station_ids[i])
  df = download_nrfa_data(stn_id, metadata)
  if (any(is.na(df$ID))) stop()
  df %>%
    group_by(ID) %>%
    write_dataset(file.path(outputroot, "nrfa-discharge-summaries"), format = "parquet")
  setTxtProgressBar(pb, i)
}
close(pb)

## ################################### ##
## 2 - Metadata                        ##
## ################################### ##

## Drop columns with complex types
df_cols <- which(sapply(metadata, FUN=function(x) inherits(x, "data.frame")))
metadata %>%
  dplyr::select(-df_cols) %>%
  mutate(across(where(is.logical), as.integer)) %>%
  write_parquet(file.path(outputroot, "nrfa-metadata.parquet"))
