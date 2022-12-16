#!/usr/bin/env Rscript

library(tidyverse)
library(arrow)
library(yaml)

options(dplyr.summarise.inform = FALSE)

if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  stations_file = args[1]
  outputroot <- args[2]
  config <- read_yaml(args[3])
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "utils.R"))

stations <- read_csv(stations_file)
station_ids <- stations %>% filter(source %in% "GRDC") %>% `$`(id)

metadata <- read_csv(file.path(config$input_data_root, config$aux_data$grdc, "grdc_stations.csv"), show_col_types = FALSE)
iso3166_codes <- read_csv(file.path(config$input_data_root, config$aux_data$grdc, "data_csv.csv"), show_col_types = FALSE) %>%
  rename(country_name=Name, country=Code)
metadata <- metadata %>% left_join(iso3166_codes, by="country")
metadata <- metadata %>% filter(grdc_no %in% station_ids)
station_ids <- metadata$grdc_no %>% as.integer()

## data_files <- list.files(
##   file.path(datadir, "original_data_2020"),
##   pattern="^[0-9]+_Q_Day.Cmd.txt$",
##   full.names = TRUE
## )

pb = txtProgressBar(min=0, max=length(station_ids), initial=0)
for (i in 1:length(station_ids)) {
  ## fn <- data_files[i]
  ## id <- gsub("([0-9]+)_(Q_Day.Cmd.txt)", "\\1", basename(fn)) %>% as.integer()
  id <- station_ids[i]
  fn <- file.path(
    config$input_data_root,
    config$aux_data$grdc,
    "original_data_2020",
    paste0(id, "_Q_Day.Cmd.txt")
  )
  if (!file.exists(fn))
    next
  x <- read_table(fn, skip=36, show_col_types = FALSE)
  x <- x %>% setNames(c("time", "Q")) %>%
    mutate(Q = ifelse(Q == -999.0, NA, Q)) %>%
    mutate(time = gsub("([0-9]{4}-[0-9]{2}-[0-9]{2})(.*)", "\\1", time)) %>%
    mutate(time = as.Date(time, format="%Y-%m-%d")) %>%
    mutate(year = format(time, "%Y") %>% as.integer) %>%
    mutate(month = format(time, "%m") %>% as.integer) %>%
    mutate(ID = id, .after=time)
  x <- summarise_discharge_data(x, metadata=NULL)
  x %>%
    group_by(ID) %>%
    write_dataset(file.path(outputroot, "grdc-discharge-summaries"), format = "parquet")
  setTxtProgressBar(pb, i)
}
close(pb)

metadata %>%
  mutate(across(where(is.logical), as.integer)) %>%
  write_parquet(file.path(outputroot, "grdc-metadata.parquet"))
