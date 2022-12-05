#!/usr/bin/env Rscript

library(tidyverse)
library(arrow)

options(dplyr.summarise.inform = FALSE)

if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  datadir = args[1]
  outputroot <- args[2]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "utils.R"))

metadata <- read_csv(file.path(datadir, "grdc_stations.csv"), show_col_types = FALSE)
iso3166_codes <- read_csv(file.path(datadir, "data_csv.csv"), show_col_types = FALSE) %>%
  rename(country_name=Name, country=Code)
metadata <- metadata %>% left_join(iso3166_codes, by="country")

data_files <- list.files(
  file.path(datadir, "original_data_2020"),
  pattern="^[0-9]+_Q_Day.Cmd.txt$",
  full.names = TRUE
)

pb = txtProgressBar(min=0, max=length(grdc_ids), initial=0)
for (i in 1:length(data_files)) {
  fn <- data_files[i]
  id <- gsub("([0-9]+)_(Q_Day.Cmd.txt)", "\\1", basename(fn)) %>% as.integer()
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
