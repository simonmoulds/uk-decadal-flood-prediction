#!/usr/bin/env Rscript

## Author : Simon Moulds
## Date   : Jan 2022

library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(scales)
library(arrow)
library(sf)
library(smoothr)
## library(gamlss)
## library(rnrfa)
library(yaml)

options(dplyr.summarise.inform = FALSE)
options(bitmapType = 'cairo') # For server

## ## FOR TESTING:
## config = read_yaml('config/config.yml')
## aggregation_period = "yr2to5_lag"
## outputroot = 'results/exp1'
## cwd = 'workflow/scripts'

config = read_yaml("config/config.yml")
aggregation_period = "yr2to5_lag"
outputroot = "results/exp2"
cwd = "workflow/scripts"

## ## Extract configuration info
## if (sys.nframe() == 0L) {
##   args = commandArgs(trailingOnly=TRUE)
##   config = read_yaml(args[1])
##   workflow = args[2]
##   aggregation_period = args[3]
##   outputroot = args[4]
##   args = commandArgs()
##   m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
##   cwd <- dirname(regmatches(args, m))
## }
source(file.path(cwd, "utils.R"))
source(file.path(cwd, "plotting.R"))
config[["modelling"]] <- parse_config_modelling(config)

skill_measure = "crpss"

output_dir = file.path(outputroot, "fig", aggregation_period)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

aggregation_period_label_list = list(
  "yr2" = "Year 2",
  "yr2to5_lag" = "Year 2-5",
  "yr6to9_lag" = "Year 6-9",
  "yr2to9_lag" = "Year 2-9"
)

obs_aggregation_period_list = list(
  "yr2" = "yr2",
  "yr2to5_lag" = "yr2to5",
  "yr6to9_lag" = "yr6to9",
  "yr2to9_lag" = "yr2to9"
)

aggregation_period_label = aggregation_period_label_list[[aggregation_period]]
obs_aggregation_period = obs_aggregation_period_list[[aggregation_period]]
obs_aggregation_period_label = aggregation_period_label

## ####################################################### ##
## ####################################################### ##
##
## Preamble
##
## ####################################################### ##
## ####################################################### ##

## Load model skill scores for observed and hindcast experiments
obs_skill_scores <- load_skill_scores(config, "observed", obs_aggregation_period)
skill_scores <- load_skill_scores(config, "hindcast", aggregation_period)
station_ids <- skill_scores$ID %>% unique()

## Load model fit metrics for hindcast experiment
fit <- load_model_fit(config, "hindcast", aggregation_period) %>% mutate(kurtosis = kurtosis - 3)

## For spatial plots:
uk_boundary =
  st_read("../data-raw/CNTR_RG_01M_2020_4326.shp") %>%
  filter(CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

europe_boundary =
  st_read("../data-raw/CNTR_RG_01M_2020_4326.shp") %>%
  filter(!CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

gauge_stns =
  catalogue() %>%
  rename(ID = id, area = "catchment-area") %>%
  filter(ID %in% station_ids) %>%
  dplyr::select(ID, name, area, latitude, longitude) %>%
  st_as_sf(coords=c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_transform(27700)

## Text size for plot labels
axis_title_size_large = 9
axis_title_size = 8
axis_title_size_small = 7
axis_label_size_large = 7
axis_label_size = 6
axis_label_size_small = 5
legend_label_size = 6
legend_title_size = 8
tag_label_size = 8
strip_label_size = 8

## Call the write plot script
if (workflow == "exp1") {
  source(file.path(cwd, "make-plots-exp1.R"))
} else if (workflow == "exp2") {
  source(file.path(cwd, "make-plots-exp2.R"))
}
