#!/usr/bin/env Rscript

library(tidyverse)
library(gamlss)
library(arrow)
library(yaml)

options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  config <- snakemake@config
  experiment <- snakemake@wildcards[["expm"]]
  aggregation_period <- snakemake@wildcards[["aggr"]]
  season <- snakemake@wildcards[["season"]]
  stn_id <- snakemake@wildcards[["stn"]]
  method <- snakemake@params[["method"]]
  outputroot <- snakemake@params[["outputroot"]]
  snakemake@source("utils.R")
} else {
  ## TESTING
  config <- read_yaml("config/config_2.yml")
  experiment <- "hindcast_Q95"
  aggregation_period <- "yr2to5_lag"
  season <- "DJFM"
  stn_id <- read_parquet("resources/station_list_with_grid_coords.parquet")$id[2]
  method <- "cv"
  outputroot <- "results/output"
  cwd = "workflow/decadal-prediction-scripts/R"
  source(file.path(cwd, "utils.R"))
}

## Only parse the sections we need here
config[["subset"]] <- parse_config_subset(config)
config[["aggregation_period"]] <- parse_config_aggregation_period(config)
config[["modelling"]] <- parse_config_modelling(config)

lead_time <- config$aggregation_period[[aggregation_period]]$lead_time
experiment_conf = config$modelling[[experiment]]

output_dir = file.path(
  outputroot, "output", experiment, "gamlss", aggregation_period, season
)
dir.create(file.path(output_dir, "prediction"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "fit"), recursive = TRUE, showWarnings = FALSE)

metadata <-
  read_parquet("resources/station_list_with_grid_coords.parquet") %>%
  filter(id %in% stn_id) %>%
  na_if(-999)

## Load input dataset
input_dir <- file.path(
  outputroot, "input", "combined",
  aggregation_period, season, stn_id
)
ds <- open_dataset(
  input_dir, partitioning = c("subset")
) %>%
  collect() %>%
  mutate(ID = stn_id)

## In `ds`, column `year` currently represents the year of initialization.
## Here we change `year` to represent the first month of the prediction
## window, consistent with that used by the the ML/AI routines.
ds <- ds %>% mutate(year = year + lead_time - 1)

## Arrange columns and remove those no longer needed
ds <- ds %>%
  arrange(ID, year, period) %>%
  dplyr::select(-lead_time)

## Identify subsets (i.e. full ensemble, best n [NAO-matched], cmip5, cmip6)
subsets = experiment_conf$subsets
ds_subsets = ds$subset %>% unique() %>% sort()
if (!all(subsets %in% ds_subsets)) {
  stop("Dataset does not contain all specified subsets")
}

## Set predictand
ds[["Q"]] = ds[[experiment_conf$predictand]]
catchment_data =
  ds %>%
  filter(ID %in% stn_id & year %in% experiment_conf$study_period)

## Normalize discharge if continuous distribution
if (!experiment_conf$model_family == "PO") {
  catchment_area <- metadata$catchment_area # km2
  catchment_data <-
    catchment_data %>%
    mutate(Q = Q * 24 * 60 * 60 / catchment_area / 1000 / 1000 * 1000) # m3/s -> mm/day
}

## Handle missing data
## FIXME - this problem arises because we no longer download precipitation
catchment_data <- catchment_data %>% dplyr::select(-any_of("P_sum"))

exclude = catchment_data$missing_pct > 30 | !complete.cases(catchment_data)
if (experiment_conf$model_family == "GA") {
  ## Cannot fit a Gamma distribution if the response variable contains zeroes
  exclude = exclude | (catchment_data[[experiment_conf$predictand]] == 0.)
}

if (!experiment_conf$model_family %in% c("GA", "PO")) {
  msg <- paste0(
    "Model family ", experiment_conf$model_family, " currently not supported"
  )
  stop(msg)
}

## Only proceed if fewer than 33% of data points are missing
if (sum(exclude) > (length(exclude) * 0.33)) {
  ## Create empty directories to satisfy snakemake
  dir.create(file.path(output_dir, "fit", stn_id), recursive = TRUE)
  dir.create(file.path(output_dir, "prediction", stn_id), recursive = TRUE)

} else {

  ## Prepare data to supply to model fitting routine
  catchment_data <- catchment_data[!exclude,]
  x <- catchment_data
  catchment_data_list <- list()
  for (m in 1:length(subsets)) {
    subset <- subsets[m]
    catchment_data_list[[m]] <- x %>% filter(subset %in% subsets[m])
  }

  ## ############################### ##
  ## Fit models
  ## ############################### ##
  for (m in 1:length(catchment_data_list)) {
    xx = catchment_data_list[[m]]
    subset = subsets[m]
    n_fold = nrow(xx)
    catchment_prediction_list = list()

    ## Fit models on all data to assess information about fit
    models = fit_models(
      experiment_conf$formulas,
      experiment_conf$sigma_formulas,
      experiment_conf$model_family,
      xx
    )
    aic = get_aic(models)
    aic <-
      aic %>%
      pivot_longer(everything(), names_to = "model", values_to = "aic")
    residual_checks <- get_residual_checks(models)

    ## Write fit summary
    residual_checks %>%
      left_join(aic, by = "model") %>%
      mutate(
        ID = stn_id,
        subset = subset
      ) %>%
      group_by(ID, model, subset) %>%
      write_dataset(
        file.path(output_dir, "fit"), format = "parquet", hive_style = FALSE
      )

    ## Now run prediction, using either forward chain or cv
    if (method == "forward") {
      out <- fit_models_forward_chain(
        xx,
        experiment_conf,
        training_period_start = 1961,
        training_period_end = 1979,
        test_period_end = 2015, #06, # FIXME
        lead_time
      )
      catchment_prediction <- out$prediction
      catchment_simulation <- out$simulation
    } else if (method == "cv") {
      catchment_prediction <- fit_models_cv(
        xx,
        experiment_conf,
        lead_time
      )
    } else {
      stop("`method` must be either 'forward' or 'cv'")
    }
    ## Write output
    catchment_prediction %>%
      mutate(
        ID = stn_id,
        date = NA,
        subset = subset
      ) %>%
      group_by(ID, model, subset) %>%
      write_dataset(
        file.path(output_dir, "prediction"), format = "parquet", hive_style = FALSE
      )
  }
}
