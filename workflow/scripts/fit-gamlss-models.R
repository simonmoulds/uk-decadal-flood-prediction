#!/usr/bin/env Rscript

library(tidyverse)
library(scales)
library(sf)
library(gamlss)
library(arrow)
library(yaml)
library(rnrfa)
library(optparse)

options(dplyr.summarise.inform = FALSE)

## ## FOR TESTING:
## config = read_yaml('config/config.yml')
## experiment = 'observed'
## aggregation_period = 'yr2to9'
## method = 'cv'
## outputroot = 'results/exp1'
## cwd = 'workflow/scripts/'

if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  config = read_yaml(args[1])
  experiment = args[2]
  aggregation_period = args[3]
  method = args[4]
  outputroot = args[5]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "utils.R"))

## Only parse the sections we need here
config[["subset"]] <- parse_config_subset(config)
config[["aggregation_period"]] <- parse_config_aggregation_period(config)
config[["modelling"]] <- parse_config_modelling(config)

output_dir = file.path(outputroot, "analysis", experiment, "gamlss", aggregation_period)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

metadata = catalogue()
experiment_conf = config$modelling[[experiment]]

## This allows us to select specific aggregation periods to model
## while specifying all possible aggregation periods in the Snakefile
if (aggregation_period %in% experiment_conf$aggregation_periods) {

  ## Load input dataset
  input_dir = file.path(outputroot, "analysis", aggregation_period, "input")
  if (!dir.exists(input_dir)) {
    next
  }
  ds <- open_dataset(input_dir) %>% collect()

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

  ## Create output directories (one for prediction, one for skill)
  subdirs = c("prediction", "skill")
  for (m in 1:length(subdirs)) {
    sd_output_dir = file.path(output_dir, subdirs[m])
    unlink(sd_output_dir, recursive = TRUE)
    dir.create(sd_output_dir, recursive = TRUE)
  }

  ## Set predictand
  ds[["Q"]] = ds[[experiment_conf$predictand]]

  ## Loop through catchments
  print(sprintf("Fitting models with %s input data", aggregation_period))
  station_ids = ds$ID %>% unique() %>% sort()
  pb = txtProgressBar(min=0, max=length(station_ids), initial=0, title=pb_title)
  for (k in 1:length(station_ids)) {

    ## ############################### ##
    ## Prepare input data
    ## ############################### ##
    stn_id = station_ids[k]
    catchment_data =
      ds %>%
      filter(ID %in% stn_id & year %in% experiment_conf$study_period)

    ## Normalize discharge
    catchment_area = metadata[["catchment-area"]][metadata$id %in% stn_id] # km2
    catchment_data =
      catchment_data %>%
      mutate(Q = Q * 24 * 60 * 60 / catchment_area / 1000 / 1000 * 1000) # m3/s -> mm/day

    ## Handle missing data
    exclude = catchment_data$missing_pct > 30 | !complete.cases(catchment_data)
    if (experiment_conf$model_family == "GA") {
      ## Cannot fit a Gamma distribution if the response variable contains zeroes
      exclude = exclude | (catchment_data[[experiment_conf$predictand]] == 0.)
    } else {
      msg <- paste0("Model family ", experiment_conf$model_family, " currently not supported")
      stop(msg)
    }

    ## If more than 33% of data points are missing then do not model this catchment
    if (sum(exclude) > (length(exclude) * 0.33)) {
      next
    }
    catchment_data = catchment_data[!exclude,]

    ## Divide input data into subsets
    x = catchment_data
    catchment_data_list = list()
    for (m in 1:length(subsets)) {
      subset = subsets[m]
      catchment_data_list[[m]] = x %>% filter(subset %in% subsets[m])
    }

    ## ############################### ##
    ## Fit models
    ## ############################### ##
    for (m in 1:length(catchment_data_list)) {
      xx = catchment_data_list[[m]]
      subset = subsets[m]
      n_fold = nrow(xx)
      catchment_prediction_list = list()

      ## ## Fit models on all data to extract AIC
      ## models = fit_models(
      ##   experiment_conf$formulas,
      ##   experiment_conf$sigma_formulas,
      ##   experiment_conf$model_family,
      ##   xx
      ## )
      ## aic = get_aic(models)

      if (method == "forward") {
        lead_time <- config$aggregation_period[[aggregation_period]]$lead_time
        catchment_prediction <- fit_models_forward_chain(
          xx,
          experiment_conf,
          training_period_start = 1961,
          training_period_end = 1979,
          test_period_end = 2006,
          lead_time
        )

      } else if (method == "cv") {
        lead_time <- config$aggregation_period[[aggregation_period]]$lead_time
        catchment_prediction <- fit_models_cv(
          xx,
          experiment_conf,
          lead_time
        )

      } else {
        stop("`method` must be either 'forward' or 'cv'")
      }
      if (is.null(catchment_prediction)) {
        next
      }

      catchment_prediction <- catchment_prediction %>% mutate(ID = stn_id)
      ## ## Make complete time series
      ## complete_ts = expand_grid(
      ##   ID = stn_id,
      ##   ## clim_season = "DJFM",
      ##   year = experiment_conf$study_period,
      ##   model = unique(catchment_prediction$model)#,
      ##   ## predictand = experiment_conf$predictand
      ## )

      ## catchment_prediction =
      ##   complete_ts %>%
      ##   left_join(
      ##     catchment_prediction,
      ##     ## by=c("ID", "clim_season", "year", "model", "predictand")
      ##     by=c("ID", "year", "model", "predictand")
      ##   ) %>%
      ##   mutate(period = aggregation_period, subset = subsets[m], .after = predictand) %>%
      ##   ## mutate(lead_time = lead_tm, .before = model) %>%
      ##   filter(year %in% experiment_conf$study_period) %>%
      ##   arrange(year, model)

      ## Write output
      catchment_prediction %>%
        mutate(
          date = NA,
          model = paste0("GAMLSS_", model, "_", subset),
          subset = subset
        ) %>%
        group_by(ID, model, subset) %>% #predictand, subset) %>%
        write_dataset(output_dir, format = "parquet")
        ## write_dataset(file.path(output_dir, "prediction"), format = "parquet")

      ## model_nms = distinct(catchment_prediction, model)$model
      ## ## Evaluate model skill
      ## skill_scores_list = list()
      ## for (p in 1:length(model_nms)) {
      ##   model_nm = model_nms[p]
      ##   pred =
      ##     catchment_prediction %>%
      ##     filter(model %in% model_nm) %>%
      ##     na.omit()
      ##   skill =
      ##     mean_square_error_skill_score(pred$obs, pred$exp) %>%
      ##     as_tibble() %>%
      ##     mutate(
      ##       ID = stn_id,
      ##       model = model_nm,
      ##       ## lead_time = lead_tm,
      ##       predictand = experiment_conf$predictand,
      ##       period = aggregation_period,
      ##       subset = subsets[m],
      ##       .before = "msss"
      ##     )
      ##   idx = length(skill_scores_list) + 1
      ##   skill_scores_list[[idx]] = skill
      ## }
      ## skill_scores = do.call("rbind", skill_scores_list)
      ## aic = aic %>% pivot_longer(everything(), names_to = "model", values_to = "aic")
      ## skill_scores = left_join(skill_scores, aic, by = "model")
      ## skill_scores %>%
      ##   group_by(ID, predictand, subset) %>%
      ##   write_dataset(file.path(output_dir, "skill"), format = "parquet")
    }
    setTxtProgressBar(pb, k)
  }
  close(pb)

} else {
  warning(paste0("Aggregation period ", aggregation_period, " not specified for experiment ", experiment))
}
