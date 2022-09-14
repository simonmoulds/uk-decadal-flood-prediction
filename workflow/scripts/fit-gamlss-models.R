#!/usr/bin/env Rscript

library(tidyverse)
library(scales)
library(sf)
library(gamlss)
library(arrow)
library(yaml)
library(optparse)

options(dplyr.summarise.inform = FALSE)

## ## FOR TESTING:
## config = read_yaml('config/config.yml')
## experiment = 'hindcast'
## method = 'forward'
## outputroot = 'results/exp2'
## cwd = 'workflow/scripts/'

if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  config = read_yaml(args[1])
  experiment = args[2]
  method = args[3]
  outputroot = args[4]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "external/R/utils.R"))
config = parse_config(config)

experiment_conf = config$modelling[[experiment]]
for (j in 1:length(experiment_conf$aggregation_periods)) {
  label = experiment_conf$aggregation_periods[j]
  input_dir = file.path(outputroot, "analysis", label, "input")
  if (!dir.exists(input_dir)) {
    next
  }
  output_dir = file.path(outputroot, "analysis", experiment, "gamlss", label)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ## Load input dataset
  ds =
    open_dataset(input_dir) %>%
    collect() %>%
    arrange(ID, year, period) %>%
    dplyr::select(-lead_time)

  ## Identify subsets
  subsets = experiment_conf$subsets
  ds_subsets = ds$subset %>% unique() %>% sort()
  if (!all(subsets %in% ds_subsets)) {
    stop("Dataset does not contain all specified subsets")
  }

  ## Create output directories
  subdirs = c("prediction", "skill") #, "catchment_plots")
  for (m in 1:length(subdirs)) {
    sd_output_dir = file.path(output_dir, subdirs[m])
    unlink(sd_output_dir, recursive = TRUE)
    dir.create(sd_output_dir, recursive = TRUE)
  }
  ## Set predictand
  ds[["Q"]] = ds[[experiment_conf$predictand]]
  station_ids = ds$ID %>% unique() %>% sort()
  ## Loop through catchments
  print(sprintf("Fitting models with %s input data", label))
  pb = txtProgressBar(min=0, max=length(station_ids), initial=0, title=pb_title)
  for (k in 1:length(station_ids)) {
    stn_id = station_ids[k]
    catchment_data =
      ds %>%
      filter(ID %in% stn_id & year %in% experiment_conf$study_period)
    ## Handle missing data in training data
    exclude = catchment_data$missing_pct > 30 | !complete.cases(catchment_data)
    if (experiment_conf$model_family == "GA") {
      ## Cannot fit a Gamma distribution if the response variable contains zeroes
      exclude = exclude | (catchment_data[[experiment_conf$predictand]] == 0.)
    }
    if (sum(exclude) > (length(exclude) * 0.33)) {
      next
    }
    catchment_data = catchment_data[!exclude,]

    x = catchment_data
    catchment_data_list = list()
    for (m in 1:length(subsets)) {
      subset = subsets[m]
      catchment_data_list[[m]] = x %>% filter(subset %in% subsets[m])
    }

    for (m in 1:length(catchment_data_list)) {
      xx = catchment_data_list[[m]]
      subset = subsets[m]
      n_fold = nrow(xx)
      catchment_prediction_list = list()

      ## Fit models on all data to extract AIC
      models = fit_models(
        experiment_conf$formulas,
        experiment_conf$sigma_formulas,
        experiment_conf$model_family,
        xx
      )
      aic = get_aic(models)

      ## training_period_start <- 1963
      ## training_period_end <- 1978
      ## max_training_period_end <- 2004 - 7
      if (method == "forward") {
        catchment_prediction <- fit_models_forward_chain(
          xx,
          experiment_conf,
          training_period_start = 1960, #1963
          training_period_end = 1979,
          test_period_end = 2005,
          window = 7
        )
      } else if (method == "cv") {
        lead_time <- config$aggregation_period[[label]]$lead_time
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

      ## Make complete time series
      complete_ts = expand_grid(
        ID = stn_id,
        ## clim_season = "DJFM",
        year = experiment_conf$study_period,
        model = unique(catchment_prediction$model),
        predictand = experiment_conf$predictand
      )

      catchment_prediction =
        complete_ts %>%
        left_join(
          catchment_prediction,
          ## by=c("ID", "clim_season", "year", "model", "predictand")
          by=c("ID", "year", "model", "predictand")
        ) %>%
        mutate(period = label, subset = subsets[m], .after = predictand) %>%
        ## mutate(lead_time = lead_tm, .before = model) %>%
        filter(year %in% experiment_conf$study_period) %>%
        arrange(year, model)

      ## Write output
      catchment_prediction %>%
        group_by(ID, predictand, subset) %>%
        write_dataset(file.path(output_dir, "prediction"), format = "parquet")

      model_nms = distinct(catchment_prediction, model)$model
      ## Evaluate model skill
      skill_scores_list = list()
      for (p in 1:length(model_nms)) {
        model_nm = model_nms[p]
        pred =
          catchment_prediction %>%
          filter(model %in% model_nm) %>%
          na.omit()
        skill =
          mean_square_error_skill_score(pred$obs, pred$exp) %>%
          as_tibble() %>%
          mutate(
            ID = stn_id,
            model = model_nm,
            ## lead_time = lead_tm,
            predictand = experiment_conf$predictand,
            period = label,
            subset = subsets[m],
            .before = "msss"
          )
        idx = length(skill_scores_list) + 1
        skill_scores_list[[idx]] = skill
      }
      skill_scores = do.call("rbind", skill_scores_list)
      aic = aic %>% pivot_longer(everything(), names_to = "model", values_to = "aic")
      skill_scores = left_join(skill_scores, aic, by = "model")
      skill_scores %>%
        group_by(ID, predictand, subset) %>%
        write_dataset(file.path(output_dir, "skill"), format = "parquet")
    }
    setTxtProgressBar(pb, k)
  }
  close(pb)
}
