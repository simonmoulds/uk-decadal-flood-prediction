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
## config = read_yaml('config/config_1.yml')
## experiment = 'observed'
## aggregation_period = 'yr2to5' #to5_lag'
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
dir.create(file.path(output_dir, "prediction"))
dir.create(file.path(output_dir, "simulation"))
dir.create(file.path(output_dir, "fit"))

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

  ## ## Create output directories (one for prediction, one for skill)
  ## subdirs = c("prediction", "skill")
  ## for (m in 1:length(subdirs)) {
  ##   sd_output_dir = file.path(output_dir, subdirs[m])
  ##   unlink(sd_output_dir, recursive = TRUE)
  ##   dir.create(sd_output_dir, recursive = TRUE)
  ## }

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

    ## Normalize discharge if continuous distribution
    if (!experiment_conf$model_family == "PO") {
      catchment_area = metadata[["catchment-area"]][metadata$id %in% stn_id] # km2
      catchment_data =
        catchment_data %>%
        mutate(Q = Q * 24 * 60 * 60 / catchment_area / 1000 / 1000 * 1000) # m3/s -> mm/day
    }

    ## Handle missing data
    ## FIXME - this problem arises because we no longer download precipitation
    if ("P_sum" %in% names(catchment_data))
      catchment_data <- catchment_data %>% dplyr::select(-P_sum)

    exclude = catchment_data$missing_pct > 30 | !complete.cases(catchment_data)
    if (experiment_conf$model_family == "GA") {
      ## Cannot fit a Gamma distribution if the response variable contains zeroes
      exclude = exclude | (catchment_data[[experiment_conf$predictand]] == 0.)
    }

    if (!experiment_conf$model_family %in% c("GA", "PO")) {
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
          file.path(output_dir, "fit"), format = "parquet"
        )

      ## Now run prediction, using either forward chain or cv
      if (method == "forward") {
        lead_time <- config$aggregation_period[[aggregation_period]]$lead_time
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
      ## Write output
      catchment_prediction %>%
        mutate(
          ID = stn_id,
          date = NA,
          subset = subset
        ) %>%
        group_by(ID, model, subset) %>% #predictand, subset) %>%
        write_dataset(
          file.path(output_dir, "prediction"), format = "parquet"
        )

      if (method == "forward") {
        catchment_simulation %>%
          mutate(
            ID = stn_id,
            date = NA,
            subset = subset
          ) %>%
        group_by(ID, model, subset) %>% #predictand, subset) %>%
        write_dataset(
          file.path(output_dir, "simulation"), format = "parquet"
        )
      }
    }
    setTxtProgressBar(pb, k)
  }
  close(pb)
} else {
  warning(paste0("Aggregation period ", aggregation_period, " not specified for experiment ", experiment))
}

## ## TESTS
## ## Compare custom implementation of crpss with that in easyVerification
## ds <- open_dataset("results/exp1/analysis/hindcast/gamlss/yr2to5_lag/prediction/") %>% collect()
## ids <- ds$ID %>% unique()
## models <- ds$model %>% unique()
## subsets <- ds$subset %>% unique()
## years <- ds$year %>% unique() %>% sort()

## ## Forecast dimensions
## n_space <- length(ids)
## n_members <- 99
## n_time <- ds$year %>% unique() %>% length()

## ## for (i in 1:length(ids)) {
## for (i in 1:length(models)) {
##   for (j in 1:length(subsets)) {
##     fcst <- array(data=NA, dim=c(n_space, n_time, n_members))
##     obs <- array(data=NA, dim=c(n_space, n_time))
##     crpss <- rep(NA, length(ids))
##     for (k in 1:length(ids)) {
##       x <-
##         ds %>%
##         filter(model %in% models[i] & subset %in% subsets[j] & ID %in% ids[k]) %>%
##         dplyr::select(year, Q_95_obs, Q01:Q99) %>%
##         arrange(year)
##       x <- tibble(year = years) %>% left_join(x, by = "year")
##       fcst_k <- x %>% dplyr::select(Q01:Q99) %>% as.matrix() %>% na.omit()
##       obs_k <- x %>% dplyr::select(Q_95_obs) %>% as.matrix() %>% na.omit()
##       ## fcst[k,,] <- fcst_k
##       ## obs[k,] <- obs_k
##       crpss[k] <- veriApply("FairCrpss", fcst=fcst_k, obs=obs_k, strategy=list(type="crossval", blocklength=7), na.rm=TRUE)$skillscore
##     }
##     ## crpss <- veriApply("FairCrpss", fcst=fcst, obs=obs, strategy=list(type="crossval", blocklength=13), na.rm=T)
##   }
## }

## ds <- open_dataset("results/exp1/analysis/hindcast/gamlss/yr2to5_lag/prediction/") %>% collect()
## ids <- ds$ID %>% unique()
## models <- c("P", "P_T")
## subsets <- c("best_n", "full")

## comp <- list()
## for (i in 1:length(models)) {
##   for (j in 1:length(subsets)) {
##     for (k in 1:length(ids)) {
##       x <- ds %>% filter(model %in% models[i], subset %in% subsets[j], ID %in% ids[k])
##       ## plot(ds$Q_95_obs)
##       ## lines(ds$Q_95_exp)
##       ## lines(ds$Q25, lty="dashed")
##       ## lines(ds$Q85, lty="dashed")
##       fcst <- x %>% dplyr::select(Q01:Q99) %>% as.matrix() %>% na.omit()
##       obs <- x %>% dplyr::select(Q_95_obs) %>% as.matrix() %>% na.omit()
##       ev_crpss <- veriApply(
##         "FairCrpss",
##         fcst=fcst,
##         obs=obs,
##         strategy=list(type="crossval", blocklength=7),
##         na.rm=TRUE
##       )
##       my_crpss <- 1 - (mean(x$crps_ens_fcst) / mean(x$crps_climat))
##       comp[[length(comp) + 1]] <- data.frame(id=ids[k], subset=subsets[j], model=models[i], ev_crpss=ev_crpss$skillscore, my_crpss=my_crpss)
##     }
##   }
## }

## df <- do.call("rbind", comp)
## df <- df %>% filter(subset %in% "best_n") %>% group_by(id) %>% filter(my_crpss==max(my_crpss))

## p = ggplot(data=df, aes(x=ev_crpss, y=my_crpss)) +
##   geom_point(pch=21) +
##   geom_abline(intercept=0, slope=1) +
##   xlim(c(-0.25, 0.5)) +
##   ylim(c(-0.25, 0.5)) +
##   coord_fixed() +
##   theme_bw()

## ggsave(filename="~/crpss_comparison.png", width=4, height=5)
