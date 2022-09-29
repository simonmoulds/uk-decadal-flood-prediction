#!/usr/bin/env Rscript

## ####################################################### ##
##
## Author : Simon Moulds
## Date   : Nov 2021 - July 2022
## Purpose:
## To implement the NAO-matching technique outlined by
## Smith et al. [S20]
##
## ####################################################### ##

library(tidyverse)
library(magrittr)
library(zoo)
library(cowplot)
library(yaml)
library(arrow)

options(bitmapType = "cairo")
options(dplyr.summarise.inform = FALSE)

## ## TESTING
## config = read_yaml('config/config.yml')
## obspath = 'results/intermediate/obs.parquet'
## fcstpath = 'results/intermediate/ensemble-forecast'
## aggr_period = 'yr2to9_lag'
## outputroot = 'results/exp2/analysis'
## cwd = 'workflow/scripts/external/R'

## Extract configuration info
if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  config = read_yaml(args[1])
  obspath = args[2]
  fcstpath = args[3]
  aggr_period = args[4]
  outputroot = args[5]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
## TODO put these functions in an R package
source(file.path(cwd, "utils.R"))
## config = parse_config_io(config)
config[["aggregation_period"]] = parse_config_aggregation_period(config)

## Variable against which to perform the mode-matching approach
## TODO put this in configuration
match_var <- "nao"
models <- c(
  config$ensemble_data$cmip5$models,
  config$ensemble_data$cmip6$models
)
climate_vars = c("nao", "ea", "amv", "european_precip", "uk_precip", "uk_temp")

## Parse aggregation period specification
## period = config$aggregation_period[[j]]
period = config$aggregation_period[[aggr_period]]
lead_tm = period$lead_time
start = min(lead_tm)
end = max(lead_tm)
study_period = period$study_period
n_years = length(study_period)
lag = period$lag
n_lag = period$n_lag
label = period$name
if (!period$hindcast)
  next

## Make output directories
outputdir = file.path(outputroot, aggr_period)
dir.create(outputdir, recursive = TRUE)

## Load observed data, save for later use
obs = get_obs(obspath, study_period, start = start, end = end)
obs = obs %>%
  pivot_longer(all_of(climate_vars), names_to = "variable", values_to = "obs")
write_parquet(
  obs,
  file.path(outputdir, "obs_study_period.parquet")
)
lead_times = lead_tm
ensemble_fcst_raw_complete = get_hindcast_data(
  fcstpath,
  ## "data",
  ## config$output_data$root,
  study_period,
  lead_times
)

## ################################### ##
## Prepare ensemble forecast           ##
## ################################### ##

## First of all we need to aggregate values if a multi-year
## period is used.
vars = c("nao", "ea", "amv", "european_precip", "uk_precip", "uk_temp")
group_vars = c("project", "mip", "source_id", "member", "init_year")
anomaly_group_vars = c("source_id", "member")
ens_mean_group_vars = c("init_year", "variable")
ensemble_fcst_raw_complete =
  ensemble_fcst_raw_complete %>%
  group_by_at(group_vars) %>%
  summarize(across(all_of(vars), mean))

## Compute anomalies of each variable
compute_anomaly = function(x) x - mean(x, na.rm = TRUE)
ensemble_fcst =
  ensemble_fcst_raw_complete %>%
  group_by_at(anomaly_group_vars) %>%
  mutate(across(all_of(vars), compute_anomaly)) %>%
  ungroup()
write_parquet(
  ensemble_fcst,
  file.path(outputdir, "ensemble_fcst.parquet")
)

## Take the ensemble mean
fcst =
  ensemble_fcst %>%
  pivot_longer(all_of(vars), names_to = "variable", values_to = "value") %>%
  group_by_at(ens_mean_group_vars) %>%
  summarize(
    ens_mean = mean(value, na.rm=TRUE),
    ens_q95 = quantile(value, probs = 0.95, na.rm = TRUE, names = FALSE),
    ens_q05 = quantile(value, probs = 0.05, na.rm = TRUE, names = FALSE)
  ) %>%
  ungroup() %>%
  left_join(obs, by = c("init_year", "variable"))

## ################################### ##
## Variance adjustment                 ##
## ################################### ##

## Here we lag the ensemble mean forecast by taking the mean
## over the current year and the previous three years (i.e.
## 4 years in total, as described in Smith et al. in Methods
## ("Lagged ensemble")).
fcst =
  fcst %>%
  group_by(variable) %>%
  mutate(
    ens_mean_lag = zoo::rollmean(
      ens_mean,
      n_lag,
      na.pad=TRUE,
      align="right"
    )
  )

## Standardise (i.e. divide by standard deviation)
fcst =
  fcst %>%
  mutate(
    obs_std = obs / sd(obs),
    ens_mean_lag_std = ens_mean_lag / sd(ens_mean_lag, na.rm=TRUE),
    ens_mean_std = ens_mean / sd(ens_mean_lag, na.rm=TRUE)
  )

## Multiply standardised values by std dev
## of observed data so that the data ranges
## are comparable.
fcst =
  fcst %>%
  mutate(
    ens_mean_lag_var_adj = ens_mean_lag_std * sd(obs),
    ens_mean_var_adj = ens_mean_std * sd(obs)
  )

## Multiply by ACC (Note that standardising
## and multiplying by ACC is equivalent to
## multiplying by RPS)
fcst =
  fcst %>%
  mutate(
    corr = corr_cross_validate(ens_mean_lag_std, obs_std, 1)
  ) %>%
  mutate(
    ens_mean_lag_var_adj = corr * ens_mean_lag_var_adj,
    ens_mean_var_adj = corr * ens_mean_var_adj,
    ens_mean_lag_std = corr * ens_mean_lag_std,
    ens_mean_std = corr * ens_mean_std
  )

write_parquet(
  fcst,
  file.path(outputdir, "ensemble_mean_fcst.parquet")
)

## ################################### ##
## NAO matching                        ##
## ################################### ##

## The ensemble mean forecast data was standardized
## in the previous section. Here we standardize the
## individual ensemble members.
group_vars = c("source_id", "member", "variable")
ensemble_fcst =
  ensemble_fcst %>%
  pivot_longer(all_of(vars), names_to = "variable", values_to = "value") %>%
  group_by_at(group_vars) %>%
  mutate(std = value / sd(value, na.rm=TRUE))

## Calculate error
nao_matched_ensemble_fcst_error = calculate_error(
  fcst,
  ensemble_fcst,
  n_years,
  n_lag,
  "nao"
)
## Join with ensemble data to get values
nao_matched_ensemble_fcst =
  nao_matched_ensemble_fcst_error %>%
  left_join(
    ensemble_fcst %>% rename(init_year_matched = init_year),
    by = c("project", "source_id", "mip", "member", "init_year_matched")
  )
write_parquet(
  nao_matched_ensemble_fcst,
  file.path(outputdir, "matched_ensemble.parquet")
)
nao_matched_ensemble_fcst_error =
  nao_matched_ensemble_fcst %>%
  group_by(project, source_id, mip, member, init_year_matched, init_year) %>%
  summarize(any_na = any(is.na(value)), error = mean(error))
write_parquet(
  nao_matched_ensemble_fcst_error,
  file.path(outputdir, "matched_ensemble_error.parquet")
)
## }
