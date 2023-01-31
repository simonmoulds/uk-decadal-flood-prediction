#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(zoo)
library(cowplot)
library(yaml)
library(arrow)

options(bitmapType = "cairo")
options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  config <- snakemake@config
  obspath <- snakemake@input[["obs"]]
  fcstpath <- snakemake@input[["fcst"]]
  aggregation_period <- snakemake@wildcards[["aggr"]]
  season <- snakemake@wildcards[["season"]]
  outputroot <- snakemake@params[["outputroot"]]
  match_var <- snakemake@params[["match_var"]]
  snakemake@source("utils.R")
} else {
  config <- read_yaml('config/config_2.yml')
  obspath <- 'results/analysis/input/yr2to9_lag/DJFM/observed.parquet'
  fcstpath <- 'results/analysis/input/yr2to9_lag/DJFM/ensemble_forecast.parquet'
  aggregation_period <- 'yr2to9_lag'
  season <- "DJFM"
  outputroot <- "results/input/nao_matching"
  match_var <- "nao"
  cwd <- "workflow/decadal-prediction-scripts/R"
  source(file.path(cwd, "utils.R"))
}

config[["aggregation_period"]] = parse_config_aggregation_period(config)
match_season <- "DJFM"
models <- c(
  config$ensemble_data$cmip5$models,
  config$ensemble_data$cmip6$models
)

## Parse aggregation period specification
period <- config$aggregation_period[[aggregation_period]]
lead_tm <- period$lead_time
start <- min(lead_tm)
end <- max(lead_tm)
study_period <- period$study_period
n_years <- length(study_period)
lag <- period$lag # Boolean: should we make a lagged ensemble?
n_lag <- period$n_lag # If `lag` is true, how many lags to include?
if (!lag) n_lag <- 1
label <- period$name

## Only run if hindcast (TODO work out better way to do this)
if (period$hindcast) {

  outputdir = file.path(outputroot, aggregation_period, season)
  obs <- read_parquet(obspath)

  ## ################################### ##
  ## Prepare ensemble forecast           ##
  ## ################################### ##

  ensemble_fcst <- read_parquet(fcstpath)

  ## Take the ensemble mean
  ens_mean_group_vars <- c("init_year", "variable")
  fcst <-
    ensemble_fcst %>%
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
  fcst <-
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
  fcst <-
    fcst %>%
    mutate(
      obs_std = obs / sd(obs, na.rm=TRUE), # FIXME
      ens_mean_lag_std = ens_mean_lag / sd(ens_mean_lag, na.rm=TRUE),
      ens_mean_std = ens_mean / sd(ens_mean_lag, na.rm=TRUE)
    )

  ## Multiply standardised values by std dev
  ## of observed data so that the data ranges
  ## are comparable.
  fcst <-
    fcst %>%
    mutate(
      ens_mean_lag_var_adj = ens_mean_lag_std * sd(obs, na.rm=TRUE),
      ens_mean_var_adj = ens_mean_std * sd(obs, na.rm=TRUE)
    )

  ## Multiply by ACC (Note that standardising
  ## and multiplying by ACC is equivalent to
  ## multiplying by RPS)
  fcst <-
    fcst %>%
    mutate(
      corr = corr_cross_validate(ens_mean_lag_std, obs_std, 1),
    ) %>%
    mutate(
      ens_mean_lag_var_adj = corr * ens_mean_lag_var_adj,
      ens_mean_var_adj = corr * ens_mean_var_adj,
      ens_mean_lag_std = corr * ens_mean_lag_std,
      ens_mean_std = corr * ens_mean_std
    )
  ## NB ens_mean_lag_std is equivalent to init_nao_em on
  ## L318 in doug_smith_code.py

  ## The following method uses cross-validation to compute the
  ## standard deviations as well as ACC [not done in Doug Smith's code]
  ## fcst =
  ##   fcst %>%
  ##   mutate(
  ##     obs_std = obs / sd(obs, na.rm=TRUE), # FIXME
  ##   )
  ## ## Multiply by ACC (Note that standardising
  ## ## and multiplying by ACC is equivalent to
  ## ## multiplying by RPS)
  ## fcst =
  ##   fcst %>%
  ##   mutate(
  ##     sd_fcst = sd_cross_validate(ens_mean_lag, 1), # TEST
  ##     sd_obs = sd_cross_validate(obs, 1),           # TEST
  ##     corr = corr_cross_validate(ens_mean_lag, obs, 1),
  ##     rps = corr * sd_obs / sd_fcst
  ##   ) %>%
  ##   mutate(
  ##     ens_mean_lag_var_adj = ens_mean_lag * rps,
  ##     ens_mean_var_adj = ens_mean * rps,
  ##     ens_mean_lag_std = corr * ens_mean_lag / sd_fcst,
  ##     ens_mean_std = corr * ens_mean / sd_fcst
  ##   )
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
  ensemble_fcst <-
    ensemble_fcst %>%
    group_by_at(group_vars) %>%
    mutate(std = value / sd(value, na.rm=TRUE))

  ## FIXME - I think std is used to calculate error and
  ## nothing else, but check this

  ## Calculate error
  nao_matched_ensemble_fcst_error = calculate_error(
    fcst,
    ensemble_fcst,
    n_years,
    n_lag,
    "nao"
  )

  ## Join with ensemble data to get values
  ## FIXME - It's not clear why we need to add the error
  ## here, as we save this separately (below)
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
  nao_matched_ensemble_fcst_error <-
    nao_matched_ensemble_fcst %>%
    group_by(project, source_id, mip, member, init_year_matched, init_year) %>%
    summarize(any_na = any(is.na(value)), error = mean(error))
  write_parquet(
    nao_matched_ensemble_fcst_error,
    file.path(outputdir, "matched_ensemble_error.parquet")
  )
}
