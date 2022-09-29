## Author : Simon Moulds
## Date   : Nov-Dec 2021

library(tidyverse)
library(zoo)
library(RcppRoll)
library(lubridate)

check_file_exists = function(fn) {
  if (file.exists(fn)) {
    return(fn)
  } else {
    stop(sprintf("File %s does not exist", fn))
  }
}

parse_config_observed <- function(config, input_data_root) {
  items = config$observed
  giss = file.path(input_data_root, items$giss$subdirectory) %>% check_file_exists()
  gpcc = file.path(input_data_root, items$gpcc$subdirectory) %>% check_file_exists()
  hadcrut4 = file.path(input_data_root, items$hadcrut4$subdirectory) %>% check_file_exists()
  hadslp2r = file.path(input_data_root, items$hadslp2r$subdirectory) %>% check_file_exists()
  ncdc = file.path(input_data_root, items$ncdc$subdirectory) %>% check_file_exists()
  list(giss = giss, gpcc = gpcc, hadcrut4 = hadcrut4, hadslp2r = hadslp2r, ncdc = ncdc)
}

parse_config_ensemble <- function(config, input_data_root) {
  ## TODO
  return(config$ensemble)
}

parse_config_aux <- function(config, input_data_root) {
  ## TODO
  return(config$aux)
}

parse_config_output <- function(config) {
  ## TODO
  return(config$output)
}

parse_config_subset <- function(config) {
  if (!"subset" %in% names(config))
    return(NULL)

  nms = sapply(config$subset, FUN=function(x) x$name)
  config_subset =
    vector(mode = "list", length = length(nms)) %>%
    setNames(nms)

  for (i in 1:length(nms)) {
    item = config$subset[[i]]
    nm = item$name
    if (is.null(item$best_n)) {
      n = NA
      best_n = FALSE
    } else {
      n = as.integer(item$best_n)
      best_n = TRUE
    }
    if (is.null(item$worst_n)) {
      n = NA
      worst_n = FALSE
    } else {
      n = as.integer(item$worst_n)
      worst_n = TRUE
    }
    if (best_n & worst_n) {
      stop("Only one of `best_n` or `worst_n` may be supplied")
    }
    full = FALSE
    if (!any(best_n, worst_n)) {
      full = TRUE
    }
    projects = ifelse(is.null(item$projects), NA, item$projects)
    models = ifelse(is.null(item$models), NA, item$models)
    config_subset[[nm]] =
      list(name = nm,
           full = full,
           best_n = best_n,
           worst_n = worst_n,
           n = n,
           projects = projects,
           models = models)
  }
  config_subset
}

parse_config_aggregation_period <- function(config) {
  nms = sapply(config$aggregation_period, FUN=function(x) x$name)
  config_aggregation_period =
    vector(mode = "list", length = length(nms)) %>%
    setNames(nms)

  for (i in 1:length(nms)) {
    nm = nms[i]
    item = config$aggregation_period[[i]]
    error_name = ifelse(is.null(item$error_name), nm, item$error_name)
    lead_time = eval(str2lang(as.character(item$lead_time)))
    study_period = eval(str2lang(as.character(item$study_period)))
    if (is.null(lead_time)) stop("`lead_time` cannot be missing")
    if (is.null(study_period)) stop("`study_period` cannot be missing")
    observed = ifelse(is.null(item$observed), FALSE, as.logical(item$observed))
    hindcast = ifelse(is.null(item$hindcast), FALSE, as.logical(item$hindcast))
    lag = ifelse(is.null(item$lag), FALSE, as.logical(item$lag))
    n_lag = ifelse(is.null(item$n_lag), 4, as.integer(item$n_lag))
    constant_lead_time = ifelse(is.null(item$constant_lead_time), TRUE, as.logical(item$constant_lead_time))
    constant_fcst_period = ifelse(is.null(item$constant_fcst_period), FALSE, as.logical(item$constant_fcst_period))
    constant_lead_time = lag & constant_lead_time
    constant_fcst_period = lag & constant_fcst_period & !constant_lead_time
    config_aggregation_period[[nm]] =
      list(name = nm,
           error_name = error_name,
           lead_time = lead_time,
           study_period = study_period,
           observed = observed,
           hindcast = hindcast,
           lag = lag,
           n_lag = n_lag,
           constant_fcst_period = constant_fcst_period,
           constant_lead_time = constant_lead_time)
  }
  config_aggregation_period
}

parse_config_modelling <- function(config) {
  datasets = sapply(config$modelling, FUN=function(x) x$name)
  keys = c(
    "name", "input_dataset", "predictand",
    "model_family", "aggregation_periods",
    "study_period", "formulas",
    "sigma_formulas", "subsets"
  )
  config_modelling =
    vector(mode = "list", length = length(datasets)) %>%
    setNames(datasets)

  for (i in 1:length(datasets)) {
    dataset = datasets[[i]]
    opt_list =
      vector(mode = "list", length=length(keys)) %>%
      setNames(keys)

    for (j in 1:length(keys)) {
      key = keys[[j]]
      opt_list[[key]] = config$modelling[[i]][[key]]
    }
    opt_list$study_period = eval(str2lang(opt_list$study_period))
    opt_list$formulas = lapply(opt_list$formulas, FUN=function(x) as.formula(x))
    model_nms = names(opt_list$formulas)
    if (is.null(opt_list$sigma_formulas)) {
      sigma_formulas = list()
      for (j in 1:length(model_nms)) {
        nm = model_nms[j]
        sigma_formulas[[nm]] = as.formula("~1")
      }
    }
    opt_list$sigma_formulas = sigma_formulas
    if (!isTRUE(all.equal(model_nms, names(opt_list$sigma_formulas)))) {
      print(model_nms)
      print(names(opt_list$sigma_formulas))
      stop("model names not the same")
    }
    config_modelling[[dataset]] = opt_list
  }
  return(config_modelling)
}

parse_config <- function(config, input_data_root) {
  observed_section <- parse_config_observed(config, input_data_root)
  ensemble_section <- parse_config_ensemble(config, input_data_root)
  aux_data_section <- parse_config_aux(config, input_data_root)
  output_section <- parse_config_output(config)
  subset_section <- parse_config_subset(config)
  aggregation_period_section <- parse_config_aggregation_period(config)
  modelling_section <- parse_config_modelling(config)
  ## TODO checks (e.g. aggregation periods used in modelling section defined)
  list(observed_data = observed_section,
       ensemble_data = ensemble_section,
       aux_data = aux_data_section,
       output_data = output_section,
       subset = subset_section,
       aggregation_period = aggregation_period_section,
       modelling = modelling_section)
}

parse_config_io <- function(config, input_data_root) {
  observed_section <- parse_config_observed(config, input_data_root)
  ensemble_section <- parse_config_ensemble(config, input_data_root)
  aux_data_section <- parse_config_aux(config, input_data_root)
  output_section = parse_config_output(config)
  ## subset_section = parse_config_subset(config)
  aggregation_period_section = parse_config_aggregation_period(config)
  ## modelling_section = parse_config_modelling(config)
  ## TODO checks (e.g. aggregation periods used in modelling section defined)
  list(observed_data = observed_section,
       ensemble_data = ensemble_section,
       aux_data = aux_data_section,
       output_data = output_section,
       ## subset = subset_section,
       aggregation_period = aggregation_period_section)
       ## modelling = modelling_section)
}

get_obs <- function(filename, study_period, start = 2, end = 9) {
  ## Read raw observed data
  obs_raw = read_parquet(filename) #file.path(dir, "obs.parquet"))
  ## Pivot from long to wide
  obs_raw = obs_raw %>% pivot_wider(names_from=variable, values_from=value)
  ## Assign a reference year to DJFM and select this season
  obs = obs_raw %>%
    mutate(season_year = ifelse(month %in% c(1,2,3), year-1, year)) %>%
    filter(month %in% c(12, 1, 2, 3)) %>%
    group_by(season_year) %>%
    filter(n() == 4) # Only complete DJFM seasons
  vars = c("nao", "ea", "amv", "european_precip", "uk_precip", "uk_temp")

  ## Compute average seasonal values
  obs = obs %>% summarize(across(all_of(vars), mean))

  ## Calculate decadal means [function defined in `utils.R`]
  ## N.B. in this data frame season year is the year in which
  ## the start of the season falls [i.e. for DJFM it is the
  ## year of December]
  obs = rolling_fun(
    yrs = obs$season_year,
    data = obs,
    cols = vars,
    funs = mean,
    start = start, end = end
  )
  ## The result of the above function is that the value for each
  ## initialization year is the mean of the observed values for
  ## 2-9 years ahead. For example, the value assigned to the 1960
  ## initialization year is the average value for years 1961 to 1968.
  ## This essentially makes the observations comparable with the
  ## forecasts.

  ## Filter study period
  obs = obs %>% filter(init_year %in% study_period)
  compute_anomaly = function(x) x - mean(x, na.rm = TRUE)
  obs =
    obs %>%
    mutate(across(all_of(vars), compute_anomaly)) %>%
    ungroup()
  ## ## Convert to long format
  ## obs = obs %>% gather(variable, obs, -init_year)
  obs
}

get_hindcast_data <- function(dataset, study_period, lead_times) {
  ensemble_fcst_raw =
    open_dataset(dataset) %>%
    mutate(lead_time = season_year - init_year) %>%
    filter(lead_time %in% lead_times) %>%
    collect()
  ## Pivot from long to wide format
  ensemble_fcst_raw = ensemble_fcst_raw %>%
    pivot_wider(names_from=variable, values_from=value)
  ## Unit conversion
  ensemble_fcst_raw =
    ensemble_fcst_raw %>%
    mutate(nao = nao / 100) %>%
    mutate(ea = ea / 100) %>%
    mutate(european_precip = european_precip * 60 * 60 * 24) %>%
    mutate(uk_precip = uk_precip * 60 * 60 * 24)
  ## ## Correct initialisation years for GFDL data, which appear to be incorrect
  ## gfdl_index = ensemble_fcst_raw$source_id %in% "GFDL-CM2p1"
  ## ensemble_fcst_raw$init_year[gfdl_index] = ensemble_fcst_raw$init_year[gfdl_index] - 1
  ensemble_fcst_raw_complete =
    ensemble_fcst_raw %>%
    filter(init_year %in% study_period)
  ensemble_fcst_raw_complete
}

download_nrfa_data <- function(stn_id, metadata) {
  ## TODO tidy up this function, giving user more control over which variables are derived
  meta = metadata %>% filter(id %in% stn_id)
  ## Gauged daily flow [m3 s-1]
  gdf = get_ts(stn_id, "gdf") %>% as_tibble(rownames="time")
  ## Catchent daily rainfall [mm]
  cdr = try(get_ts(stn_id, "cdr") %>% as_tibble(rownames="time"))
  if (inherits(cdr, "try-error"))
    cdr <- gdf %>% rename(cdr = gdf) %>% mutate(cdr = NA)

  ## Create complete time series, in case the
  ## raw time series has missing values.
  start_date = gdf$time[1]
  end_date = gdf$time[nrow(gdf)]
  complete_ts = seq.POSIXt(
    as.POSIXct(start_date, tz="GMT", format="%Y-%m-%d"),
    as.POSIXct(end_date, tz="GMT", format="%Y-%m-%d"),
    by="1 day"
  ) %>% as.Date() %>% as.character()
  gdf =
    tibble(time = complete_ts) %>%
    left_join(gdf, by="time")
  availability = sum(!is.na(gdf$gdf)) / nrow(gdf) * 100
  df =
    gdf %>%
    left_join(cdr, by="time") %>%
    mutate(ID=stn_id, .after=time) %>%
    mutate(time = as.Date(time))
  ## TODO Filter out years with fewer than 330 days of records
  df =
    df %>%
    mutate(year = format(time, "%Y") %>% as.integer) %>%
    mutate(month = format(time, "%m") %>% as.integer)
  ## Neri et al [https://doi.org/10.1002/joc.5915]:
  ## "To avoid double counting the same event, we only consider
  ## one event in a window of +/- 5 days + logarithm of the
  ## drainage area"
  catchment_area = meta[["catchment_area"]]
  window_size = ((5 + log(catchment_area * 0.386102)) * 2) %>% round()
  ## Solari et al [https://doi.org/10.1002/2016WR019426]:
  ## "As the moving window travels through the series, each
  ## time that the data maximum in the window is located at
  ## its center, the maximum is regarded as a peak"
  df =
    df %>% # Neri et al
    mutate(pot = roll_max(gdf, n=window_size, align="center", fill=NA)) %>%
    mutate(is_peak = Vectorize(isTRUE)(pot == gdf)) %>%
    mutate(peak = ifelse(is_peak, pot, NA))
  ## Unsure whether we need to make this unique or not?
  peaks = df$pot[df$is_peak] ##%>% unique()
  n_years = length(df$year %>% unique())
  peaks_sorted = sort(peaks, decreasing = TRUE)
  threshold_1 = peaks_sorted[n_years]
  threshold_2 = peaks_sorted[n_years * 2]
  threshold_3 = peaks_sorted[n_years * 3]
  threshold_4 = peaks_sorted[n_years * 4]

  df =
    df %>%
    mutate(pot_1 = ifelse(Vectorize(isTRUE)(peak >= threshold_1), 1, 0)) %>%
    mutate(pot_2 = ifelse(Vectorize(isTRUE)(peak >= threshold_2), 1, 0)) %>%
    mutate(pot_3 = ifelse(Vectorize(isTRUE)(peak >= threshold_3), 1, 0)) %>%
    mutate(pot_4 = ifelse(Vectorize(isTRUE)(peak >= threshold_4), 1, 0))

  ## Add climate season label to rows
  df =
    df %>%
    mutate(
      clim_season = case_when(
        month %in% c(12, 1, 2, 3) ~ "DJFM",
        month %in% c(4, 5) ~ "AM",
        month %in% c(6, 7, 8, 9) ~ "JJAS",
        month %in% c(10, 11) ~ "ON"
      )
    ) %>%
    mutate(season_year = ifelse(month %in% c(1, 2, 3), year - 1, year))

  ## Summarize to get flood counts
  df =
    df %>%
    group_by(ID, clim_season, season_year) %>%
    summarize(
      missing_pct = (sum(is.na(gdf)) / n()) * 100,
      Q_max = max(gdf, na.rm = TRUE),
      Q_mean = mean(gdf, na.rm = TRUE),
      Q_05 = quantile(gdf, probs = 0.05, na.rm = TRUE, names = FALSE),
      Q_50 = quantile(gdf, probs = 0.50, na.rm = TRUE, names = FALSE),
      Q_90 = quantile(gdf, probs = 0.90, na.rm = TRUE, names = FALSE),
      Q_95 = quantile(gdf, probs = 0.95, na.rm = TRUE, names = FALSE),
      P_sum = sum(cdr, na.rm = TRUE),
      POT_1 = sum(pot_1, na.rm = TRUE),
      POT_2 = sum(pot_2, na.rm = TRUE),
      POT_3 = sum(pot_3, na.rm = TRUE),
      POT_4 = sum(pot_4, na.rm = TRUE),
    ) %>%
    ungroup() %>%
    mutate(across(Q_max:POT_4, ~ifelse(is.finite(.), ., NA)))
}
## ## Set some default formatting options for ggplot
## plot_format_objects <- list(
##   scale_x_continuous(breaks=seq(1960, 2010, 10), limits=c(1960, 2005)),
##   labs(x="Start of 8-year period", y="NAO anomaly (hPa)")
## )

rolling_fun <- function(yrs, data, cols, funs, start=2, end=9) {
  ## Compute rolling n-year means.
  ##
  ## Args:
  ##   yrs   : integer. Year index
  ##   data  : data.frame.
  ##   cols  : character.
  ##   funs  : function or list of functions.
  ##   start : integer. Index of start point (where index
  ##           of the point for which the rolling mean is
  ##           being computed is 1).
  ##   end   : integer. Index of end point.
  ##
  ## Return:
  ##   Data frame
  if (is.list(funs) & (!isTRUE(all(cols == names(funs)))))
    stop()
  ## Preallocate output
  out = lapply(cols, FUN=function(x) rep(NA, length(yrs)))
  names(out) = cols
  out = c(list(init_year = yrs), out)

  ## Also collect the start year of the time window
  period_start = rep(NA, length(yrs))
  for (i in 1:(length(yrs)-end+1)) {
  ## for (i in 1:(length(yrs)-end)) {
    start_index = i + start - 1
    end_index = i + end - 1
    for (j in 1:length(cols)) {
      nm = cols[j]
      x = data[[nm]]
      if (is.function(funs)) {
        fun = funs
      } else {
        fun = funs[[nm]]
      }
      ## fun = ifelse(is.function(funs), funs, funs[[nm]])
      out[[nm]][i] = fun(x[start_index:end_index], na.rm=TRUE)
    }
    ## period_start[i] = yrs[start_index]
  }
  out = out %>% as.data.frame() %>% as_tibble()
  out
}

corr_cross_validate <- function(fcst, obs, leave_out_add=0) {
  ## Compute cross validated anomaly correlation.
  ##
  ## Args:
  ##   fcst : numeric. Forecast data
  ##   obs  : numeric. Observed data
  ##   n    : integer. Number of additional
  ##          points to leave out
  ##
  ## Return:
  ##   Numeric
  ntimes = length(obs)
  index = seq(1, ntimes) # for selection
  corr = rep(0., ntimes)
  for (i in 1:ntimes) {
    if (leave_out_add > 0) {
      leave_out_start = max(c(1, i - leave_out_add))
      leave_out_end = min(c(ntimes, i + leave_out_add))
      leave_out_index = seq(leave_out_start, leave_out_end)
    } else {
      leave_out_index = c()
    }
    keep_index = !(index %in% leave_out_index)
    corr_this = cor.test(
      fcst[keep_index],
      obs[keep_index],
      method="pearson", alternative="greater"
    )
    corr[i] = unname(corr_this$estimate)
  }
  corr
}

calculate_error <- function(fcst, ensemble_fcst, n_years, n_forecast, match_var) {
  ## Match ensemble members.
  ##
  ## Select n members by comparing with variance-adjusted
  ## ensemble mean data.
  ##
  ## Args:
  ##   fcst : data.frame. Ensemble mean forecast data.
  ##   ensemble_fcst : data.frame. Ensemble fcst data.
  ##   n_years       : integer. Number of years in the study period.
  ##   n_forecast    : integer. Number of forecasts to use to create
  ##                   lagged forecast.
  ##   match_var     : character. Variable to use in matching
  ##                   algorithm.
  ##   n_select      : integer. Number of members to select.
  ##   best          : bool. Whether to select the `n_select` best
  ##                   performing or the `n_select` worst performing.
  ##
  ## Return:
  ##   Data frame.

  ## Select only the forecast data for the variable against
  ## which we will perform the matching
  match_fcst =
    fcst %>%
    filter(variable %in% match_var) %>%
    rename(init_year_lag = init_year)

  ## When the ensemble data is lagged we need to calculate
  ## the absolute error between the variance-adjusted ensemble
  ## mean for time point i and the individual ensemble members
  ## which contribute to the lagged value, i.e. time point i,
  ## i-1, ..., i-n, where n is the total lag (3 in our case).
  ## Here we create a data frame which we can join with
  ## `ensemble_fcst` to account for this.
  init_year_lag =
    data.frame(init_year = study_period) %>%
    slice(
      rep(1:n_forecast, n_years - n_forecast + 1)
      + (rep(1:(n_years - n_forecast + 1), each = n_forecast) - 1)
    ) %>%
    as_tibble()
  ## init_year_lag

  ## Now we create a temporary grouping variable to compute the
  ## year to which the lagged values contribute (e.g. the lagged
  ## value for year 1963 is the mean of values for years 1960,
  ## 1961, 1962, 1963). As we lag backwards this is simply the latest
  ## year in each group. We call this year `init_year_lag`
  init_year_lag =
    init_year_lag %>%
    mutate(tmp_group = rep(1:(n() / n_forecast), each = n_forecast)) %>%
    group_by(tmp_group) %>%
    mutate(init_year_lag = max(init_year)) %>%
    ungroup() %>%
    dplyr::select(-tmp_group)

  ## ## Have a look at a few examples:
  ## init_year_lag %>% filter(init_year_lag %in% 1965)
  ## init_year_lag %>% filter(init_year_lag %in% 1983)
  ## init_year_lag %>% filter(init_year_lag %in% 2005)

  ## Join the data frame with ensemble_fcst, which approximately
  ## quadruples its size (values are duplicated when they are
  ## included in a different lag year (i.e. `init_year_lag`))
  ensemble_fcst_lag =
    ensemble_fcst %>%
    left_join(init_year_lag, by = "init_year") %>%
    arrange(source_id, member, init_year_lag, init_year)

  ## Join the datasets
  ensemble_fcst_lag =
    ensemble_fcst_lag %>%
    left_join(match_fcst, by = c("variable", "init_year_lag"))

  ## Calculate mean absolute error for the matching variable
  ensemble_fcst_lag =
    ensemble_fcst_lag %>%
    filter(variable %in% match_var) %>%
    mutate(error = abs(std - ens_mean_lag_std))

  ## OLD:

  ## ## For each lagged year, select the n best performing members:
  ## nao_matched_ensemble_fcst =
  ##   ensemble_fcst_lag %>%
  ##   group_by(source_id, member, init_year_lag) %>%
  ##   mutate(min_error_year = init_year[which.min(error)]) %>%
  ##   filter(init_year %in% min_error_year) %>%
  ##   ## summarise(error = min(error, na.rm=TRUE)) #%>%
  ##   ungroup() %>%
  ##   group_by(init_year_lag) ## %>%
  ## ## slice_min(error, n=n_select)

  ## if (best) {
  ##   nao_matched_ensemble_fcst =
  ##     nao_matched_ensemble_fcst %>%
  ##     slice_min(error, n=n_select)
  ## } else {
  ##   nao_matched_ensemble_fcst =
  ##     nao_matched_ensemble_fcst %>%
  ##     slice_max(error, n=n_select)
  ## }

  ## ## Tidy up
  ## nao_matched_ensemble_fcst =
  ##   nao_matched_ensemble_fcst %>%
  ##   dplyr::select(project, source_id, mip, member, init_year, init_year_lag)

  ## ## Join with original ensemble data
  ## nao_matched_ensemble_fcst =
  ##   nao_matched_ensemble_fcst %>%
  ##   left_join(ensemble_fcst) %>%
  ##   rename(init_year_matched = init_year) %>%
  ##   rename(init_year = init_year_lag)

  ## ## OLD 2:

  ## ## For each lagged year, select the n best performing members:
  ## nao_matched_ensemble_fcst =
  ##   ensemble_fcst_lag %>%
  ##   group_by(source_id, member, init_year_lag) %>%
  ##   mutate(min_error_year = init_year[which.min(error)]) %>%
  ##   filter(init_year %in% min_error_year) %>%
  ##   ## summarise(error = min(error, na.rm=TRUE)) #%>%
  ##   ungroup() ## %>%
  ##   ## group_by(init_year_lag) ## %>%
  ## ## slice_min(error, n=n_select)

  ## nao_matched_ensemble_fcst =
  ##   nao_matched_ensemble_fcst %>%
  ##   dplyr::select(project, source_id, mip, member, init_year, init_year_lag, error)

  ## ## Join with original ensemble data
  ## nao_matched_ensemble_fcst =
  ##   nao_matched_ensemble_fcst %>%
  ##   ## left_join(ensemble_fcst) %>% # join on init_year
  ##   rename(init_year_matched = init_year) %>%
  ##   rename(init_year = init_year_lag)

  ## ## Return matched forecast
  ## nao_matched_ensemble_fcst

  ## NEW:

  ## For each lagged year, select the n best performing members:
  nao_matched_ensemble_fcst =
    ensemble_fcst_lag %>%
    ## group_by(source_id, member, init_year_lag) %>%
    ## mutate(min_error_year = init_year[which.min(error)]) %>%
    ## filter(init_year %in% min_error_year) %>%
    ## ## summarise(error = min(error, na.rm=TRUE)) #%>%
    ungroup() ## %>%
    ## group_by(init_year_lag) ## %>%
  ## slice_min(error, n=n_select)

  nao_matched_ensemble_fcst =
    nao_matched_ensemble_fcst %>%
    dplyr::select(project, source_id, mip, member, init_year, init_year_lag, error)

  ## Join with original ensemble data
  nao_matched_ensemble_fcst =
    nao_matched_ensemble_fcst %>%
    ## left_join(ensemble_fcst) %>% # join on init_year
    rename(init_year_matched = init_year) %>%
    rename(init_year = init_year_lag)

  ## Return matched forecast
  nao_matched_ensemble_fcst
}

## create_annual_ensemble_forecast = function(ensemble_fcst_error,
##                                            ensemble_fcst_raw,
##                                            vars = climate_vars,
##                                            model_select = NA,
##                                            project_select = NA,
##                                            full = TRUE,
##                                            best_n = FALSE,
##                                            worst_n = FALSE,
##                                            n_select = 20,
##                                            lead_times = 2,
##                                            lag = TRUE) {

##   stop("Not yet implemented")
##   if (length(lead_times) == 1) {
##     if (lag) {
##       lead_times = seq(lead_times, lead_times + 3)
##       ensemble_fcst_raw = ensemble_fcst_raw %>% filter(lead_time %in% lead_times)
##       ensemble_fcst_raw = ensemble_fcst_raw %>% group_by(season_year)
##     }
##   }
## }

## create_multiyear_ensemble_forecast = function() {}
create_ensemble_forecast <- function(ensemble_fcst_error,
                                    ensemble_fcst,
                                    vars = climate_vars,
                                    model_select = NA,
                                    project_select = NA,
                                    full = TRUE,
                                    best_n = FALSE,
                                    worst_n = FALSE,
                                    n_select = 20) {
                                    ## lead_times = c(2:9),
                                    ## lag = TRUE,
                                    ## n_lag = 4) {

  ## Filter by models
  models = ensemble_fcst$source_id %>% unique() %>% toupper()
  if (!is.na(model_select)) {
    model_select = toupper(model_select)
    if (!all(model_select %in% models)) {
      stop(paste0("Invalid model specification. Valid models are:\n", paste0(models, collapse = ", ")))
    }
    ensemble_fcst =
      ensemble_fcst %>%
      filter(toupper(source_id) %in% model_select)
  }

  ## Filter by projects (i.e. CMIP5/6)
  projects = ensemble_fcst$project %>% unique() %>% toupper()
  if (!is.na(project_select)) {
    project_select = toupper(project_select)
    if (!all(project_select %in% projects)) {
      stop(paste0("Invalid project specification. Valid projects are:\n", paste0(projects, collapse = ", ")))
    }
    ensemble_fcst =
      ensemble_fcst %>%
      filter(toupper(project) %in% project_select)
  }

  # Only allow models without NA values [notable at present is CESM1-1-CAM5]
  nao_matched_ensemble_fcst_subset =
    ensemble_fcst_error %>%
    filter(!any_na)

  ## Select n best (worst) performing members
  if (best_n | worst_n) {
    ## `init_year_matched` is the actual initialization year of the member
    ## `init_year` is the initialization year of the current forecast
    slice_fun = slice_min
    if (worst_n) {
      slice_fun = slice_max
    }
    nao_matched_ensemble_fcst_subset =
      nao_matched_ensemble_fcst_subset %>%
      group_by(init_year) %>%
      slice_fun(error, n = n_select)
  }
  nao_matched_ensemble_fcst_subset =
    nao_matched_ensemble_fcst_subset %>%
    dplyr::select(-any_na, -error)

  ## Join with forecast data
  ensemble_fcst = ensemble_fcst %>% dplyr::select(-std) %>% pivot_wider(names_from = variable, values_from = value)
  ensemble_fcst =
    nao_matched_ensemble_fcst_subset %>%
    left_join(ensemble_fcst, by = c("project", "mip", "source_id", "member", "init_year_matched", "init_year"))

  ## Compute ensemble mean for each year
  ensemble_group_vars = c("init_year")
  ensemble_fcst =
    ensemble_fcst %>%
    group_by_at(ensemble_group_vars) %>%
    summarize(across(all_of(vars), mean, na.rm = TRUE))
  ensemble_fcst
}

## create_ensemble_forecast_old = function(ensemble_fcst_error,
##                                     ensemble_fcst_raw,
##                                     vars = climate_vars,
##                                     model_select = NA,
##                                     project_select = NA,
##                                     full = TRUE,
##                                     best_n = FALSE,
##                                     worst_n = FALSE,
##                                     n_select = 20,
##                                     lead_times = c(2:9),
##                                     lag = TRUE,
##                                     n_lag = 4) {

##   ## ## Filter by lead time
##   ## ensemble_fcst_raw =
##   ##   ensemble_fcst_raw %>%
##   ##   filter(lead_time %in% lead_times)
##   ## anomaly_group_vars = c("source_id", "member", "lead_time")
##   ## anomalyfun = function(x) x - mean(x, na.rm = TRUE)
##   ## ensemble_fcst_raw =
##   ##   ensemble_fcst_raw %>%
##   ##   group_by_at(anomaly_group_vars) %>%
##   ##   mutate(across(all_of(vars), anomalyfun))

##   ## Filter by models
##   models = ensemble_fcst_raw$source_id %>% unique() %>% toupper()
##   if (!is.na(model_select)) {
##     model_select = toupper(model_select)
##     if (!all(model_select %in% models)) {
##       stop(paste0("Invalid model specification. Valid models are:\n", paste0(models, collapse = ", ")))
##     }
##     ensemble_fcst_raw = ensemble_fcst_raw %>% filter(toupper(source_id) %in% model_select)
##   }

##   ## Filter by projects (i.e. CMIP5/6)
##   projects = ensemble_fcst_raw$project %>% unique() %>% toupper()
##   if (!is.na(project_select)) {
##     project_select = toupper(project_select)
##     if (!all(project_select %in% projects)) {
##       stop(paste0("Invalid project specification. Valid projects are:\n", paste0(projects, collapse = ", ")))
##     }
##     ensemble_fcst_raw = ensemble_fcst_raw %>% filter(toupper(project) %in% project_select)
##   }

##   # Only allow models without NA values [notable at present is CESM1-1-CAM5]
##   nao_matched_ensemble_fcst_subset =
##     ensemble_fcst_error %>%
##     filter(!any_na)

##   ## Select n best (worst) performing members
##   if (best_n | worst_n) {
##     ## `init_year_matched` is the actual initialization year of the member
##     ## `init_year` is the initialization year of the current forecast
##     slice_fun = slice_min
##     if (worst_n) {
##       slice_fun = slice_max
##     }
##     nao_matched_ensemble_fcst_subset =
##       nao_matched_ensemble_fcst_subset %>%
##       group_by(init_year) %>%
##       slice_fun(error, n = n_select)
##   }
##   nao_matched_ensemble_fcst_subset =
##     nao_matched_ensemble_fcst_subset %>%
##     dplyr::select(-any_na, -error)

##   stop()

##   ## Aggregate over lead times (multi-year forecast only)
##   if (length(lead_times) > 1) {
##     group_vars = c(
##       "project", "mip", #"experiment",
##       "source_id", "member", "init_year"
##     )
##     ## Average over multi-year period
##     ensemble_fcst_raw =
##       ensemble_fcst_raw %>%
##       filter(lead_time %in% lead_times) %>%
##       group_by_at(group_vars) %>%
##       summarize(across(all_of(vars), mean))
##     ## Compute anomaly
##     anomaly_group_vars = c("source_id", "member")
##   } else {
##     lead_times = seq(lead_tm, lead_tm + n_lag - 1)
##     ensemble_fcst_raw =
##       ensemble_fcst_raw %>%
##       filter(lead_time %in% lead_times)
##     anomaly_group_vars = c("source_id", "member", "lead_time")
##   }
##   anomalyfun = function(x) x - mean(x, na.rm = TRUE)
##   ensemble_fcst_raw =
##     ensemble_fcst_raw %>%
##     group_by_at(anomaly_group_vars) %>%
##     mutate(across(all_of(vars), anomalyfun))

##   stop()

##   ## Rename init_year prior to merging with subset
##   ensemble_fcst_raw =
##     ensemble_fcst_raw %>%
##     rename(init_year_matched = init_year)

##   ## Join with NAO-matched members
##   if (lag) {
##     ensemble_fcst_raw =
##       nao_matched_ensemble_fcst_subset %>%
##       ## filter(lag %in% 1) %>%
##       left_join(
##         ensemble_fcst_raw,
##         by=c("project", "mip", "source_id", "member", "init_year_matched")
##       )
##   } else {
##     ensemble_fcst_raw =
##       nao_matched_ensemble_fcst_subset %>%
##       filter(lag %in% 1) %>%
##       left_join(
##         ensemble_fcst_raw,
##         by=c("project", "mip", "source_id", "member", "init_year_matched")
##       )
##   }

##   ## ## Recalculate seasonal_averageson_year if output is not aggregated
##   ## if (!aggregate) {
##   ##   ensemble_fcst_raw =
##   ##     ensemble_fcst_raw %>%
##   ##     mutate(season_year = init_year + lead_time)
##   ## }

##   ensemble_group_vars = c("init_year")

##   ## Compute ensemble mean for each year
##   ensemble_fcst_raw =
##     ensemble_fcst_raw %>%
##     group_by_at(ensemble_group_vars) %>%
##     summarize(across(all_of(vars), mean, na.rm = TRUE))

##   ensemble_fcst_raw
## }

compute_quantiles <- function(newdata,
                             data,
                             ...,
                             quantiles = c(0.5, 0.25, 0.75, 0.05, 0.95),
                             model_family) {
  ## Create centiles for a list of models
  ##
  ## Args:
  ##   data:         data.frame. Data to use for model prediction
  ##   ...:          model objects
  ##   quantiles:    numeric. quantiles to compute.
  ##   model_family: character. Currently only GA and PO are supported.
  ##
  ## Returns:
  ##   Tibble.

  ## Extract model names from dots
  dots = match.call(expand.dots=FALSE)$...
  model_nms = names(dots)
  models = list(...)
  names(models) = model_nms
  id_cols = c("ID", "clim_season", "year", "lead_time")
  quantile_names = paste0("Q", formatC(quantiles * 100, width=2, flag=0))
  computed_quantiles = list()
  for (i in 1:length(model_nms)) {
    nm = model_nms[i]
    model = models[[nm]]
    tbl =
      rep(NA, length(quantile_names) + 1) %>%
      setNames(c("model", quantile_names)) %>%
      as.list() %>% as_tibble()
    tbl[["model"]] = nm
    if (!is.null(model)) {
      mu = predict(
        models[[nm]],
        newdata = newdata,
        what = "mu",
        type = "response",
        data = data
      )
      if (model_family %in% c("GA")) { # TODO add more cases
        sigma = predict(
          models[[nm]],
          newdata = newdata,
          what = "sigma",
          type = "response",
          data = data
        )
      }
      for (j in 1:length(quantile_names)) {
        q = quantiles[j]
        qnm = quantile_names[j]
        if (model_family == "GA") {
          tbl[[qnm]] = qGA(q, mu, sigma)
        } else if (model_family == "PO") {
          tbl[[qnm]] = qPO(q, mu)
        }
      }
    }
    tbl = cbind(
      tbl,
      newdata %>% dplyr::select(any_of(id_cols))
    )
    computed_quantiles[[i]] = tbl
  }
  computed_quantiles =
    do.call("rbind", computed_quantiles) %>%
    as_tibble()
  computed_quantiles
}

fit_models_cv <- function(x,
                          config_section,
                          lead_time,
                          ...) {
  n_fold <- nrow(x)
  for (p in 1:n_fold) {
    ## lead_time <- config_section[[label]]$lead_time
    buffer <- length(lead_time) - 1
    idx = seq_len(n_fold)
    test_idx <- p
    remove_idx <- seq(p - buffer, p + buffer)
    remove_idx <- remove_idx[remove_idx > 0]
    train_idx = idx[!idx %in% remove_idx]
    train_data = x[train_idx,]
    test_data = x[test_idx,]

    ## N.B. because we use a function to fit models the data.frame
    ## used in the original fit is not stored in the model call
    ## correctly. This means that when we call `predict.gamlss`
    ## we have to supply the data.frame used for the original
    ## fit (i.e. train_data)
    models = fit_models(
      config_section$formulas,
      config_section$sigma_formulas,
      config_section$model_family,
      train_data
    )
    prediction = do.call(
      "compute_quantiles",
      c(list(newdata = test_data, data = train_data), models, list(model_family=config_section$model_family))
    )
    ## Create output data frame
    obs = test_data[["Q"]] #config_section$predictand]]
    ## exp = test_data[["Q50"]]
    obs_column_name = paste0(config_section$predictand, '_obs')
    exp_column_name = paste0(config_section$predictand, '_exp')
    prediction =
      prediction %>%
      mutate(
        !!obs_column_name := obs,
        !!exp_column_name := Q50
      ) %>%
      mutate(
        date = NA,
        model = paste0("GAMLSS_", model)
      ) %>%
      dplyr::select(obs_column_name, exp_column_name, year, model)
    ## obs = test_data[[config_section$predictand]]
    ## prediction =
    ##   prediction %>%
    ##   mutate(
    ##     predictand = config_section$predictand,
    ##     obs = obs,
    ##     exp = Q50 #!!ensym(config_section$predictand)
    ##     ## exp = !!ensym(config_section$predictand)
    ##   )
    pred_idx = length(catchment_prediction_list) + 1
    catchment_prediction_list[[pred_idx]] = prediction
  }
  if (length(catchment_prediction_list) == 0) {
    catchment_prediction <- NULL
  } else {
    catchment_prediction = do.call("rbind", catchment_prediction_list)
  }
  catchment_prediction
}

fit_models_forward_chain <- function(x,
                                     config_section,
                                     training_period_start,
                                     training_period_end,
                                     test_period_end,
                                     lead_time,
                                     ...) {
  catchment_prediction_list = list()
  buffer <- length(lead_time)
  max_training_period_end = test_period_end - buffer #window
  while (training_period_end <= max_training_period_end) {
    training_years <- seq(training_period_start, training_period_end)
    test_year <- max(training_years) + buffer
    train_data <- x %>% filter(year %in% training_years)
    test_data <- x %>% filter(year %in% test_year)
    ## Only continue if there is sufficient training/testing data
    if (nrow(train_data) == 0 | nrow(test_data) == 0) {
      training_period_end <- training_period_end + 1
      next
    }
    ## N.B. because we use a function to fit models the data.frame
    ## used in the original fit is not stored in the model call
    ## correctly. This means that when we call `predict.gamlss`
    ## we have to supply the data.frame used for the original
    ## fit (i.e. train_data)
    models = fit_models(
      config_section$formulas,
      config_section$sigma_formulas,
      config_section$model_family,
      train_data
    )
    prediction = do.call(
      "compute_quantiles",
      c(list(newdata = test_data, data = train_data), models, list(model_family=config_section$model_family))
    )
    ## Create output data frame
    obs = test_data[["Q"]] #config_section$predictand]]
    ## exp = test_data[["Q50"]]
    obs_column_name = paste0(config_section$predictand, '_obs')
    exp_column_name = paste0(config_section$predictand, '_exp')
    prediction =
      prediction %>%
      mutate(
        !!obs_column_name := obs,
        !!exp_column_name := Q50
      ) %>%
      ## mutate(
      ##   date = NA,
      ##   model = paste0("GAMLSS_", model)
      ## ) %>%
      dplyr::select(all_of(c(obs_column_name, exp_column_name)), year, model)

    pred_idx = length(catchment_prediction_list) + 1
    catchment_prediction_list[[pred_idx]] = prediction
    training_period_end <- training_period_end + 1
  }
  if (length(catchment_prediction_list) == 0) {
    catchment_prediction <- NULL
  } else {
    catchment_prediction = do.call("rbind", catchment_prediction_list)
  }
  catchment_prediction
}

fit_models <- function(formulas,
                       sigma_formulas,
                       model_family,
                       data,
                       ...) {

  model_names = names(formulas)
  if (!isTRUE(all.equal(model_names, names(sigma_formulas)))) {
    stop()
  }
  fitted_models = list()
  for (i in 1:length(model_names)) {
    model_nm = model_names[i]
    model = try(
      gamlss(
        formula = formulas[[i]],
        sigma.formula = sigma_formulas[[i]],
        family = model_family,
        trace = FALSE,
        data = data,
        ...
      )
    )
    if (inherits(model, "try-error")) model = NULL
    fitted_models[[model_nm]] = model
  }
  fitted_models
}

get_aic <- function(model_list) {
  aic = sapply(
    models,
    FUN=function(x) ifelse(is.null(x), NA, AIC(x))
  ) %>% as_tibble_row()
  aic
}

mean_square_error_skill_score <- function(obs, exp) {
  ## MSSS
  mse = mean((exp - obs) ^ 2)
  mse_ref = mean((mean(obs) - obs) ^ 2)
  msss = 1 - (mse / mse_ref)
  ## ACC
  acc = cor(obs, exp, method = "pearson")
  ## correlation
  r = cor(obs, exp, method = "pearson")
  ## potential skill [= coefficient of determination]
  ps = r ^ 2
  ## slope reliability
  srel = (r - (sd(exp) / sd(obs))) ^ 2
  ## standardized mean error
  sme = ((mean(exp) - mean(obs)) / sd(obs)) ^ 2
  ## msss = ps - srel - sme
  list(msss = msss, ps = ps, srel = srel, sme = sme, acc = acc)
}
