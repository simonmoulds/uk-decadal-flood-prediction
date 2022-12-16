## Author : Simon Moulds
## Date   : Nov-Dec 2021

library(tidyverse)
library(zoo)
library(RcppRoll)
library(lubridate)
library(scoringRules)
library(SpecsVerification)

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
  return(config$ensemble)
}

parse_config_aux <- function(config, input_data_root) {
  return(config$aux)
}

parse_config_output <- function(config) {
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

get_obs <- function(filename,
                    study_period,
                    start = 2,
                    end = 9) {
  ## Read raw observed data
  obs_raw = read_parquet(filename)
  ## Pivot from long to wide
  obs_raw = obs_raw %>% pivot_wider(names_from=variable, values_from=value)
  ## Assign a reference year to DJFM and select this season
  obs = obs_raw %>%
    mutate(season_year = ifelse(month %in% c(1,2,3), year-1, year)) %>%
    filter(month %in% c(12, 1, 2, 3)) %>%
    group_by(season_year) %>%
    filter(n() == 4) # Only complete DJFM seasons

  ## ## TESTING
  ## winter <- which.min(months) != 1
  ## next_year_months <- seq(min(months), months[length(months)])
  ## obs = obs_raw %>%
  ##   filter(month %in% c(12, 1, 2, 3)) #%>%
  ##   mutate(season_year = ifelse(month %in% c(1,2,3), year-1, year)) %>%
  ##   group_by(season_year) %>%
  ##   filter(n() == 4) # Only complete DJFM seasons
  ## ## END TESTING

  ## Compute average seasonal values
  vars = c("nao", "ea", "amv", "european_precip", "uk_precip", "uk_temp")
  obs = obs %>% summarize(across(starts_with(vars), mean))

  obs_antecedent <- obs_raw %>%
    filter(month %in% c(9, 10, 11)) %>%
    mutate(season_year = year) %>%
    group_by(season_year) #%>%
    ## filter(n() == 3) # Only complete SON seasons

  antecedent_vars <- c("european_precip", "uk_temp", "uk_precip")
  obs_antecedent <- obs_antecedent %>%
    summarize(across(starts_with(antecedent_vars), mean)) %>%
    rename_at(vars(starts_with(antecedent_vars)), function(x) paste0(x, "_antecedent"))

  obs <- obs %>% left_join(obs_antecedent, by="season_year")
  all_vars <- names(obs)
  all_vars <- all_vars[!all_vars %in% "season_year"]
  ## obs = obs %>% summarize(across(all_of(vars), mean))
  ## all_antecedent_vars <- names(obs_antecedent)
  ## all_antecedent_vars <- all_antecedent_vars[!all_antecedent_vars %in% "season_year"]

  ## Calculate decadal means [function defined in `utils.R`]
  ## N.B. in this data frame season year is the year in which
  ## the start of the season falls [i.e. for DJFM it is the
  ## year in which December falls
  obs = rolling_fun(
    yrs = obs$season_year,
    data = obs,
    ## cols = vars,
    cols = all_vars,
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
    ## mutate(across(all_of(vars), compute_anomaly)) %>%
    mutate(across(all_of(all_vars), compute_anomaly)) %>%
    ungroup()
  ## ## Convert to long format
  ## obs = obs %>% gather(variable, obs, -init_year)
  obs
}

get_obs_new <- function(dataset,
                        study_period,
                        start = 2,
                        end = 9,
                        vars = c("nao"),
                        months = c(12, 1, 2, 3)) {

  ## Read raw observed data
  ## obs_raw = read_parquet(filename)
  ## Pivot from long to wide
  obs_raw = dataset %>% pivot_wider(names_from=variable, values_from=value)
  ## Select season
  obs <- obs_raw %>% filter(month %in% months)
  if (which.min(months) != 1) {
    ## This adjusts season_year if the season covers two years (e.g.DJF)
    next_year_months <- seq(min(months), months[length(months)])
    obs <-
      obs %>%
      mutate(season_year = ifelse(month %in% next_year_months, year-1, year))
  } else {
    obs <- obs %>% mutate(season_year = year)
  }
  ## Restrict to complete seasons
  obs <- obs %>% group_by(season_year) %>% filter(n() == length(months))

  ## Compute average seasonal values
  obs = obs %>% summarize(across(starts_with(vars), mean))

  ## Calculate decadal means [function defined in `utils.R`]
  ## N.B. in this data frame season year is the year in which
  ## the start of the season falls [i.e. for DJFM it is the
  ## year in which December falls
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

  ## Convert to long format
  obs =
    obs %>%
    pivot_longer(starts_with(vars), names_to = "variable", values_to = "obs")
  obs
}

ensemble_fcst_unit_conversion <- function(x) {
  x <- x %>%
    mutate(across(matches("^nao$"), function(x) x / 100)) %>%
    mutate(across(matches("^ea$"), function(x) x / 100)) %>%
    mutate(across(contains("precip"), function(x) x * 60 * 60 * 24))
  x
}

get_hindcast_data_new <- function(dataset,
                                  study_period,
                                  lead_times,
                                  vars,
                                  months) {

  ## ens_fcst_raw = dataset %>% pivot_wider(names_from=variable, values_from=value)
  ## Select season
  ens_fcst <- dataset %>% filter(month %in% months)
  if (which.min(months) != 1) {
    ## This adjusts season_year if the season covers two years (e.g.DJF)
    next_year_months <- seq(min(months), months[length(months)])
    ens_fcst <-
      ens_fcst %>%
      mutate(season_year = ifelse(month %in% next_year_months, year-1, year))
  } else {
    ens_fcst <- ens_fcst %>% mutate(season_year = year)
  }
  ens_fcst <- ens_fcst %>% collect()
  ens_fcst <- ens_fcst %>% pivot_wider(names_from=variable, values_from=value)
  ## Restrict to complete seasons
  group_vars <- c("season_year", "project", "mip", "source_id", "member", "init_year")
  ens_fcst <- ens_fcst %>%
    group_by_at(group_vars) %>%
    filter(n() == length(months)) %>%
    arrange(init_year, member, source_id, year, month)

  ## ## Compute average seasonal values
  ## ens_fcst <- ens_fcst %>% group_by_at(group_vars) %>% summarize(value = mean(value), n = n())

  ens_fcst <- ens_fcst %>% summarize(across(starts_with(vars), mean))
  ens_fcst <-
    ens_fcst %>%
    mutate(lead_time = season_year - init_year + 1) %>%
    filter(lead_time %in% lead_times)

  ## Select data for the study period
  ens_fcst_complete <- ens_fcst %>% filter(init_year %in% study_period)

  ## Do unit conversion
  ens_fcst_complete <- ens_fcst_complete %>% ensemble_fcst_unit_conversion()

  ## Aggregate over lead times
  group_vars = c("project", "mip", "source_id", "member", "init_year")
  ens_fcst_complete <-
    ens_fcst_complete %>%
    group_by_at(group_vars) %>%
    summarize(across(starts_with(vars), mean))

  ## Complete anomalies over the study period
  compute_anomaly <- function(x) x - mean(x, na.rm = TRUE)
  anomaly_group_vars <- c("source_id", "member")
  ens_fcst_complete <-
    ens_fcst_complete %>%
    group_by_at(anomaly_group_vars) %>%
    mutate(across(starts_with(vars), compute_anomaly)) %>%
    ungroup()

  ## Pivot longer
  ens_fcst_complete <-
    ens_fcst_complete %>%
    pivot_longer(-all_of(group_vars), names_to="variable", values_to="value")
  ens_fcst_complete
}

get_hindcast_data <- function(dataset,
                              study_period,
                              lead_times,
                              all_ids = TRUE,
                              id=NA) {

  ensemble_fcst_raw = open_dataset(dataset)
  if (!all_ids & !is.na(id)) {
    ensemble_fcst_raw <-
      ensemble_fcst_raw %>%
      filter(ID %in% id)
  }
  ensemble_fcst_raw <-
    ensemble_fcst_raw %>%
    ## open_dataset(dataset) %>%
    mutate(lead_time = season_year - init_year) %>%
    filter(lead_time %in% lead_times) %>%
    collect()
  ## Pivot from long to wide format
  ensemble_fcst_raw <- ensemble_fcst_raw %>%
    pivot_wider(names_from=variable, values_from=value)
  ## Unit conversion
  ensemble_fcst_raw <-
    ensemble_fcst_raw %>%
    mutate(across(matches("^nao$"), function(x) x / 100)) %>%
    mutate(across(matches("^ea$"), function(x) x / 100)) %>%
    mutate(across(starts_with("european_precip"), function(x) x * 60 * 60 * 24)) %>%
    mutate(across(starts_with("uk_precip"), function(x) x * 60 * 60 * 24))
    ## mutate(european_precip = european_precip * 60 * 60 * 24) %>%
    ## mutate(uk_precip = uk_precip * 60 * 60 * 24)
  ## ## Correct initialisation years for GFDL data, which appear to be incorrect
  ## gfdl_index = ensemble_fcst_raw$source_id %in% "GFDL-CM2p1"
  ## ensemble_fcst_raw$init_year[gfdl_index] = ensemble_fcst_raw$init_year[gfdl_index] - 1
  ensemble_fcst_raw_complete =
    ensemble_fcst_raw %>%
    filter(init_year %in% study_period)
  ensemble_fcst_raw_complete
}

summarise_discharge_data <- function(x, metadata) {

  ## ## Previous POT analysis
  ## ## Neri et al [https://doi.org/10.1002/joc.5915]:
  ## ## "To avoid double counting the same event, we only consider
  ## ## one event in a window of +/- 5 days + logarithm of the
  ## ## drainage area"
  ## catchment_area = metadata[["catchment_area"]]
  ## window_size = ((5 + log(catchment_area * 0.386102)) * 2) %>% round()
  ## ## Solari et al [https://doi.org/10.1002/2016WR019426]:
  ## ## "As the moving window travels through the series, each
  ## ## time that the data maximum in the window is located at
  ## ## its center, the maximum is regarded as a peak"
  ## df =
  ##   df %>% # Neri et al
  ##   mutate(pot = roll_max(gdf, n=window_size, align="center", fill=NA)) %>%
  ##   mutate(is_peak = Vectorize(isTRUE)(pot == gdf)) %>%
  ##   mutate(peak = ifelse(is_peak, pot, NA))
  ## ## Unsure whether we need to make this unique or not?
  ## peaks = df$pot[df$is_peak] ##%>% unique()
  ## n_years = length(df$year %>% unique())
  ## peaks_sorted = sort(peaks, decreasing = TRUE)
  ## threshold_1 = peaks_sorted[n_years]
  ## threshold_2 = peaks_sorted[n_years * 2]
  ## threshold_3 = peaks_sorted[n_years * 3]
  ## threshold_4 = peaks_sorted[n_years * 4]

  ## df =
  ##   df %>%
  ##   mutate(pot_1 = ifelse(Vectorize(isTRUE)(peak >= threshold_1), 1, 0)) %>%
  ##   mutate(pot_2 = ifelse(Vectorize(isTRUE)(peak >= threshold_2), 1, 0)) %>%
  ##   mutate(pot_3 = ifelse(Vectorize(isTRUE)(peak >= threshold_3), 1, 0)) %>%
  ##   mutate(pot_4 = ifelse(Vectorize(isTRUE)(peak >= threshold_4), 1, 0))

  ## Add climate season label to rows
  x <- x %>%
    mutate(
      clim_season = case_when(
        month %in% c(12, 1, 2, 3) ~ "DJFM",
        month %in% c(4, 5) ~ "AM",
        month %in% c(6, 7, 8, 9) ~ "JJAS",
        month %in% c(10, 11) ~ "ON"
      )
    ) %>%
    mutate(season_year = ifelse(month %in% c(1, 2, 3), year - 1, year))

  ## Quantiles per season
  thresholds <- x %>%
    group_by(ID, clim_season) %>%
    summarize(
      Q_99_threshold = quantile(Q, 0.99, na.rm = TRUE),
      Q_95_threshold = quantile(Q, 0.95, na.rm = TRUE)
    )

  x <- x %>% left_join(thresholds, by=c("ID", "clim_season"))

  ## Summarize to get flood counts
  x <- x %>%
    group_by(ID, clim_season, season_year) %>%
    summarize(
      missing_pct = (sum(is.na(Q)) / n()) * 100,
      Q_max = max(Q, na.rm = TRUE),
      Q_mean = mean(Q, na.rm = TRUE),
      Q_05 = quantile(Q, probs = 0.05, na.rm = TRUE, names = FALSE),
      Q_50 = quantile(Q, probs = 0.50, na.rm = TRUE, names = FALSE),
      Q_90 = quantile(Q, probs = 0.90, na.rm = TRUE, names = FALSE),
      Q_95 = quantile(Q, probs = 0.95, na.rm = TRUE, names = FALSE),
      ## P_sum = sum(cdr, na.rm = TRUE),
      ## POT_1 = sum(pot_1, na.rm = TRUE),
      ## POT_2 = sum(pot_2, na.rm = TRUE),
      ## POT_3 = sum(pot_3, na.rm = TRUE),
      ## POT_4 = sum(pot_4, na.rm = TRUE),
      POT_1 = sum(Q > Q_99_threshold, na.rm = TRUE),
      POT_2 = sum(Q > Q_95_threshold, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(across(Q_max:POT_2, ~ifelse(is.finite(.), ., NA)))
    ## mutate(across(Q_max:POT_4, ~ifelse(is.finite(.), ., NA)))
  x
}

download_nrfa_data <- function(stn_id, metadata) {
  ## TODO tidy up this function, giving user more control over which variables are derived
  meta = metadata %>% filter(id %in% stn_id)
  ## Gauged daily flow [m3 s-1]
  gdf = get_ts(stn_id, "gdf") %>% as_tibble(rownames="time")
  ## ## Catchent daily rainfall [mm]
  ## cdr = try(get_ts(stn_id, "cdr") %>% as_tibble(rownames="time"))
  ## if (inherits(cdr, "try-error"))
  ##   cdr <- gdf %>% rename(cdr = gdf) %>% mutate(cdr = NA)

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
  x <- gdf %>%
    rename(Q = gdf) %>%
    ## left_join(cdr, by="time") %>%
    mutate(ID=stn_id, .after=time) %>%
    mutate(time = as.Date(time))
  ## TODO Filter out years with fewer than 330 days of records
  x <- x %>%
    mutate(year = format(time, "%Y") %>% as.integer) %>%
    mutate(month = format(time, "%m") %>% as.integer)

  x <- summarise_discharge_data(x, meta)
  x
}

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
      out[[nm]][i] = fun(x[start_index:end_index], na.rm=TRUE)
    }
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
    corr_this = try(cor.test(
      fcst[keep_index],
      obs[keep_index],
      method="pearson", alternative="greater"
    ), silent=TRUE)
    if (isTRUE(inherits(corr_this, "try-error"))) {
      corr[i] = NA
    } else {
      corr[i] = unname(corr_this$estimate)
    }
  }
  corr
}

sd_cross_validate <- function(x, leave_out_add=1) {
  ntimes = length(x)
  index = seq(1, ntimes) # for selection
  sdev = rep(0., ntimes)
  for (i in 1:ntimes) {
    if (leave_out_add > 0) {
      leave_out_start = max(c(1, i - leave_out_add))
      leave_out_end = min(c(ntimes, i + leave_out_add))
      leave_out_index = seq(leave_out_start, leave_out_end)
    } else {
      leave_out_index = c()
    }
    keep_index = !(index %in% leave_out_index)
    sd_this = try(sd(x[keep_index], na.rm=TRUE), silent=TRUE)
    if (isTRUE(inherits(sd_this, "try-error"))) {
      sdev[i] = NA
    } else {
      sdev[i] = sd_this
    }
  }
  sdev
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

  ## Join the data frame with ensemble_fcst, which approximately
  ## quadruples its size (values are duplicated when they are
  ## included in a different lag year (i.e. `init_year_lag`))
  ensemble_fcst_lag =
    ensemble_fcst %>%
    filter(variable %in% match_var) %>%
    left_join(init_year_lag, by = "init_year") %>%
    arrange(source_id, member, init_year_lag, init_year)

  ## Join the datasets
  ensemble_fcst_lag =
    ensemble_fcst_lag %>%
    left_join(match_fcst, by = c("variable", "init_year_lag"))

  ## Calculate mean absolute error for the matching variable
  ensemble_fcst_lag =
    ensemble_fcst_lag %>%
    ## filter(variable %in% match_var) %>%
    mutate(error = abs(std - ens_mean_lag_std))

  ## Prepare output data
  nao_matched_ensemble_fcst = ensemble_fcst_lag %>% ungroup()
  nao_matched_ensemble_fcst =
    nao_matched_ensemble_fcst %>%
    dplyr::select(project, source_id, mip, member, init_year, init_year_lag, error) %>%
    rename(init_year_matched = init_year) %>%
    rename(init_year = init_year_lag)
  return(nao_matched_ensemble_fcst)
}

parse_grid_cell <- function(x) {
  grid_lat = str_extract(x, "(N|S)\\d+\\.*\\d*")
  grid_lon = str_extract(x, "(E|W)\\d+\\.*\\d*")
  north = toupper(str_sub(grid_lat, 1, 1)) == "N"
  east = toupper(str_sub(grid_lon, 1, 1)) == "E"
  grid_lat = str_sub(grid_lat, 2, -1) %>% as.numeric()
  grid_lon = str_sub(grid_lon, 2, -1) %>% as.numeric()
  grid_lat = ifelse(north, grid_lat, grid_lat * -1)
  grid_lon = ifelse(east, grid_lon, grid_lon * -1)
  lapply(seq_len(length(x)), FUN=function(i) c(grid_lon[i], grid_lat[i]))
}

select_nearest_grid_cell <- function(lat, lon, grid_coords) {
  ## grid_coords <- c(
  ##   "uk_precip_field_N50.0_W0.0",
  ##   "uk_precip_field_N50.0_W10.0",
  ##   "uk_precip_field_N50.0_W5.0",
  ##   "uk_precip_field_N55.0_W0.0",
  ##   "uk_precip_field_N55.0_W10.0",
  ##   "uk_precip_field_N55.0_W5.0",
  ##   "uk_precip_field_N60.0_W0.0",
  ##   "uk_precip_field_N60.0_W10.0",
  ##   "uk_precip_field_N60.0_W5.0")
  ## lat = 52.6
  ## lon = -2.6
  coords <- parse_grid_cell(grid_coords)
  dists <- sapply(coords, FUN=function(x) geosphere::distHaversine(c(lon, lat), x))
  i <- which.min(dists)
  grid_coords[i]
}

create_ensemble_forecast <- function(ensemble_fcst_error,
                                     ensemble_fcst,
                                     vars,
                                     model_select = NA,
                                     project_select = NA,
                                     full = TRUE,
                                     best_n = FALSE,
                                     worst_n = FALSE,
                                     n_select = 20) {

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
    summarize(across(starts_with(vars), mean, na.rm = TRUE))
    ## summarize(across(all_of(vars), mean, na.rm = TRUE))
  ensemble_fcst
}


make_prediction <- function(newdata,
                            data,
                            ...,
                            quantiles = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.98, 0.99),
                            model_family,
                            obs,
                            compute_skill = FALSE,
                            ref_prediction = NA) {
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
  if (!model_family %in% c("GA", "PO")) {
    msg <- paste0("Model family ", model_family, " currently not supported")
    stop(msg)
  }

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
    tbl <- matrix(NA, nrow=nrow(newdata), ncol=length(quantile_names) + 3) %>%
      as_tibble() %>%
      setNames(c("date", "model", "observations", quantile_names))
    ## tbl =
    ##   rep(NA, length(quantile_names) + 1) %>%
    ##   setNames(c("model", quantile_names)) %>%
    ##   as.list() %>% as_tibble()
    tbl[["date"]] = NA
    tbl[["model"]] = nm
    tbl[["observations"]] = obs
    if (!is.na(model)) {
      mu = predict(
        model,
        newdata = newdata,
        what = "mu",
        type = "response",
        data = data
      )
      if (model_family == "GA") {
        sigma = predict(
          model,
          newdata = newdata,
          what = "sigma",
          type = "response",
          data = data
        )
        shape = 1 / (sigma ** 2)
        scale = mu / shape

      } else {
        sigma = NA
        shape = NA
        scale = NA
      }
      tbl[["mu"]] = mu
      tbl[["sigma"]] = sigma
      tbl[["shape"]] = shape
      tbl[["scale"]] = scale

      if (compute_skill) {
        if (is.na(ref_prediction)) {
          stop("`ref_prediction` cannot be be NA if `compute_skill` is TRUE")
        }
        ## See http://www.gamlss.com/wp-content/uploads/2013/01/gamlss-manual.pdf [A.5.1]
        ## We need scale/shape to implement scoringRules, but the GAMLSS implementation
        ## of the Gamma distn is a reparameterization which sets sigma**2=1/shape and
        ## mu=shape*scale [TODO - check this with Louise]
        ## Compute CRPS & CRPSS
        if (model_family == "GA") {
          tbl[["crps_fcst"]] <- crps_gamma(obs, shape=shape, scale=scale)
          ## Also try computing CRPS with ensemble approach
          ens_fcst <- qGA(seq(0.01, 0.99, by = 0.01), mu, sigma)

        } else {
          tbl[["crps_fcst"]] <- NA
          ens_fcst <- qPO(seq(0.01, 0.99, by = 0.01), mu)
        }
        ## tbl[["crps_fcst"]] = crps_gamma(obs, shape=shape, scale=scale)
        ## ## Also try computing CRPS with ensemble approach
        ## ens_fcst = qGA(seq(0.01, 0.99, by = 0.01), mu, sigma)
        ens_fcst_mat = t(matrix(ens_fcst))
        if (nrow(ens_fcst_mat) != length(obs)) {
          stop()
        }
        tbl[["crps_ens_fcst"]] = EnsCrps(ens_fcst_mat, obs)
        ## Ensemble climatology forecast
        ref_prediction_mat = t(matrix(ref_prediction))
        if (nrow(ref_prediction_mat) != length(obs)) {
          stop()
        }
        tbl[["crps_climat"]] = EnsCrps(ref_prediction_mat, obs)
        tbl[["aic"]] = AIC(models[[nm]])
      }
      ## Predict quantiles using mu and sigma
      for (j in 1:length(quantile_names)) {
        q = quantiles[j]
        qnm = quantile_names[j]
        if (model_family == "GA") {
          tbl[[qnm]] = qGA(q, mu, sigma)
        } else {
          tbl[[qnm]] = qPO(q, mu)
        }
      }

    } else {
      tbl[["mu"]] = NA
      tbl[["sigma"]] = NA
      tbl[["shape"]] = NA
      tbl[["scale"]] = NA
      if (compute_skill) {
        tbl[["crps_fcst"]] = NA
        tbl[["crps_ens_fcst"]] = NA
        tbl[["crps_climat"]] = NA
        tbl[["aic"]] = NA
      }
      for (qnm in quantile_names) {
        tbl[[qnm]] = NA
      }
    }
    ## Merge with ID columns from new data
    tbl = cbind(
      newdata %>% dplyr::select(any_of(id_cols)),
      tbl
    )
    computed_quantiles[[i]] = tbl
  }
  ## Combine predictions from all models
  computed_quantiles =
    do.call("rbind", computed_quantiles) %>%
    as_tibble()
  return(computed_quantiles)
}

fit_models_cv <- function(x,
                          config_section,
                          lead_time,
                          ...) {

  catchment_prediction_list = list()
  catchment_simulation_list = list()

  n_fold <- nrow(x)
  for (p in 1:n_fold) {
    ## The buffer is needed to prevent contamination between the training set and test set
    buffer <- length(lead_time) - 1
    idx = seq_len(n_fold)
    test_idx <- p
    remove_idx <- seq(p - buffer, p + buffer)
    remove_idx <- remove_idx[remove_idx > 0]
    train_idx = idx[!idx %in% remove_idx]
    train_data = x[train_idx,]
    test_data = x[test_idx,]
    year <- test_data$year # Should be length 1

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

    ## ############### ##
    ## Make prediction ##
    ## ############### ##

    ## Assume climatology is the mean value of Q
    ens_climatology_prediction <- train_data[["Q"]]
    climatology_prediction <- ens_climatology_prediction %>% mean(na.rm = TRUE)
    obs = test_data[["Q"]]
    prediction = do.call(
      "make_prediction",
      c(list(newdata = test_data, data = train_data),
        models,
        list(model_family = config_section$model_family,
             obs = obs,
             compute_skill = TRUE,
             ref_prediction = ens_climatology_prediction))
    )

    ## Create output data frame
    ## N.B. Q50 is the median of the probabilistic forecast, not the streamflow quantile
    obs_column_name = paste0(config_section$predictand, '_obs')
    exp_column_name = paste0(config_section$predictand, '_exp')
    climat_column_name = paste0(config_section$predictand, '_climat')
    prediction =
      prediction %>%
      mutate(
        !!obs_column_name := observations,
        !!exp_column_name := Q50,
        !!climat_column_name := climatology_prediction
      ) #%>%
      #mutate(date = NA, year = year)

    ## Add data to list
    pred_idx = length(catchment_prediction_list) + 1
    catchment_prediction_list[[pred_idx]] = prediction

    ## ## ############### ##
    ## ## Make simulation ##
    ## ## ############### ##

    ## obs = train_data[["Q"]]
    ## simulation = do.call(
    ##   "make_prediction",
    ##   c(list(newdata = train_data, data = train_data),
    ##     models,
    ##     list(model_family = config_section$model_family))

    ## obs_column_name = paste0(config_section$predictand, '_obs')
    ## exp_column_name = paste0(config_section$predictand, '_exp')
    ## climat_column_name = paste0(config_section$predictand, '_climat')
    ## simulation =
    ##   simulation %>%
    ##   mutate(
    ##     !!obs_column_name := obs,
    ##     !!exp_column_name := Q50,
    ##     !!climat_column_name := climatology_prediction
    ##   ) %>%
    ##   mutate(date = NA, year = year) %>%
    ## mutate(test_year = test_year)

    ## sim_idx = length(simulation_prediction_list) + 1
    ## catchment_simulation_list[[sim_idx]] = simulation

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
  catchment_simulation_list = list()

  buffer <- length(lead_time)
  max_training_period_end = test_period_end - buffer #window
  while (training_period_end <= max_training_period_end) {
    training_years <- seq(training_period_start, training_period_end)
    test_year <- max(training_years) + buffer
    train_data <- x %>% filter(year %in% training_years)
    test_data <- x %>% filter(year %in% test_year)

    ## Only continue if there is sufficient training/testing data
    ## if (nrow(train_data) == 0 | nrow(test_data) == 0) {
    if (nrow(train_data) < 10 | nrow(test_data) == 0) {
      training_period_end <- training_period_end + 1
      next
    }

    ## START SECTION - common with fit_models_cv
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
    ## ## Assume climatology is the mean value of Q
    ## ens_climatology_prediction <- train_data[["Q"]]
    ## climatology_prediction <- ens_climatology_prediction %>% mean(na.rm = TRUE)
    ## obs = test_data[["Q"]]
    ## prediction = do.call(
    ##   "make_prediction",
    ##   c(list(newdata = test_data, data = train_data),
    ##     models,
    ##     list(model_family = config_section$model_family,
    ##          obs = obs,
    ##          ref_prediction = ens_climatology_prediction))
    ## )
    ## ## Create output data frame
    ## obs_column_name = paste0(config_section$predictand, '_obs')
    ## exp_column_name = paste0(config_section$predictand, '_exp')
    ## climat_column_name = paste0(config_section$predictand, '_climat')
    ## ## N.B. Q50 is the median of the probabilistic forecast, not the streamflow quantile
    ## prediction =
    ##   prediction %>%
    ##   mutate(
    ##     !!obs_column_name := obs,
    ##     !!exp_column_name := Q50,
    ##     !!climat_column_name := climatology_prediction
    ##   ) %>%
    ##   mutate(
    ##     date = NA,
    ##     year = year
    ##     ## model = paste0("GAMLSS_", model)
    ##   ) #%>%
    ##   ## dplyr::select(all_of(c(obs_column_name, exp_column_name)), year, model)

    ## pred_idx = length(catchment_prediction_list) + 1
    ## catchment_prediction_list[[pred_idx]] = prediction

    ## ############### ##
    ## Make prediction ##
    ## ############### ##

    ## Assume climatology is the mean value of Q
    ens_climatology_prediction <- train_data[["Q"]]
    climatology_prediction <- ens_climatology_prediction %>% mean(na.rm = TRUE)
    obs = test_data[["Q"]]
    prediction = do.call(
      "make_prediction",
      c(list(newdata = test_data, data = train_data),
        models,
        list(model_family = config_section$model_family,
             obs = obs,
             compute_skill = TRUE,
             ref_prediction = ens_climatology_prediction))
    )

    ## Create output data frame
    ## N.B. Q50 is the median of the probabilistic forecast, not the streamflow quantile
    obs_column_name = paste0(config_section$predictand, '_obs')
    exp_column_name = paste0(config_section$predictand, '_exp')
    climat_column_name = paste0(config_section$predictand, '_climat')
    prediction =
      prediction %>%
      mutate(
        !!obs_column_name := observations,
        !!exp_column_name := Q50,
        !!climat_column_name := climatology_prediction
      ) #%>%
      #mutate(date = NA, year = year)

    ## Add data to list
    pred_idx = length(catchment_prediction_list) + 1
    catchment_prediction_list[[pred_idx]] = prediction

    ## ############### ##
    ## Make simulation ##
    ## ############### ##

    obs = train_data[["Q"]]
    simulation = do.call(
      "make_prediction",
      c(list(newdata = train_data, data = train_data),
        models,
        list(model_family = config_section$model_family, obs = obs))
    )

    obs_column_name = paste0(config_section$predictand, '_obs')
    exp_column_name = paste0(config_section$predictand, '_exp')
    simulation =
      simulation %>%
      mutate(
        !!obs_column_name := observations,
        !!exp_column_name := Q50,
      ) %>%
      mutate(test_year = test_year)

    sim_idx = length(catchment_simulation_list) + 1
    catchment_simulation_list[[sim_idx]] = simulation

    ## Update training period end for next iteration
    training_period_end <- training_period_end + 1
  }
  if (length(catchment_prediction_list) == 0) {
    catchment_prediction <- NULL
    catchment_simulation <- NULL
  } else {
    catchment_prediction = do.call("rbind", catchment_prediction_list)
    catchment_simulation = do.call("rbind", catchment_simulation_list)
  }
  list(prediction = catchment_prediction, simulation = catchment_simulation)
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
  if (!model_family %in% c("GA", "PO")) {
    stop(paste0("Model family '", model_family, "' not yet supported"))
  }
  fitted_models = list()
  for (i in 1:length(model_names)) {
    model_nm = model_names[i]
    if (model_family == "GA") {
      model <- try(
        gamlss(
          formula = formulas[[i]],
          sigma.formula = sigma_formulas[[i]],
          family = GA,
          trace = FALSE,
          data = data,
          ...
        ), silent=TRUE
      )
    } else {
      model <- try(
        gamlss(
          formula = formulas[[i]],
          family = PO,
          trace = FALSE,
          data = data,
          ...
        ), silent=TRUE
      )
    }
    if (inherits(model, "try-error")) {
      model = NA
    }
    fitted_models[[model_nm]] = model
  }
  fitted_models
}

get_aic <- function(model_list) {
  aic = sapply(
    model_list,
    FUN=function(x) ifelse(is.null(x), NA, AIC(x))
  ) %>% as_tibble_row()
  aic
}

check_residuals <- function (x) {
  if (!is.gamlss(x))
    stop(paste("This is not an gamlss object", "\n", ""))
  if (is.null(x$residuals)) #
    stop(paste("There are no quantile residuals in the object"))
  residx <- resid(x) # get the residuals
  w <- x$weights
  qq <- as.data.frame(qqnorm(residx, plot = FALSE))
  Filliben <- cor(qq$y,qq$x)
  m.1 <- mean(residx)
  m.2 <- var(residx) # cov.wt(mr,w)$cov
  n.obs <- sum(w)
  m.3 <- sum((residx-m.1)**3)/n.obs
  m.4 <- sum((residx-m.1)**4)/n.obs
  b.1 <- m.3^2/m.2^3
  sqrtb.1 <- sign(m.3)*sqrt(abs(b.1))
  b.2 <- m.4/m.2^2
  return(list(mean=m.1, variance=m.2, skewness=sqrtb.1, kurtosis=b.2, filliben=Filliben, nobs=n.obs))
}

get_residual_checks <- function(model_list) {
  rows <- list()
  for (i in 1:length(model_list)) {
    nm <- names(model_list)[i]
    rows[[i]] <- check_residuals(model_list[[i]]) %>% as_tibble_row() %>% mutate(model = nm)
  }
  return(do.call("rbind", rows))
}

myfun <- function(x) {
  ## Select the best model based on AIC value
  x_best =
    x %>%
    group_by(ID, subset, period) %>%
    filter(crps_ens_fcst==min(crps_ens_fcst))
    ## filter(aic == min(aic))
    ## filter(!!sym(skill_measure) == max(!!sym(skill_measure)))
  skill =
    gauge_stns %>%
    left_join(x_best) %>%
    mutate(skill = !!sym(skill_measure))
  skill
}

load_model_predictions <- function(config, experiment, aggregation_period) {
  predictions = open_dataset(
    file.path(outputroot, "analysis", experiment, "gamlss", aggregation_period, "prediction")
  ) %>%
    collect() #%>%
    ## filter(subset %in% c("full", "best_n")) #%>%
    ## mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))
  model_levels <- unique(predictions$model)
  model_labels <- model_levels %>% gsub("_", "", .)
  predictions <-
    predictions %>%
    mutate(model = factor(model, levels = model_levels, labels = model_labels))
  predictions
}

load_model_simulations <- function(config, experiment, aggregation_period) {
  simulations <- open_dataset(
    file.path(outputroot, "analysis", experiment, "gamlss", aggregation_period, "simulation")
  ) %>%
    collect() %>%
    filter(subset %in% c("full", "best_n")) #%>%
    ## mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))
  model_levels <- unique(simulations$model)
  model_labels <- model_levels %>% gsub("_", "", .)
  simulations <-
    simulations %>%
    mutate(model = factor(model, levels = model_levels, labels = model_labels))
  simulations
}

load_model_fit <- function(config, experiment, aggregation_period) {
  fit <- open_dataset(
    file.path(outputroot, "analysis", experiment, "gamlss", aggregation_period, "fit")
  ) %>% collect()
  fit <- fit %>% mutate(period = aggregation_period)
  model_levels <- unique(fit$model)
  model_labels <- model_levels %>% gsub("_", "", .)
  fit <- fit %>% mutate(model = factor(model, levels = model_levels, labels = model_labels))
  fit
}

load_skill_scores <- function(config, experiment, aggregation_period) {
  ds <- open_dataset(
    file.path(outputroot, "analysis", experiment, "gamlss", aggregation_period, "prediction")
  ) %>% collect()
  skill <- ds %>%
    group_by(ID, model, subset) %>%
    summarize(
      crps_fcst = mean(crps_fcst),
      crps_ens_fcst = mean(crps_ens_fcst),
      crps_climat = mean(crps_climat),
      aic=mean(aic)
    ) %>%
    ## mutate(crpss = 1 - (crps_fcst / crps_climat)) %>%
    mutate(crpss = 1 - (crps_ens_fcst / crps_climat)) %>% # TODO check
    mutate(period = aggregation_period)

  model_levels <- unique(skill$model)
  model_labels <- model_levels %>% gsub("_", "", .)
  skill <- skill %>% mutate(model = factor(model, levels = model_levels, labels = model_labels))

  return(skill)
}


## mean_square_error_skill_score <- function(obs, exp) {
##   ## MSSS
##   mse = mean((exp - obs) ^ 2)
##   mse_ref = mean((mean(obs) - obs) ^ 2)
##   msss = 1 - (mse / mse_ref)
##   ## ACC
##   acc = cor(obs, exp, method = "pearson")
##   ## correlation
##   r = cor(obs, exp, method = "pearson")
##   ## potential skill [= coefficient of determination]
##   ps = r ^ 2
##   ## slope reliability
##   srel = (r - (sd(exp) / sd(obs))) ^ 2
##   ## standardized mean error
##   sme = ((mean(exp) - mean(obs)) / sd(obs)) ^ 2
##   ## msss = ps - srel - sme
##   list(msss = msss, ps = ps, srel = srel, sme = sme, acc = acc)
## }

## NOT USED:
##
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
