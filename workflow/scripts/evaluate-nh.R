#!/usr/bin/env Rscript

## Author : Simon Moulds
## Date   : Jan 2022

library(tidyverse)
library(lubridate)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(scales)
library(arrow)
library(sf)
library(yaml)

options(dplyr.summarise.inform = FALSE)

## ## Extract configuration info
## if (sys.nframe() == 0L) {
##   args = commandArgs(trailingOnly=TRUE)
##   config = read_yaml(args[1])
##   outputroot = args[2]
##   ## args = commandArgs()
##   ## m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
##   ## cwd <- dirname(regmatches(args, m))
## }

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

## ## For testing:
## outputroot <- "results"
## nh_inputdir <- file.path(outputroot, "analysis", "yr2to9_lag", "nh-output", "time_series")
## plot_outputdir <- file.path(outputroot, "analysis", "yr2to9_lag", "nh-output", "plots")
## gamlss_inputdir <- file.path(outputroot, "analysis", "yr2to9_lag", "gamlss-output", "prediction")
## metadata <- read_parquet(file.path(outputroot, "nrfa-metadata.parquet"))

mean_square_error_skill_score = function(obs, exp, na.rm = TRUE) {
  if (na.rm) {
    idx <- is.na(obs) | is.na(exp)
    obs <- obs[!idx]
    exp <- exp[!idx]
  }
  ## mse = mean((exp - obs) ^ 2)
  ## mse_ref = mean((mean(obs) - obs) ^ 2)
  ## msss = 1 - (mse / mse_ref)
  acc = cor(obs, exp, method = "pearson")
  ## correlation
  r = cor(obs, exp, method = "pearson")
  ## potential skill [= coefficient of determination]
  ps = r ^ 2
  ## slope reliability
  srel = (r - (sd(exp) / sd(obs))) ^ 2
  ## standardized mean error
  sme = ((mean(exp) - mean(obs)) / sd(obs)) ^ 2
  msss = ps - srel - sme
  list(msss = msss, ps = ps, srel = srel, sme = sme, acc = acc)
}

plot_outputdir <- 'results/test_plots'
if (dir.exists(plot_outputdir)) {
  unlink(plot_outputdir, recursive = TRUE)
}
dir.create(plot_outputdir, recursive = TRUE)

## ## fs <- list.files('results', pattern = "test_[0-9]{4}.csv", full.names = TRUE)
## ## y <- read_csv(fs, show_col_types = FALSE) %>%
## y1 <- read_csv('results/test_xgboost.csv', show_col_types = FALSE) %>%
##   rename(ID = basin) %>%
##   mutate(year = year(date)) %>%
##   arrange(ID, year) %>%
##   mutate(model = "xgboost")

## y2 <- read_csv('results/test_tabnet.csv', show_col_types = FALSE) %>%
##   rename(ID = basin) %>%
##   mutate(year = year(date)) %>%
##   arrange(ID, year) %>%
##   mutate(model = "tabnet")

## y <- rbind(y1, y2) %>% mutate(model = factor(model, levels = c("xgboost", "tabnet"), labels = c("XGBoost", "TabNet")))

y <- read_csv('results/test.csv', show_col_types = FALSE) %>%
  rename(ID = basin) %>%
  mutate(year = year(date)) %>%
  arrange(ID, year) %>%
  mutate(model = "tabnet") %>%
  mutate(model = factor(model, levels = c("tabnet"), labels = c("TabNet")))

## y <- read_csv('results/test.csv', show_col_types = FALSE) %>% rename(ID = basin) %>% mutate(year = year(date))
## y <- open_dataset(gamlss_inputdir) %>% collect()
station_ids <- y$ID %>% unique()
n_stations <- length(station_ids)

msss_list <- list()
acc_list <- list()
for (i in 1:n_stations) {
  stn <- station_ids[i]
  yy <- y %>% filter(ID %in% stn)
  yy_xgb <- yy %>% filter(model %in% "XGBoost")
  yy_tab <- yy %>% filter(model %in% "TabNet")
  skill_xgb <- mean_square_error_skill_score(yy_xgb$obs, yy_xgb$exp)
  skill_tab <- mean_square_error_skill_score(yy_tab$obs, yy_tab$exp)
  msss_list[[i]] <- tibble(id = stn, xgb = skill_xgb$msss, tab = skill_tab$msss)
  acc_list[[i]] <- tibble(id = stn, xgb = skill_xgb$acc, tab = skill_tab$acc)

  p <- ggplot() +
    theme_bw() +
    geom_line(
      aes(y=exp, x=year, colour = model), data=yy, # colour = "magenta"
    ) +
    ## scale_fill_manual(values = cbbPalette) +
    ## scale_color_manual(values = cbbPalette) +
    ## scale_fill_discrete(values = cbbPalette) +
    ## scale_colour_discrete(values = cbbPalette) +
    ## facet_wrap(. ~ ID, ncol = 1) + #, labeller = label_parsed) +
    ylab(expression(Streamflow~(mm~d^{-1}))) +
    xlab("") +
    ## scale_x_continuous(breaks = pretty_breaks()) +
    ## ## scale_y_continuous(expression(Streamflow~(m^{3}~s^{-1})), breaks = pretty_breaks()) +
    ## ## N.B. use alpha to create another legend
    ## ## https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
    ## ## geom_point(
    geom_line(
      aes(y=obs, x=year), #, alpha="Observed"),
      color = "black",
      data=yy,
      size = 1 #0.2
    ) +
    ## scale_alpha_manual(name=NULL, values=1, breaks="Observed") +
    ## ggtitle(sprintf("ID = %d", stn_id)) +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          legend.title = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = strip_label_size),
          ## legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_label_size),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = axis_label_size_small),
          axis.text.x = element_text(size = axis_label_size_small))
  ggsave(file.path(plot_outputdir, paste0("ts_plot_", stn, ".png")), width = 5, height = 5, units = "in")
}

## names(msss_list) <- station_ids
## names(acc_list) <- station_ids
msss <- do.call("rbind", msss_list)
acc <- do.call("rbind", acc_list)

## msss <- as.numeric(msss_list)
## acc <- as.numeric(acc_list)

## msss_list <- list()
## for (i in 1:n_stations) {
##   stn <- station_ids[i]
##   catchment_area <- metadata %>% filter(id %in% stn) %>% `$`(catchment_area)

##   ## Read GAMLSS model prediction
##   yy <- y %>% filter(ID %in% stn) %>% rename(gamlss_exp = exp) %>% dplyr::select(year, obs, gamlss_exp)
##   yy$obs <- yy$obs * 86400 / catchment_area * 1000 / 1000 / 1000
##   yy$gamlss_exp <- yy$gamlss_exp * 86400 / catchment_area * 1000 / 1000 / 1000

##   ## Read LSTM prediction
##   xx <- read_csv(file.path(nh_inputdir, paste0(stn, ".csv")), show_col_types = FALSE)
##   xx <- xx %>%
##     mutate(year = lubridate::year(date)) %>%
##     rename(obs = Q95_obs, nh_exp = Q95_sim) %>%
##     dplyr::select(year, obs, nh_exp, -date, -time_step) %>%
##     mutate(year = year - 1) # FIXME
##     ## gather(-season_year, key = "key", value = "value")

##   ## Compute MSSS
##   idx <- is.na(yy$gamlss_exp) | is.na(xx$nh_exp)
##   gamlss_msss <- mean_square_error_skill_score(yy$obs[!idx], yy$gamlss_exp[!idx])$msss
##   lstm_msss <- mean_square_error_skill_score(xx$obs[!idx], xx$nh_exp[!idx])$msss
##   msss_list[[i]] <- tibble(ID = stn, GAMLSS=gamlss_msss, LSTM=lstm_msss)

##   ## Make plot
##   xx <- xx %>%
##     left_join(yy %>% dplyr::select(-obs), by = c("year")) %>%
##     gather(-year, key = "key", value = "value") %>%
##     filter(!key %in% "obs")

##   xx <- xx %>% mutate(key = factor(key, levels = c("nh_exp", "gamlss_exp"), labels = c("LSTM", "GAMLSS")))
##   p <- ggplot() +
##     theme_bw() +
##     geom_line(
##       aes(y=value, x=year, colour = key), data=xx
##     ) +
##     ## scale_fill_manual(values = cbbPalette) +
##     ## scale_color_manual(values = cbbPalette) +
##     ## scale_fill_discrete(values = cbbPalette) +
##     ## scale_colour_discrete(values = cbbPalette) +
##     ## facet_wrap(. ~ ID, ncol = 1) + #, labeller = label_parsed) +
##     ylab(expression(Streamflow~(mm~d^{-1}))) +
##     xlab("") +
##     ## scale_x_continuous(breaks = pretty_breaks()) +
##     ## ## scale_y_continuous(expression(Streamflow~(m^{3}~s^{-1})), breaks = pretty_breaks()) +
##     ## ## N.B. use alpha to create another legend
##     ## ## https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
##     ## ## geom_point(
##     geom_line(
##       aes(y=obs, x=year), #, alpha="Observed"),
##       color = "black",
##       data=yy,
##       size = 1 #0.2
##     ) +
##     ## scale_alpha_manual(name=NULL, values=1, breaks="Observed") +
##     ## ggtitle(sprintf("ID = %d", stn_id)) +
##     theme(legend.position = "bottom",
##           legend.direction = "vertical",
##           legend.title = element_blank(),
##           strip.background = element_blank(),
##           panel.grid = element_blank(),
##           strip.text = element_text(size = strip_label_size),
##           ## legend.title = element_text(size = legend_title_size),
##           legend.text = element_text(size = legend_label_size),
##           axis.title.y = element_blank(),
##           axis.text.y = element_text(size = axis_label_size_small),
##           axis.text.x = element_text(size = axis_label_size_small))
##   ggsave(file.path(plot_outputdir, paste0("ts_plot_", stn, ".png")), width = 5, height = 5, units = "in")
## }

## msss <- do.call("rbind", msss_list)

## length(which(msss$LSTM > msss$GAMLSS))
## length(which(msss$LSTM < 0))
## msss %>% arrange(desc(GAMLSS))
## msss %>% arrange(desc(LSTM))
## msss %>% arrange(LSTM)

## msss <- msss %>% filter(LSTM > 0)
## plot(msss$GAMLSS, msss$LSTM, xlim = c(0, 1), ylim = c(0, 1))
## boxplot(msss$LSTM)
## fs <- list.files("results/analysis/yr2to9_lag/nh-output/time_series", full.names = T)
## data_list <- list()
## for (i in 1:length(fs)) {
##   data_list[[i]] = read_csv(fs[i], show_col_types = FALSE)
## }
## x <- do.call("rbind", data_list)
## x <- x %>% arrange(ID)

## station_ids <- x$ID %>% unique()
## n_stations <- length(station_ids)

## msss_list <- list()
## for (i in 1:n_stations) {
##   stn <- station_ids[i]
##   catchment_area <- metadata %>% filter(id %in% stn) %>% `$`(catchment_area)

##   ## Read GAMLSS model prediction
##   yy <- y %>% filter(ID %in% stn) %>% rename(gamlss_exp = exp) %>% dplyr::select(year, obs, gamlss_exp)
##   yy$obs <- yy$obs * 86400 / catchment_area * 1000 / 1000 / 1000
##   yy$gamlss_exp <- yy$gamlss_exp * 86400 / catchment_area * 1000 / 1000 / 1000

##   ## Read LSTM prediction
##   xx <- x %>% filter(ID %in% stn) %>% arrange(date)
##   ## xx <- read_csv(file.path(nh_inputdir, paste0(stn, ".csv")), show_col_types = FALSE)
##   xx <- xx %>%
##     mutate(year = lubridate::year(date)) %>%
##     rename(obs = Q95_obs, nh_exp = Q95_sim) %>%
##     dplyr::select(year, obs, nh_exp, -date, -time_step) %>%
##     mutate(year = year - 1) # FIXME
##     ## gather(-season_year, key = "key", value = "value")

##   yy <- yy %>% left_join(xx %>% dplyr::select(-obs), by = "year")

##   ## Compute MSSS
##   idx <- is.na(yy$gamlss_exp) | is.na(yy$nh_exp)
##   gamlss_skill <- mean_square_error_skill_score(yy$obs[!idx], yy$gamlss_exp[!idx]) %>%
##     as_tibble() %>%
##     mutate(ID = stn, model = "GAMLSS", .before = msss)
##   lstm_skill <- mean_square_error_skill_score(yy$obs[!idx], yy$nh_exp[!idx]) %>%
##     as_tibble() %>%
##     mutate(ID = stn, model = "LSTM", .before = msss)
##   skill <- rbind(gamlss_skill, lstm_skill)
##   msss_list[[i]] <- skill #tibble(ID = stn, GAMLSS=gamlss_msss, LSTM=lstm_msss)

##   ## Make plot
##   ## xx <- xx %>%
##   ##   left_join(yy %>% dplyr::select(-obs), by = c("year")) %>%
##   plot_data <- yy %>%
##     dplyr::select(-obs) %>%
##     gather(-year, key = "key", value = "value") %>%
##     filter(!key %in% "obs")

##   plot_data <- plot_data %>%
##     mutate(key = factor(key, levels = c("nh_exp", "gamlss_exp"), labels = c("LSTM", "GAMLSS")))

##   p <- ggplot() +
##     theme_bw() +
##     geom_line(
##       aes(y=value, x=year, colour = key), data=plot_data
##     ) +
##     ## scale_fill_manual(values = cbbPalette) +
##     ## scale_color_manual(values = cbbPalette) +
##     ## scale_fill_discrete(values = cbbPalette) +
##     ## scale_colour_discrete(values = cbbPalette) +
##     ## facet_wrap(. ~ ID, ncol = 1) + #, labeller = label_parsed) +
##     ylab(expression(Streamflow~(mm~d^{-1}))) +
##     xlab("") +
##     ## scale_x_continuous(breaks = pretty_breaks()) +
##     ## ## scale_y_continuous(expression(Streamflow~(m^{3}~s^{-1})), breaks = pretty_breaks()) +
##     ## ## N.B. use alpha to create another legend
##     ## ## https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
##     ## ## geom_point(
##     geom_line(
##       aes(y=obs, x=year), #, alpha="Observed"),
##       color = "black",
##       data=yy,
##       size = 1 #0.2
##     ) +
##     ## scale_alpha_manual(name=NULL, values=1, breaks="Observed") +
##     ## ggtitle(sprintf("ID = %d", stn_id)) +
##     theme(legend.position = "bottom",
##           legend.direction = "vertical",
##           legend.title = element_blank(),
##           strip.background = element_blank(),
##           panel.grid = element_blank(),
##           strip.text = element_text(size = strip_label_size),
##           ## legend.title = element_text(size = legend_title_size),
##           legend.text = element_text(size = legend_label_size),
##           axis.title.y = element_blank(),
##           axis.text.y = element_text(size = axis_label_size_small),
##           axis.text.x = element_text(size = axis_label_size_small))
##   ggsave(file.path(plot_outputdir, paste0("ts_plot_", stn, ".png")), width = 5, height = 5, units = "in")
## }

## msss <- do.call("rbind", msss_list)
## msss_lstm <- msss %>% filter(model %in% "LSTM")
## msss_lstm %>% arrange(desc(msss))
## msss_gamlss <- msss %>% filter(model %in% "GAMLSS")
## ## length(which(msss$LSTM > msss$GAMLSS))
## length(which(msss_lstm$msss > 0))
## ## msss %>% arrange(desc(GAMLSS))
## ## msss %>% arrange(LSTM)

## ## msss_lstm_0 <- msss %>% filter(model %in% "LSTM" & msss > 0) %>% arrange(desc(msss))

## ## For spatial plots:
## library(sf)
## uk_boundary =
##   st_read("../data-raw/CNTR_RG_01M_2020_4326.shp") %>%
##   filter(CNTR_NAME %in% "United Kingdom") %>%
##   st_transform(crs = 27700)

## europe_boundary =
##   st_read("../data-raw/CNTR_RG_01M_2020_4326.shp") %>%
##   filter(!CNTR_NAME %in% "United Kingdom") %>%
##   st_transform(crs = 27700)

## library(rnrfa)
## gauge_stns =
##   catalogue() %>%
##   rename(ID = id, area = "catchment-area") %>%
##   filter(ID %in% station_ids) %>%
##   dplyr::select(ID, name, area, latitude, longitude) %>%
##   st_as_sf(coords=c("longitude", "latitude")) %>%
##   st_set_crs(4326) %>%
##   st_transform(27700)

## rdbu_pal = RColorBrewer::brewer.pal(9, "RdBu")
## p = ggplot() +
##   geom_sf(
##     data = europe_boundary,
##     color=NA,
##     fill="lightgrey"
##   ) +
##   geom_sf(
##     data = uk_boundary,
##     lwd = 0.25
##   ) +
##   geom_sf(
##     data = gauge_stns %>% left_join(msss_lstm, by = "ID"),
##     aes(fill = msss),
##     shape = 21,
##     size = 1.5,
##     lwd = 0.1,
##     alpha = 0.8
##   ) +
##   ## facet_wrap(. ~ period, ncol = 1) +
##   coord_sf(
##     xlim = c(-8, 2),
##     ylim = c(50, 59),
##     default_crs = st_crs(4326)
##   ) +
##   ## scale_shape_manual(values = c(21, 24, 22)) +
##   scale_fill_stepsn(
##     colours = rdbu_pal,
##     ## breaks = seq(-0.8, 0.8, 0.2),
##     ## values = scales::rescale(c(-0.8, 0, 0.8)),
##     ## limits = c(-0.3, 0.9)
##     breaks = seq(-0.2, 0.8, 0.2),
##     values = scales::rescale(c(-0.2, 0, 0.8)),
##     limits = c(-0.1, 0.9)
##   ) +
##   theme_bw() +
##   theme(
##     strip.background = element_blank(),
##     ## legend.position = "bottom",
##     ## legend.box = "vertical",
##     ## legend.justification = "left",
##     ## legend.box.just = "left",
##     legend.title = element_text(size = legend_title_size),
##     legend.text = element_text(size = legend_label_size),
##     strip.text = element_blank(),
##     panel.grid.major = element_line(size = 0.25),
##     axis.text = element_text(size = axis_label_size_small)
##   ) +
##   guides(
##     shape = guide_legend(
##       title = "Model",
##       title.position = "top",
##       order = 1
##     ),
##     fill = guide_colorbar(
##       title="MSSS",
##       title.position="top",
##       frame.colour = "black",
##       ticks.colour = "black",
##       frame.linewidth = 0.25,
##       ticks.linewidth = 0.25,
##       barwidth = 0.75,
##       barheight = 10,
##       order = 2
##     )
##   )
