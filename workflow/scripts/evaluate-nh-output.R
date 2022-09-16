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

source("workflow/scripts/external/R/utils.R")

## ## Extract configuration info
## if (sys.nframe() == 0L) {
##   args = commandArgs(trailingOnly=TRUE)
##   config = read_yaml(args[1])
##   outputroot = args[2]
##   ## args = commandArgs()
##   ## m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
##   ## cwd <- dirname(regmatches(args, m))
## }

## For testing:
outputroot <- "results/exp2/analysis/hindcast"
## nh_inputdir <- file.path(outputroot, "analysis", "yr2to9_lag", "nh-output", "time_series")
## plot_outputdir <- file.path(outputroot, "analysis", "yr2to9_lag", "nh-output", "plots")
gamlss_datadir <- file.path(outputroot, "gamlss", "yr2", "prediction")
lstm_datadir <- file.path(outputroot, "lstm", "yr2", "prediction")
xgboost_datadir <- file.path(outputroot, "xgboost", "yr2", "prediction")
tabnet_datadir <- file.path(outputroot, "tabnet", "yr2", "prediction")
## metadata <- read_parquet(file.path(outputroot, "nrfa-metadata.parquet"))

## if (dir.exists(plot_outputdir)) {
##   unlink(plot_outputdir, recursive = TRUE)
## }
## dir.create(plot_outputdir, showWarnings = FALSE)

## Join data frames and reshape
## FIXME - consistent output format
gamlss_output <- open_dataset(gamlss_datadir) %>% collect()
## FIXME - consistent column names
lstm_output <- open_dataset(lstm_datadir) %>% collect() %>% rename(Q95_exp = Q95_sim)
xgboost_output <- open_dataset(xgboost_datadir) %>% collect() %>% mutate(model = "XGBoost")
tabnet_output <- open_dataset(tabnet_datadir) %>% collect() %>% mutate(model = "TabNet")

output <- rbind(xgboost_output, tabnet_output) %>% mutate(date = as.Date(date)) %>% arrange(ID, date, model)

compute_skill_scores <- function(x, ...) {
  station_ids = unique(x$ID) %>% sort()
  skill_scores = list()
  pb <- txtProgressBar(min = 0, max = length(station_ids), initial = 0)
  for (i in 1:length(station_ids)) {
    id <- station_ids[i]
    xx <- x %>% filter(ID %in% id) %>% arrange(year)
    skill <- mean_square_error_skill_score(xx$Q95_obs, xx$Q95_exp) %>% as_tibble()
    skill <- skill %>% mutate(ID = id, .before = names(.)[1])
    skill_scores[[i]] <- skill
    setTxtProgressBar(pb, i)
  }
  close(pb)
  skill_scores <- do.call("rbind", skill_scores)
}

## FIXME - important to have GAMLSS as single-site comparison
## gamlss_skill_scores
## lstm_skill_scores <- compute_skill_scores(lstm_output)
xgboost_skill_scores <- compute_skill_scores(xgboost_output) %>% mutate(model = "XGBoost", .before = ID)
tabnet_skill_scores <- compute_skill_scores(tabnet_output) %>% mutate(model = "TabNet", .before = ID)
skill_scores <- rbind(xgboost_skill_scores, tabnet_skill_scores)

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

## Spatial data for plotting
uk_boundary =
  st_read("../data-raw/CNTR_RG_01M_2020_4326.shp") %>%
  filter(CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

europe_boundary =
  st_read("../data-raw/CNTR_RG_01M_2020_4326.shp") %>%
  filter(!CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

library(rnrfa)
gauge_stns =
  catalogue() %>%
  rename(ID = id, area = "catchment-area") %>%
  filter(ID %in% station_ids) %>%
  dplyr::select(ID, name, area, latitude, longitude) %>%
  st_as_sf(coords=c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_transform(27700)

plotfun1 <- function(x, ...) {
  ## Boxplot comparison
  p = ggplot(x, aes(x = model, y=acc)) +
    ## stat_boxplot(coef = NULL) +
    geom_boxplot(
      lwd = 0.25,
      outlier.size = 0.25
    ) +
    scale_x_discrete(name = "") +
    ## scale_y_continuous(
    ##   name=legend_title,
    ##   ## breaks=seq(-0.2, 1, by=0.2),
    ##   ## limits=c(-0.25, 1.1)
    ##   breaks=seq(-1, 1, by=0.2),
    ##   limits=c(-1, 1)
    ## ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = axis_title_size_small),
      axis.text = element_text(size = axis_label_size_small)
    )
  p
}
p1 <- plotfun1(skill_scores)

plotfun2 <- function(x, id, ...) {
  xx <- x %>% filter(ID %in% id)
  models <- unique(xx$model)
  yy <- xx %>% filter(model %in% models[1])
  ## Time series plot
  p <- ggplot() +
    theme_bw() +
    geom_line(
      aes(y=Q95_exp, x=year, colour = model), data=xx
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
      aes(y=Q95_obs, x=year), #, alpha="Observed"),
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
  p
  ## ggsave(file.path(plot_outputdir, paste0("ts_plot_", stn, ".png")), width = 5, height = 5, units = "in")
}

## Top skill scores
skill_scores %>% arrange(desc(msss))
p1 <- plotfun2(output, 18003)
p2 <- plotfun2(output, 46005)
p3 <- plotfun2(output, 79005)
p4 <- plotfun2(output, 73005)
p5 <- plotfun2(output, 76007)
p6 <- plotfun2(output, 80001)

plotfun3 <- function(x, ...) {
  ## Spatial plot (ACC)
  rdbu_pal = RColorBrewer::brewer.pal(9, "RdBu")
  p = ggplot() +
    geom_sf(
      data = europe_boundary,
      color=NA,
      fill="lightgrey"
    ) +
    geom_sf(
      data = uk_boundary,
      lwd = 0.25
    ) +
    geom_sf(
      data = gauge_stns %>% left_join(x, by = "ID"),
      aes(fill = skill),
      shape = 21,
      size = 1.5,
      lwd = 0.1,
      alpha = 0.8
    ) +
    ## facet_wrap(. ~ period, ncol = 1) +
    coord_sf(
      xlim = c(-8, 2),
      ylim = c(50, 59),
      default_crs = st_crs(4326)
    ) +
    ## scale_shape_manual(values = c(21, 24, 22)) +
    scale_fill_stepsn(
      colours = rdbu_pal,
      ## breaks = seq(-0.8, 0.8, 0.2),
      ## values = scales::rescale(c(-0.8, 0, 0.8)),
      ## limits = c(-0.3, 0.9)
      breaks = seq(-0.2, 0.8, 0.2),
      values = scales::rescale(c(-0.2, 0, 0.8)),
      limits = c(-0.1, 0.9)
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      ## legend.position = "bottom",
      ## legend.box = "vertical",
      ## legend.justification = "left",
      ## legend.box.just = "left",
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_label_size),
      strip.text = element_blank(),
      panel.grid.major = element_line(size = 0.25),
      axis.text = element_text(size = axis_label_size_small)
    ) #+
    ## guides(
    ##   shape = guide_legend(
    ##     title = "Model",
    ##     title.position = "top",
    ##     order = 1
    ##   ),
    ##   fill = guide_colorbar(
    ##     title="MSSS",
    ##     title.position="top",
    ##     frame.colour = "black",
    ##     ticks.colour = "black",
    ##     frame.linewidth = 0.25,
    ##     ticks.linewidth = 0.25,
    ##     barwidth = 0.75,
    ##     barheight = 10,
    ##     order = 2
    ##   )
    ## )
  p
}

skill_scores_i = skill_scores %>% filter(model %in% "TabNet") %>% mutate(skill = acc)
p1 <- plotfun3(skill_scores_i)
skill_scores_i = skill_scores %>% filter(model %in% "XGBoost") %>% mutate(skill = acc)
p2 <- plotfun3(skill_scores_i)

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

## x <- rbind(xgboost_output, tabnet_output)
## x <- x %>% mutate(date = as.Date(date))

## y <- open_dataset(gamlss_inputdir) %>% collect()
## station_ids <- y$ID %>% unique()
## n_stations <- length(station_ids)

## ## msss_list <- list()
## ## for (i in 1:n_stations) {
## ##   stn <- station_ids[i]
## ##   catchment_area <- metadata %>% filter(id %in% stn) %>% `$`(catchment_area)

## ##   ## Read GAMLSS model prediction
## ##   yy <- y %>% filter(ID %in% stn) %>% rename(gamlss_exp = exp) %>% dplyr::select(year, obs, gamlss_exp)
## ##   yy$obs <- yy$obs * 86400 / catchment_area * 1000 / 1000 / 1000
## ##   yy$gamlss_exp <- yy$gamlss_exp * 86400 / catchment_area * 1000 / 1000 / 1000

## ##   ## Read LSTM prediction
## ##   xx <- read_csv(file.path(nh_inputdir, paste0(stn, ".csv")), show_col_types = FALSE)
## ##   xx <- xx %>%
## ##     mutate(year = lubridate::year(date)) %>%
## ##     rename(obs = Q95_obs, nh_exp = Q95_sim) %>%
## ##     dplyr::select(year, obs, nh_exp, -date, -time_step) %>%
## ##     mutate(year = year - 1) # FIXME
## ##     ## gather(-season_year, key = "key", value = "value")

## ##   ## Compute MSSS
## ##   idx <- is.na(yy$gamlss_exp) | is.na(xx$nh_exp)
## ##   gamlss_msss <- mean_square_error_skill_score(yy$obs[!idx], yy$gamlss_exp[!idx])$msss
## ##   lstm_msss <- mean_square_error_skill_score(xx$obs[!idx], xx$nh_exp[!idx])$msss
## ##   msss_list[[i]] <- tibble(ID = stn, GAMLSS=gamlss_msss, LSTM=lstm_msss)

## ##   ## Make plot
## ##   xx <- xx %>%
## ##     left_join(yy %>% dplyr::select(-obs), by = c("year")) %>%
## ##     gather(-year, key = "key", value = "value") %>%
## ##     filter(!key %in% "obs")

## ##   xx <- xx %>% mutate(key = factor(key, levels = c("nh_exp", "gamlss_exp"), labels = c("LSTM", "GAMLSS")))
## ##   p <- ggplot() +
## ##     theme_bw() +
## ##     geom_line(
## ##       aes(y=value, x=year, colour = key), data=xx
## ##     ) +
## ##     ## scale_fill_manual(values = cbbPalette) +
## ##     ## scale_color_manual(values = cbbPalette) +
## ##     ## scale_fill_discrete(values = cbbPalette) +
## ##     ## scale_colour_discrete(values = cbbPalette) +
## ##     ## facet_wrap(. ~ ID, ncol = 1) + #, labeller = label_parsed) +
## ##     ylab(expression(Streamflow~(mm~d^{-1}))) +
## ##     xlab("") +
## ##     ## scale_x_continuous(breaks = pretty_breaks()) +
## ##     ## ## scale_y_continuous(expression(Streamflow~(m^{3}~s^{-1})), breaks = pretty_breaks()) +
## ##     ## ## N.B. use alpha to create another legend
## ##     ## ## https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
## ##     ## ## geom_point(
## ##     geom_line(
## ##       aes(y=obs, x=year), #, alpha="Observed"),
## ##       color = "black",
## ##       data=yy,
## ##       size = 1 #0.2
## ##     ) +
## ##     ## scale_alpha_manual(name=NULL, values=1, breaks="Observed") +
## ##     ## ggtitle(sprintf("ID = %d", stn_id)) +
## ##     theme(legend.position = "bottom",
## ##           legend.direction = "vertical",
## ##           legend.title = element_blank(),
## ##           strip.background = element_blank(),
## ##           panel.grid = element_blank(),
## ##           strip.text = element_text(size = strip_label_size),
## ##           ## legend.title = element_text(size = legend_title_size),
## ##           legend.text = element_text(size = legend_label_size),
## ##           axis.title.y = element_blank(),
## ##           axis.text.y = element_text(size = axis_label_size_small),
## ##           axis.text.x = element_text(size = axis_label_size_small))
## ##   ggsave(file.path(plot_outputdir, paste0("ts_plot_", stn, ".png")), width = 5, height = 5, units = "in")
## ## }

## ## msss <- do.call("rbind", msss_list)

## ## length(which(msss$LSTM > msss$GAMLSS))
## ## length(which(msss$LSTM < 0))
## ## msss %>% arrange(desc(GAMLSS))
## ## msss %>% arrange(desc(LSTM))
## ## msss %>% arrange(LSTM)

## ## msss <- msss %>% filter(LSTM > 0)
## ## plot(msss$GAMLSS, msss$LSTM, xlim = c(0, 1), ylim = c(0, 1))
## ## boxplot(msss$LSTM)
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

