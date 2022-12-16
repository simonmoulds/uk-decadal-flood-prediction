#!/usr/bin/env Rscript

## Author : Simon Moulds
## Date   : Jan 2022

library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(scales)
library(arrow)
library(sf)
library(smoothr)
## library(gamlss)
## library(rnrfa)
library(yaml)

options(dplyr.summarise.inform = FALSE)
options(bitmapType = 'cairo') # For server

## ## FOR TESTING:
## config = read_yaml('config/config.yml')
## aggregation_period = "yr2to5_lag"
## outputroot = 'results/exp1'
## cwd = 'workflow/scripts'

config = read_yaml("config/config_2.yml")
aggregation_period = "yr2"
outputroot = "results/exp2"
cwd = "workflow/scripts"

## ## Extract configuration info
## if (sys.nframe() == 0L) {
##   args = commandArgs(trailingOnly=TRUE)
##   config = read_yaml(args[1])
##   aggregation_period = args[2]
##   outputroot = args[3]
##   args = commandArgs()
##   m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
##   cwd <- dirname(regmatches(args, m))
## }
source(file.path(cwd, "utils.R"))
source(file.path(cwd, "plotting.R"))
config[["modelling"]] <- parse_config_modelling(config)

skill_measure = "crpss"
all_skill_measures <- c("crps_fcst", "crps_ens_fcst", "crps_climat", "crpss", "aic")

output_dir = file.path(outputroot, "fig", aggregation_period)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

aggregation_period_label_list = list(
  "yr2" = "Year 2",
  "yr2to5_lag" = "Year 2-5",
  "yr6to9_lag" = "Year 6-9",
  "yr2to9_lag" = "Year 2-9"
)

obs_aggregation_period_list = list(
  "yr2" = "yr2",
  "yr2to5_lag" = "yr2to5",
  "yr6to9_lag" = "yr6to9",
  "yr2to9_lag" = "yr2to9"
)

aggregation_period_label = aggregation_period_label_list[[aggregation_period]]
obs_aggregation_period = obs_aggregation_period_list[[aggregation_period]]
obs_aggregation_period_label = aggregation_period_label

## ####################################################### ##
## ####################################################### ##
##
## Preamble
##
## ####################################################### ##
## ####################################################### ##

## In this analysis we have many different models
## so we need to load data slighly differently
models <- sapply(config$modelling, FUN=function(x) x$name)
predictand <- sapply(config$modelling, FUN=function(x) x$predictand)
aggregation_period <- sapply(config$modelling, FUN=function(x) x$aggregation_period)

all_skill_scores <- list()
all_predictions <- list()
for (i in 1:length(models)) {
  skill_scores_i <- load_skill_scores(config, models[i], aggregation_period[i])
  skill_scores_i <- skill_scores_i %>% mutate(experiment = models[i])
  predictions_i <- load_model_predictions(config, models[i], aggregation_period[i])
  cor <- predictions_i %>%
    group_by(ID, subset, model) %>%
    summarize(R=cor.test(observations, Q50)$estimate)
  skill_scores_i <- skill_scores_i %>% left_join(cor, by = c("ID", "subset", "model"))
  predictions_i <- predictions_i %>% mutate(experiment = models[i])
  all_skill_scores[[i]] <- skill_scores_i
  all_predictions[[i]] <- predictions_i
}

names(all_skill_scores) <- models
names(all_predictions) <- models

## skill_scores <- do.call("rbind", all_skill_scores)
## predictions <- do.call("rbind", all_predictions)
station_ids <- skill_scores$ID %>% unique()

## Compute R value

## obs_depvar = "Q_95_centred_obs"
## exp_depvar = "Q_95_centred_exp"
## obs_expm = "observed_Q95"
## exp_expm = "hindcast_Q95"
## obs_depvar = "Q_max_centred_obs"
## exp_depvar = "Q_max_centred_exp"
## obs_expm = "observed_QMAX"
## exp_expm = "hindcast_QMAX"

## ## Load model skill scores for observed and hindcast experiments
## obs_skill_scores <- load_skill_scores(config, obs_expm, obs_aggregation_period)
## skill_scores <- load_skill_scores(config, exp_expm, aggregation_period)

## ## Load model fit metrics for hindcast experiment
## fit <- load_model_fit(config, "hindcast", aggregation_period) %>% mutate(kurtosis = kurtosis - 3)

## For spatial plots:
uk_boundary =
  st_read("../data-raw/CNTR_RG_01M_2020_4326.shp") %>%
  filter(CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

europe_boundary =
  st_read("../data-raw/CNTR_RG_01M_2020_4326.shp") %>%
  filter(!CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

gauge_stns =
  catalogue() %>%
  rename(ID = id, area = "catchment-area") %>%
  filter(ID %in% station_ids) %>%
  dplyr::select(ID, name, area, latitude, longitude) %>%
  st_as_sf(coords=c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_transform(27700)

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

## ####################################################### ##
## ####################################################### ##
##
## Analysis
##
## ####################################################### ##
## ####################################################### ##

## ## Find out which observed stations are positive
## obs_skill_subset <-
##   skill_scores %>%
##   filter(subset %in% "observed") %>%
##   mutate(model = droplevels(model)) %>%
##   group_by(ID, experiment) %>%
##   filter(crps_ens_fcst==min(crps_ens_fcst))

## ## obs_ids <- obs_skill_subset %>% filter(crpss>0) %>% `$`(ID)

## ## ## Median skill using observed predictors
## ## obs_skill_subset %>% `$`(crpss) %>% median()
## ## obs_skill_subset %>% filter(crpss > 0) %>% `$`(crpss) %>% median()

## full_ens_ids <-
##   skill_scores %>%
##   ## filter(model %in% c("P", "PT") & subset %in% "full") %>%
##   filter(subset %in% "full") %>%
##   mutate(model = droplevels(model)) %>%
##   group_by(ID) %>%
##   filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
##   filter(crpss>0) %>% `$`(ID)

## nao_matched_ids <-
##   skill_scores %>%
##   ## filter(model %in% c("P", "PT") & subset %in% "best_n") %>%
##   filter(subset %in% "best_n") %>%
##   mutate(model = droplevels(model)) %>%
##   group_by(ID) %>%
##   filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
##   filter(crpss>0) %>% `$`(ID)

## ## Percent full ensemble IDs with +ve skill in observed IDs with +ve skill
## sum(obs_ids %in% full_ens_ids) / length(obs_ids)
## sum(obs_ids %in% nao_matched_ids) / length(obs_ids)

## ## Overall statistic
## ## % stations with +ve MSSS
## stat <-
##   skill_scores %>%
##   filter(period %in% aggregation_period & model %in% c("P", "PT")) %>%
##   group_by(ID, subset) %>%
##   filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
##   ungroup() %>%
##   group_by(subset) %>%
##   summarize(pct_positive = sum(crpss > 0) / n() * 100)

## ## Time series plot
## skill_scores_subset =
##   skill_scores %>%
##   filter(subset %in% c("best_n", "full")) %>%
##   filter(period %in% aggregation_period) %>%
##   mutate(skill = !!sym(skill_measure)) %>%
##   dplyr::select(-all_of(all_skill_measures)) %>%
##   pivot_wider(names_from = subset, values_from = skill) %>%
##   mutate(skill_diff = best_n - full)

## skill_scores_subset <-
##   skill_scores_subset %>% ungroup() %>%
##   left_join(
##     (skill_scores %>%
##      filter(subset %in% "best_n") %>%
##      filter(period %in% aggregation_period) %>%
##      mutate(skill = !!sym(skill_measure)) %>%
##      dplyr::select(model, ID, period, crps_ens_fcst, skill)),
##     by = c("model", "ID", "period")
##   )

## skill =
##   skill_scores_subset %>%
##   ## filter(model %in% c("P", "PT")) %>%
##   mutate(model = droplevels(model)) %>%
##   group_by(ID) %>%
##   filter(crps_ens_fcst == min(crps_ens_fcst)) %>%
##   arrange(desc(best_n))

## ## skill =
## ##   skill_scores_subset %>%
## ##   filter(model %in% c("P", "PT")) %>%
## ##   group_by(ID) %>%
## ##   filter(skill_diff == max(skill_diff)) %>%
## ##   arrange(desc(best_n))

## ids_best = head(skill, n=20)$ID
## models_best = as.character(head(skill, n=20)$model)

## predictions <- load_model_predictions(config, exp_expm, aggregation_period)
## simulations <- load_model_simulations(config, exp_expm, aggregation_period)

## ## predictions = open_dataset(
## ##   file.path(outputroot, "analysis", "hindcast", "gamlss", aggregation_period, "prediction")
## ## ) %>%
## ##   collect() %>%
## ##   filter(subset %in% c("full", "best_n")) %>%
## ##   mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))

## ## model_levels <- unique(predictions$model)
## ## model_labels <- model_levels %>% gsub("_", "", .)
## ## predictions <- predictions %>% mutate(model = factor(model, levels = model_levels, labels = model_labels))
## predictions <- predictions %>% mutate(obs = !!sym(obs_depvar), exp = !!sym(exp_depvar))
## simulations <- simulations %>% mutate(obs = !!sym(obs_depvar), exp = !!sym(exp_depvar))

## ## library(viridis)
## ## p1 <- predictions %>%
## ## p1 <- simulations %>%
## ##   filter(test_year %in% 1985) %>%
## ##   filter(ID %in% ids_best[1] & model %in% models_best[1]) %>%
## ##   myplotfun66()

## myplotfun_ts <- function(simulation, prediction, cutoff) {
##   p <- ggplot(data=simulation, aes(x=year))+
##     theme_bw()+ xlab("")+
##     ylab(bquote('DJFM MAX ('*m^3*'/'*s*')'))+ #'/'*km^2*')'))+ #instantaneous peak flow
##     geom_ribbon(aes(ymin=Q01, ymax=Q99),  fill="ivory2")+
##     geom_line(aes(y=Q01),col="grey50", size=0.5)+#lower is darker
##     geom_line(aes(y=Q95),col="grey50", size=0.5)+
##     geom_line(aes(y=Q98),col="red", size=0.5)+
##     geom_line(aes(y=Q99),col="grey50", size=0.5)+
##     geom_line(aes(y=Q50),col="blue",size=1)+
##     geom_point(aes(y=obs),pch=21, color="black", fill="grey50", size=2)+

##     # Add the predictions
##     geom_ribbon(data=prediction, aes(x=year, ymin=Q01, ymax=Q99),  fill="cornsilk")+
##     geom_line(data=prediction, aes(y=Q01),col="grey50", size=0.5, linetype="dashed")+#lower is darker
##     geom_line(data=prediction, aes(y=Q95),col="grey50", size=0.5, linetype="dashed")+
##     geom_line(data=prediction, aes(y=Q98),col="red", size=0.5, linetype="dashed")+
##     geom_line(data=prediction, aes(y=Q99),col="grey50", size=0.5, linetype="dashed")+
##     geom_line(data=prediction, aes(y=Q50),col="blue",size=1)+
##     geom_point(data=prediction, aes(y=obs),pch=21, color="black", fill="grey50", size=2)+
##     scale_x_continuous(lim=c(1886,2024),
##                        breaks=seq(1850,2020,10),
##                        expand=c(0,0)) +
##     coord_cartesian(xlim=c(1967,2015))
##     ## coord_cartesian(ylim=c(0,100), xlim=c(1967,2005))+

##     ## ## Add some annotations
##     ## annotate(geom ='text', label = "Red line indicates nonstationary 50-y flood",
##     ##          x = 1970, y = 80,  size=4.5, col="black", hjust = 0)+#hjust 0 is left align
##     ## annotate(geom ='text', label = "in every year",
##     ##          x = 1970, y = 75,  size=4.5, col="black", hjust = 0)+
##     ## theme(legend.position = c(0.2, 0.8))
##   p
## }

## ## acc = predictions %>% group_by(ID, subset, model) %>% summarize(acc=cor.test(observations, Q50)$estimate) %>% arrange(desc(acc))
## ## ids_best = acc$ID[1:5]
## ## models_best = acc$model[1:5]

## cutoff <- 1991
## sim <- simulations %>%
##   filter(test_year %in% cutoff) %>%
##   filter(ID %in% ids_best[1] & model %in% models_best[1] & subset %in% "best_n")
## pred <- predictions %>%
##   filter(ID %in% ids_best[1] & model %in% models_best[1] & subset %in% "best_n") %>%
##   filter(year >= cutoff)

## p <- myplotfun_ts(sim, pred, cutoff)

## skill_scores_subset =
##   skill_scores %>%
##   ## filter(model %in% c("P", "PT")) %>% #, "NAOPT")) %>%
##   mutate(model = droplevels(model)) %>%
##   filter(subset %in% "best_n", period %in% aggregation_period) %>%
##   mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
##   mutate(skill = !!sym(skill_measure))

## skill = myfun(skill_scores_subset)

## ## Let's look at ACC
## acc <- predictions %>% group_by(ID, subset, model) %>% summarize(acc=cor.test(observations, Q50)$estimate) %>% arrange(desc(acc))
## skill <- skill %>% left_join(acc, by = c("ID", "subset", "model")) %>% mutate(skill=acc)

myplotfun1 <- function(x, legend_title = "ACC") {
  rdbu_pal = RColorBrewer::brewer.pal(12, "RdBu")
  p =
    ggplot() +
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
      data = x,
      aes(fill = skill),
      shape=21,
      size = 2,
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
    ## scale_fill_stepsn(
    ##   colours = rev(rdbu_pal),
    ##   breaks = seq(-0.2, 0.8, 0.2),
    ##   values = scales::rescale(c(-0.2, 0, 0.8)),
    ##   limits = c(-0.1, 0.9)
    ## ) +
    scale_fill_stepsn(
      ## colours = rev(rdbu_pal)[c(4:5, 7:11)], #[3:9],
      ## breaks = seq(-0.4, 0.8, 0.2),
      ## limits = c(-0.3, 0.9)
      colours = rev(rdbu_pal)[c(4:6, 7:11)], #[3:9],
      breaks = seq(-0.2, 0.6, 0.1),
      limits = c(-0.3, 0.7),
      labels=mylabelfun
    ) +
    theme_bw() +
    ## theme(
    ##   strip.background = element_blank(),
    ##   legend.position = "bottom",
    ##   legend.box = "vertical",
    ##   legend.justification = "left",
    ##   legend.box.just = "left",
    ##   legend.title = element_text(size = legend_title_size),
    ##   legend.text = element_text(size = legend_label_size),
    ##   strip.text = element_blank(),
    ##   panel.grid.major = element_line(size = 0.25),
    ##   axis.text = element_text(size = axis_label_size_small),
    ## ) +
    guides(
      ## shape = guide_legend(
      ##   title = "Model",
      ##   title.position = "top",
      ##   order = 1
      ## ),
      ## fill = guide_colorbar(
      ##   title=legend_title,
      ##   title.position="top",
      ##   frame.colour = "black",
      ##   ticks.colour = "black",
      ##   frame.linewidth = 0.25,
      ##   ticks.linewidth = 0.25,
      ##   barwidth = 12,
      ##   barheight = 0.75,
      ##   order = 2
      ## )
      fill = guide_colorbar(
        title=legend_title,
        title.position="top",
        frame.colour = "black",
        ticks.colour = "black",
        frame.linewidth = 0.25,
        ticks.linewidth = 0.25,
        ## barwidth = 12,
        ## barheight = 0.75,
        barwidth = 6,
        barheight = 0.4,
        order = 2
      ) ##+
  ## theme(legend.spacing.y = unit(0.2, "cm"),
  ##       legend.margin = margin(0, 0, -0.5, 0, unit="cm"))
    )

  ## ## Add model counts to top right corner of plot
  ## n_pt <- table(x$model)[["PT"]]
  ## n_p <- table(x$model)[["P"]]
  ## labs <- c(paste0("italic(n)==", n_p), paste0("italic(n)==", n_pt))
  ## d <- data.frame(x = c(0, 0), y = c(59, 58.5), lab = labs, model = c("P", "PT"))
  p <- p +
    ## geom_point(
    ##   data = d,
    ##   aes(x, y, shape = model),
    ##   size = 1.5,
    ##   lwd = 0.1,
    ##   show.legend = FALSE
    ## ) +
    ## geom_text(
    ##   data = d,
    ##   aes(x, y, label = lab),
    ##   parse = TRUE,
    ##   hjust = 0,
    ##   nudge_x = 0.3,
    ##   size = 2
    ## ) +
    ## scale_shape_manual(values = c(21, 24)) +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(size = axis_label_size),
          axis.text.y = element_text(size = axis_label_size),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.justification = "left",
          legend.box.just = "left",
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_label_size),
          strip.text = element_blank(),
          panel.grid.major = element_line(size = 0.25))
  p
}

skill1 <-
  ## all_skill_scores$hindcast_POT %>%
  all_skill_scores$hindcast_POT %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst))
skill1 <- gauge_stns %>% left_join(skill1)

skill2 <-
  all_skill_scores$hindcast_Q95 %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst))
skill2 <- gauge_stns %>% left_join(skill2)

p1 <- myplotfun1(skill1 %>% mutate(skill = R), "R")
p2 <- myplotfun1(skill2 %>% mutate(skill = R), "R")

p1 <- p1 + labs(title="DJFM POT [Q99]") + theme(plot.title = element_text(size = 8, margin = margin(0,0,1,0, unit="pt")))
p2 <- p2 + labs(title="DJFM Q95") + theme(plot.title = element_text(size = 8, margin = margin(0,0,1,0, unit="pt")))
p <- p1 + p2 + plot_layout(ncol=2) + plot_layout(guides="collect") & theme(legend.position = "bottom")
ggsave("results/exp2/fig/yr2/fig_map.png", width=5, height=5, units="in")

stop()

## p1 <- p1 +
##   labs(title = paste0("ID=", ids_best[1])) + #, " (", models_best[1], ")")) +
##   theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

## p2 <- predictions %>%
##   filter(ID %in% ids_best[2] & model %in% models_best[2]) %>%
##   myplotfun66()
## p2 <- p2 +
##   labs(title = paste0("ID=", ids_best[2])) + #, " (", models_best[2], ")")) +
##   theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

## p3 <- predictions %>%
##   filter(ID %in% ids_best[3] & model %in% models_best[3]) %>%
##   myplotfun66()
## p3 <- p3 +
##   labs(title = paste0("ID=", ids_best[3])) + #, " (", models_best[3], ")")) +
##   theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

## p4 <- predictions %>%
##   filter(ID %in% ids_best[4] & model %in% models_best[4]) %>%
##   myplotfun66()
## p4 <- p4 +
##   labs(title = paste0("ID=", ids_best[4])) + #, " (", models_best[4], ")")) +
##   theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

## p5 <- predictions %>%
##   filter(ID %in% ids_best[5] & model %in% models_best[5]) %>%
##   myplotfun66()
## p5 <- p5 +
##   labs(title = paste0("ID=", ids_best[5])) + #, " (", models_best[5], ")")) +
##   theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

## format_p_value <- function(p_value) {
##   if (p_value < 0.01) {
##     return("(P < 0.01)")
##   } else {
##     return(paste0("(P = ", sprintf(p_value, fmt = "%#.2f"), ")"))
##   }
## }

## make_annotation <- function(crpss, id, mod) {
##   crpss_full <- crpss %>% filter(ID %in% id & subset %in% "full" & model %in% mod) %>% `$`(crpss)
##   crpss_matched <- crpss %>% filter(ID %in% id & subset %in% "best_n" & model %in% mod) %>% `$`(crpss)
##   annotation = paste0(
##     paste0(
##       "CRPSS (Full) = ", sprintf(crpss_full, fmt = "%#.2f"), "\n"#, " ", format_p_value(acc_full$acc_p), "\n"
##     ),
##     paste0(
##       "CRPSS (NAO-matched) = ", sprintf(crpss_matched, fmt = "%#.2f")#, "\n"#, " ", format_p_value(acc_matched$acc_p), "\n"
##     )
##   )
##   annotation
## }

## get_y_range <- function(p) {
##   yrange <- layer_scales(p)$y$range$range
##   return(yrange)
## }
## get_y_position <- function(p, rel_pos) {
##   yrange <- get_y_range(p)
##   return(yrange[1] + diff(yrange) * rel_pos)
## }

## annotate_plot <- function(p, id, mod) {
##   annotation_size = 4
##   annotation_rel_pos = 1
##   yrange <- get_y_range(p)
##   yrange[2] <- yrange[2] * 1.025
##   p <- p +
##     ylim(yrange) +
##     annotate(
##       geom = "text",
##       x = 1960,
##       y = yrange[1] + diff(yrange) * annotation_rel_pos,
##       ## label = make_annotation(acc, ids_best[1], models_best[1]),
##       label = make_annotation(skill_scores, id, mod), #ids_best[1], models_best[1]),
##       hjust=0,
##       vjust=1,
##       size = annotation_size / ggplot2::.pt
##     ) +
##   ## p <- p +
##   ##   ylim(range) +
##     annotate(
##       geom = "text",
##       x = 2005,
##       y = yrange[1] + diff(yrange) * annotation_rel_pos,
##       label = mod,
##       hjust = 1,
##       vjust = 1,
##       size = annotation_size / ggplot2::.pt,
##       fontface = "bold"
##     )
##   p
## }
## p1 <- p1 %>% annotate_plot(ids_best[1], models_best[1])
## p2 <- p2 %>% annotate_plot(ids_best[2], models_best[2])
## p3 <- p3 %>% annotate_plot(ids_best[3], models_best[3])
## p4 <- p4 %>% annotate_plot(ids_best[4], models_best[4])
## p5 <- p5 %>% annotate_plot(ids_best[5], models_best[5])

## ## Plot location of the selected gauge stations
## gauge_stns_subset = gauge_stns %>% filter(ID %in% ids_best)
## p6 = myplotfun777(gauge_stns_subset)

## p1 = p1 + theme(plot.margin = margin(0, 0, 0, 0))
## p2 = p2 + theme(plot.margin = margin(0, 0, 0, 0.25, unit="cm"),
##                 axis.title.y = element_blank())
## p3 = p3 + theme(plot.margin = margin(0, 0, 0, 0),
##                 axis.title.y = element_blank())
## p4 = p4 + theme(plot.margin = margin(0, 0, 0, 0),
##                 axis.title.y = element_text(size = axis_title_size))
## p5 = p5 + theme(plot.margin = margin(0, 0, 0, 0.25, unit="cm"),
##                 axis.title.y = element_blank())
## p6 = p6 + theme(plot.margin = margin(0, 0, 0, 0),
##                 axis.title = element_blank(),
##                 axis.text = element_text(size = axis_label_size_small))

## p = p1 + p2 + p3 + p4 + p5 + p6 +
##   plot_layout(ncol = 3, nrow = 2, widths = c(2, 2, 2)) & theme(legend.position = "bottom")
## p = p + plot_layout(guides = "collect") & theme(legend.margin = margin(0, 0, 0, 0), legend.key.size = unit(1, 'lines'), legend.box.spacing = unit(0, 'cm'))
## ggsave(file.path(output_dir, "fig_ts.png"), plot = p, width = 6, height = 4.25, units = "in")

## ## #################################### ##
## ## Load data
## ## #################################### ##

## outputroot <- file.path(outputroot, "analysis", "hindcast")
## gamlss_datadir <- file.path(outputroot, "gamlss", aggregation_period, "prediction")
## ## lstm_datadir <- file.path(outputroot, "lstm", aggregation_period, "prediction")
## xgboost_datadir <- file.path(outputroot, "xgboost", aggregation_period, "prediction")
## tabnet_datadir <- file.path(outputroot, "tabnet", aggregation_period, "prediction")
## ## metadata <- read_parquet(file.path(outputroot, "nrfa-metadata.parquet"))

## ## Join data frames and reshape
## gamlss_output <- open_dataset(gamlss_datadir) %>%
##   collect() %>%
##   rename(Q95_exp = Q_95_exp, Q95_obs = Q_95_obs) %>%
##   dplyr::select(Q95_obs, Q95_exp, year, ID, date, model)

## ## lstm_output <- open_dataset(lstm_datadir) %>%
## ##   collect() %>%
## ##   rename(Q95_exp = Q95_sim) %>%
## ##   mutate(model = "LSTM")

## xgboost_output <- open_dataset(xgboost_datadir) %>%
##   collect() %>%
##   mutate(model = "XGBoost")

## tabnet_output <- open_dataset(tabnet_datadir) %>%
##   collect() %>%
##   mutate(model = "TabNet")

## #################################### ##
## Get skill scores for each model
## #################################### ##

## compute_skill_scores <- function(x, ...) {
##   models <- unique(x$model) %>% sort()
##   station_ids <- unique(x$ID) %>% sort()
##   skill_scores <- list()
##   for (i in 1:length(models)) {
##     xx <- x %>% filter(model %in% models[i])
##     model_skill_scores <- list()
##     print(paste0("Computing skill scores for model ", models[i]))
##     pb <- txtProgressBar(min = 0, max = length(station_ids), initial = 0)
##     for (j in 1:length(station_ids)) {
##       id <- station_ids[j]
##       yy <- xx %>% filter(ID %in% id) %>% arrange(year)
##       skill <- mean_square_error_skill_score(yy$Q95_obs, yy$Q95_exp) %>% as_tibble()
##       skill <- skill %>% mutate(ID = id, .before = names(.)[1])
##       model_skill_scores[[j]] <- skill
##       setTxtProgressBar(pb, j)
##     }
##     close(pb)
##     model_skill_scores <- do.call("rbind", model_skill_scores) %>% mutate(model = models[i], .before = ID)
##     skill_scores[[i]] <- model_skill_scores
##   }
##   skill_scores <- do.call("rbind", skill_scores)
##   skill_scores
## }

## gamlss_skill_scores <-
##   compute_skill_scores(gamlss_output) %>%
##   arrange(ID)

## gamlss_skill_scores_best <-
##   gamlss_skill_scores %>%
##   group_by(ID) %>%
##   filter(acc == max(acc))

## gamlss_output <-
##   gamlss_skill_scores_best %>%
##   dplyr::select(model, ID) %>%
##   left_join(gamlss_output)

## output <- rbind(
##   xgboost_output,
##   tabnet_output,
##   ## lstm_output,
##   gamlss_output
## ) %>%
##   mutate(date = as.Date(date)) %>%
##   arrange(ID, date, model) %>%
##   mutate(submodel = model) %>% # So we don't lose the information completely
##   mutate(model = ifelse(str_detect(model, "GAMLSS_(.*)_(full|best_n)"), "GAMLSS", model))

## ## lstm_skill_scores <- compute_skill_scores(lstm_output)
## xgboost_skill_scores <- compute_skill_scores(xgboost_output)
## tabnet_skill_scores <- compute_skill_scores(tabnet_output)

## ## Combine skill scores in a single data frame
## skill_scores <- rbind(
##   xgboost_skill_scores,
##   tabnet_skill_scores,
##   ## lstm_skill_scores,
##   gamlss_skill_scores_best
## )

## ## Provide a catch-all model name for GAMLSS
## skill_scores <-
##   skill_scores %>%
##   mutate(submodel = model) %>% # So we don't lose the information completely
##   mutate(model = ifelse(str_detect(model, "GAMLSS_(.*)_(full|best_n)"), "GAMLSS", model))

## ## Only include stations where all models are available
## models <- c("GAMLSS", "TabNet", "XGBoost")
## n_models <- length(models)
## common_ids <-
##   skill_scores %>%
##   group_by(ID) %>%
##   summarize(n = n()) %>%
##   filter(n == 3) %>% `$`(ID)

## skill_scores <- skill_scores %>% filter(ID %in% common_ids)

## ## Add the best model
## skill_scores_best <- skill_scores %>% group_by(ID) %>% filter(acc == max(acc)) %>% mutate(model = "Best")
## skill_scores <- rbind(skill_scores, skill_scores_best)

## station_ids <- unique(skill_scores$ID) %>% sort()
## skill_scores <- skill_scores %>% mutate(model = factor(model, levels = c("GAMLSS", "XGBoost", "TabNet", "Best")))

## ## ################################### ##
## ## Plot 1: Boxplot comparison of models
## ## ################################### ##

## ## Need gauge station locations
## gauge_stns =
##   catalogue() %>%
##   rename(ID = id, area = "catchment-area") %>%
##   filter(ID %in% station_ids) %>%
##   dplyr::select(ID, name, area, latitude, longitude) %>%
##   st_as_sf(coords=c("longitude", "latitude")) %>%
##   st_set_crs(4326) %>%
##   st_transform(27700)

## plotfun1 <- function(x, legend_title = "ACC") {
##   rdbu_pal = RColorBrewer::brewer.pal(9, "RdBu")
##   p =
##     ggplot() +
##     geom_sf(
##       data = europe_boundary,
##       color=NA,
##       fill="lightgrey"
##     ) +
##     geom_sf(
##       data = uk_boundary,
##       lwd = 0.25
##     ) +
##     geom_sf(
##       data = x,
##       aes(fill = skill, shape = model),
##       size = 1.5,
##       lwd = 0.1,
##       alpha = 0.8
##     ) +
##     ## facet_wrap(. ~ period, ncol = 1) +
##     coord_sf(
##       xlim = c(-8, 2),
##       ylim = c(50, 59),
##       default_crs = st_crs(4326)
##     ) +
##     scale_shape_manual(values = c(21, 24, 22)) +
##     scale_fill_stepsn(
##       colours = rev(rdbu_pal),
##       ## breaks = seq(-0.8, 0.8, 0.2),
##       ## values = scales::rescale(c(-0.8, 0, 0.8)),
##       ## limits = c(-0.3, 0.9)
##       breaks = seq(-0.2, 0.8, 0.2),
##       values = scales::rescale(c(-0.2, 0, 0.8)),
##       limits = c(-0.1, 0.9)
##     ) +
##     theme_bw() +
##     theme(
##       strip.background = element_blank(),
##       legend.position = "bottom",
##       legend.box = "vertical",
##       legend.justification = "left",
##       legend.box.just = "left",
##       legend.title = element_text(size = legend_title_size),
##       legend.text = element_text(size = legend_label_size),
##       strip.text = element_blank(),
##       panel.grid.major = element_line(size = 0.25),
##       axis.text = element_text(size = axis_label_size_small),
##     ) +
##     guides(
##       shape = guide_legend(
##         title = "Model",
##         title.position = "top",
##         order = 1
##       ),
##       fill = guide_colorbar(
##         title=legend_title,
##         ## title="MSSS",
##         title.position="top",
##         frame.colour = "black",
##         ticks.colour = "black",
##         frame.linewidth = 0.25,
##         ticks.linewidth = 0.25,
##         barwidth = 12,
##         barheight = 0.75,
##         order = 2
##       )
##     )
##   p
## }

## skill_scores <- skill_scores %>% mutate(skill = acc)
## skill_scores_subset <- skill_scores %>% filter(!model %in% "Best") %>% group_by(ID) %>% filter(skill == max(skill))
## skill_scores_subset <- gauge_stns %>% left_join(skill_scores_subset, by = "ID")
## p1 <- plotfun1(skill_scores_subset)

## ## Add number of stations which perform best by model
## n_best_model <- table(skill_scores_subset$model)
## n_gamlss <- n_best_model[["GAMLSS"]]
## n_xgboost <- n_best_model[["XGBoost"]]
## n_tabnet <- n_best_model[["TabNet"]]
## labs <- c(paste0("italic(n)==", n_gamlss), paste0("italic(n)==", n_xgboost), paste0("italic(n)==", n_tabnet))
## d <- data.frame(x = c(0, 0, 0), y = c(59, 58.5, 58), lab = labs, model = c("GAMLSS", "XGBoost", "TabNet"))
## p1 <- p1 +
##   geom_point(data = d, aes(x, y, shape = model), size = 1, lwd = 0.1, show.legend = FALSE) +
##   geom_text(data = d, aes(x, y, label = lab), parse = TRUE, hjust = 0, nudge_x = 0.3, size = 2) +
##   scale_shape_manual(values = c(21, 24, 22)) +
##   theme(axis.title = element_blank())

## rdbu_pal = brewer.pal(9, "RdBu")
## p1 <- p1 +
##   scale_fill_stepsn(
##     colours = rev(rdbu_pal)[3:9],
##     breaks = seq(-0.4, 0.8, 0.2),
##     limits = c(-0.3, 0.9)
##   )
## ## p1

## ## skill_scores_i = skill_scores %>% filter(model %in% "GAMLSS") %>% mutate(skill = acc)
## ## p1 <- plotfun1(skill_scores_i)
## ## skill_scores_i = skill_scores %>% filter(model %in% "TabNet") %>% mutate(skill = acc)
## ## p2 <- plotfun1(skill_scores_i)
## ## skill_scores_i = skill_scores %>% filter(model %in% "XGBoost") %>% mutate(skill = acc)
## ## p3 <- plotfun1(skill_scores_i)

## ## plotfun1 <- function(x, ...) {
## ##   ## Boxplot comparison
## ##   p = ggplot(x, aes(x = model, y=acc)) +
## ##     ## stat_boxplot(coef = NULL) +
## ##     geom_boxplot(
## ##       lwd = 0.25,
## ##       outlier.size = 0.25
## ##     ) +
## ##     scale_x_discrete(name = "") +
## ##     ## scale_y_continuous(
## ##     ##   name=legend_title,
## ##     ##   ## breaks=seq(-0.2, 1, by=0.2),
## ##     ##   ## limits=c(-0.25, 1.1)
## ##     ##   breaks=seq(-1, 1, by=0.2),
## ##     ##   limits=c(-1, 1)
## ##     ## ) +
## ##     theme_bw() +
## ##     theme(
## ##       panel.grid = element_blank(),
## ##       axis.title = element_text(size = axis_title_size_small),
## ##       axis.text = element_text(size = axis_label_size_small)
## ##     )
## ##   p
## ## }
## ##

## plotfun2 <- function(x, legend_title = "ACC") {
##   p <- ggplot(x, aes(x = model, y = skill)) + #interaction(subset, model), y = skill)) +
##     ## geom_boxplot(aes(colour = subset), lwd = 0.25, outlier.size = 0.25) +
##     geom_boxplot(lwd = 0.25, outlier.size = 0.25) +
##     ## geom_line(aes(group = ID), lwd = 0.15, alpha = 0.35, colour = "darkgrey") +
##     ## scale_colour_discrete(
##     ##   name = "",
##     ##   labels = c("Full ensemble", "NAO-matched ensemble"),
##     ##   type = c("#FC8D62", "#66C2A5")
##     ##   ## type = rev(cbbPalette)
##     ## ) +
##     ## geom_signif(
##     ##   ## comparisons = list(c("P", "PT"), c("PT", "NAOPT"), c("P", "NAOPT")),
##     ##   comparisons = list(c("best_n", "full")),
##     ##   map_signif_level = TRUE,
##     ##   step_increase = 0.06,
##     ##   tip_length = 0.01,
##     ##   size = 0.25,
##     ##   textsize = 1.5
##     ## ) +
##     ## scale_x_discrete(name = "") +
##     theme_bw() +
##     theme(
##       panel.grid = element_blank(),
##       axis.title = element_text(size = axis_title_size_small),
##       axis.text = element_text(size = axis_label_size_small)
##     )
##   p
## }

## ## Signif [one-sided Wilcoxon signed rank test]
## format_p <- function(x) {
##   if (x <= 0.001) return("***")
##   if (x <= 0.01) return("**")
##   if (x <= 0.05) return("*")
##   if (x > 0.05) return("NS.")
## }

## print("e")
## skill_measure = "acc"
## x <- skill_scores %>%
##   filter(model %in% "GAMLSS") %>%
##   `[[`(skill_measure)
## y <- skill_scores %>%
##   filter(model %in% "XGBoost") %>%
##   `[[`(skill_measure)
## pval1 = wilcox.test(x, y, paired = TRUE, alternative = "two.sided")$p.value %>% format_p

## x <- skill_scores %>%
##   filter(model %in% "GAMLSS") %>%
##   `[[`(skill_measure)
## y <- skill_scores %>%
##   filter(model %in% "TabNet") %>%
##   `[[`(skill_measure)
## pval2 = wilcox.test(y, x, paired = TRUE, alternative = "two.sided")$p.value %>% format_p

## x <- skill_scores %>%
##   filter(model %in% "XGBoost") %>%
##   `[[`(skill_measure)
## y <- skill_scores %>%
##   filter(model %in% "TabNet") %>%
##   `[[`(skill_measure)
## pval3 = wilcox.test(y, x, paired = TRUE, alternative = "two.sided")$p.value %>% format_p

## print("f")
## skill_scores_median <-
##   skill_scores %>%
##   group_by(model) %>%
##   summarize(md = median(skill), n = n()) %>%
##   mutate(md_label = sprintf(md, fmt = "%#.3f"))

## print("g")
## p2 <- plotfun2(skill_scores)
## print("h")
## p2 <- p2 +
##   geom_signif(
##     y_position = c(0.84, 0.92, 0.76), xmin = c(1, 1, 2), xmax = c(2, 3, 3),
##     annotation = c(pval1, pval2, pval3), #, format_p(pval3)),
##     step_increase = 0.06,
##     tip_length = 0.02,
##     size = 0.25,
##     textsize = 2) +
##   geom_text(data = skill_scores_median, aes(x = model, y = md, label = md_label), nudge_y = 0.05, size = 2) +
##   geom_text(data = skill_scores_median, aes(x = model, y = -0.6, label = n), nudge_y = -0.05, size = 2, fontface = "italic") +
##   scale_x_discrete(name = "", labels = levels(skill_scores$model)) +
##   scale_y_continuous(
##     name="ACC",
##     breaks=seq(-0.8, 1, by=0.2),
##     limits=c(-0.8, 1)
##   ) +
##   theme(legend.position = "bottom",
##         legend.title = element_blank(),
##         legend.text = element_text(size = legend_label_size),
##         legend.direction = "vertical",
##         legend.margin = margin(0, 0, 0, 0, unit = "cm"),
##         legend.box = "vertical",
##         legend.box.margin = margin(-0.5, 0, 0, 0, unit = "cm"),
##         aspect.ratio = 1
##         )

## print("i")
## p <- p1 + p2 + plot_layout(nrow = 1, ncol = 2)
## p <- p + plot_layout(widths = c(2, 2))
## p$patches$plots[[1]] =
##   p$patches$plots[[1]] +
##   labs(tag = "a") +
##   theme(plot.tag.position = c(0.1, 1.01),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p =
##   p +
##   labs(tag = "b") +
##   theme(plot.tag.position = c(0.125, 0.87),
##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))

## ggsave(file.path(plot_outputdir, "fig1.png"), plot = p, width = 6, height = 6, units = "in")

## plotfun2 <- function(x, id, ...) {
##   xx <- x %>% filter(ID %in% id)
##   models <- unique(xx$model)
##   yy <- xx %>% filter(model %in% models[1])
##   ## Time series plot
##   p <- ggplot() +
##     theme_bw() +
##     geom_line(
##       aes(y=Q95_exp, x=year, colour = model), data=xx
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
##       aes(y=Q95_obs, x=year), #, alpha="Observed"),
##       color = "black",
##       data=yy,
##       size = 0.6 #0.2
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
##   p
##   ## ggsave(file.path(plot_outputdir, paste0("ts_plot_", stn, ".png")), width = 5, height = 5, units = "in")
## }

## ## Top skill scores
## skill_scores_subset <- skill_scores %>% filter(!model %in% "Best") %>% arrange(desc(acc))
## ids_best = head(skill_scores_subset$ID, n=5)

## p1 <- plotfun2(output, ids_best[1])
## p2 <- plotfun2(output, ids_best[2])
## p3 <- plotfun2(output, ids_best[3])
## p4 <- plotfun2(output, ids_best[4])
## p5 <- plotfun2(output, ids_best[5])
## ## p6 <- plotfun2(output, 84019)

## format_p_value = function(p_value) {
##   if (p_value < 0.01) {
##     return("(P < 0.01)")
##   } else {
##     return(paste0("(P = ", sprintf(p_value, fmt = "%#.2f"), ")"))
##   }
## }

## ## Regionalization:
## ## ================
## ##
## ## See: https://stackoverflow.com/q/67579831 for extracting contours from raster
## ##      [raster::rasterToContour / grDevices::contourLines]
## ## See: https://swilke-geoscience.net/post/2020-09-10-kriging_with_r/kriging/ for kriging using gstat
## ##
## ## Basic idea:
## ## ===========
## ## * Use ordinary kriking to interpolate skill
## ## * Extract zero skill contour
## ## * Remove section of raster with zero skill
## ## * Interpolate over the resulting grid

## ## make_annotation = function(acc, id) {
## ##   acc_full <- acc %>% filter(ID %in% id & subset %in% "Full ensemble")
## ##   acc_matched <- acc %>% filter(ID %in% id & subset %in% "NAO-matched ensemble")
## ##   annotation = paste0(
## ##     paste0(
## ##       "ACC (Full) = ", sprintf(acc_full$acc, fmt = "%#.2f"), " ", format_p_value(acc_full$acc_p), "\n"
## ##     ),
## ##     paste0(
## ##       "ACC (NAO-matched) = ", sprintf(acc_matched$acc, fmt = "%#.2f"), " ", format_p_value(acc_matched$acc_p), "\n"
## ##     )
## ##   )
## ##   annotation
## ## }

## ## get_y_range <- function(p) {
## ##   yrange <- layer_scales(p)$y$range$range
## ##   return(yrange)
## ## }
## ## get_y_position <- function(p, rel_pos) {
## ##   yrange <- get_y_range(p)
## ##   return(yrange[1] + diff(yrange) * rel_pos)
## ## }

## ## annotation_size = 4
## ## annotation_rel_pos = 1
## ## yrange <- get_y_range(p1)
## ## yrange[2] <- yrange[2] * 1.025
## ## p1 <- p1 +
## ##   ylim(yrange) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960,
## ##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
## ##     label = make_annotation(acc, ids_best[1]),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = annotation_size / ggplot2::.pt
## ##   )

## ## yrange <- get_y_range(p2)
## ## yrange[2] <- yrange[2] * 1.025
## ## p2 <- p2 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960,
## ##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
## ##     label = make_annotation(acc, ids_best[2]),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = annotation_size / ggplot2::.pt
## ##   )

## ## yrange <- get_y_range(p3)
## ## yrange[2] <- yrange[2] * 1.025
## ## p3 <- p3 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960,
## ##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
## ##     label = make_annotation(acc, ids_best[3]),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = annotation_size / ggplot2::.pt
## ##   )

## ## yrange <- get_y_range(p4)
## ## yrange[2] <- yrange[2] * 1.025
## ## p4 <- p4 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960,
## ##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
## ##     label = make_annotation(acc, ids_best[4]),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = annotation_size / ggplot2::.pt
## ##   )

## ## yrange <- get_y_range(p5)
## ## yrange[2] <- yrange[2] * 1.025
## ## p5 <- p5 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960,
## ##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
## ##     label = make_annotation(acc, ids_best[5]),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = annotation_size / ggplot2::.pt
## ##   )
## myplotfun3 <- function(x) {
##   p =
##     ggplot() +
##     geom_sf(
##       data = europe_boundary,
##       color=NA,
##       fill="lightgrey"
##     ) +
##     geom_sf(
##       data = uk_boundary,
##       lwd = 0.25
##     ) +
##     geom_sf(
##       data = x,
##       ## shape = 21,
##       size = 0,
##       ## lwd = 0.1
##       ## alpha = 0.8
##     ) +
##     coord_sf(
##       xlim = c(-8, 2),
##       ylim = c(50, 60),
##       default_crs = st_crs(4326)
##     ) +
##     theme_bw() #+
##   p
## }

## gauge_stns_subset = gauge_stns %>% filter(ID %in% ids_best)
## p6 = myplotfun3(gauge_stns_subset)
## p6 =
##   p6 +
##   ggrepel::geom_label_repel(
##              data=gauge_stns_subset,
##              aes(label=ID, geometry=geometry),
##              stat="sf_coordinates",
##              min.segment.length = 0,
##              size = 1.75)

## p1 = p1 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = strip_label_size),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = legend_label_size),
##                 axis.title.y = element_text(size = axis_title_size, angle = 90),
##                 axis.text.y = element_text(size = axis_label_size_small),
##                 axis.text.x = element_text(size = axis_label_size_small))
## p2 = p2 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = strip_label_size),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = legend_label_size),
##                 ## axis.title.y = element_text(size = axis_title_size, angle = 90),
##                 axis.title.y = element_blank(),
##                 axis.text.y = element_text(size = axis_label_size_small),
##                 axis.text.x = element_text(size = axis_label_size_small))
## p3 = p3 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = strip_label_size),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = legend_label_size),
##                 ## axis.title.y = element_text(size = axis_title_size, angle = 90),
##                 axis.title.y = element_blank(),
##                 axis.text.y = element_text(size = axis_label_size_small),
##                 axis.text.x = element_text(size = axis_label_size_small))
## p4 = p4 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = strip_label_size),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = legend_label_size),
##                 axis.title.y = element_text(size = axis_title_size, angle = 90),
##                 ## axis.title.y = element_text(size = axis_title_size),
##                 axis.text.y = element_text(size = axis_label_size_small),
##                 axis.text.x = element_text(size = axis_label_size_small))
## p5 = p5 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = strip_label_size),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = legend_label_size),
##                 axis.title.y = element_blank(),
##                 axis.text.y = element_text(size = axis_label_size_small),
##                 axis.text.x = element_text(size = axis_label_size_small))
## p6 = p6 + theme(panel.grid.major = element_line(size = 0.25),
##                 axis.title = element_blank(),
##                 axis.text = element_text(size = axis_label_size_small))

## p = p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, nrow = 2, widths = c(2, 2, 2)) & theme(legend.position = "bottom")
## p = p + plot_layout(guides = "collect")

## ggsave(file.path(plot_outputdir, "fig2.png"), plot = p, width = 5, height = 5, units = "in")
