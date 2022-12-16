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
## config = read_yaml('config/config_1.yml')
## aggregation_period = "yr2to5_lag"
## outputroot = 'results/exp1'
## cwd = 'workflow/scripts'

## Extract configuration info
if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  config = read_yaml(args[1])
  aggregation_period = args[2]
  outputroot = args[3]
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
source(file.path(cwd, "utils.R"))
source(file.path(cwd, "plotting.R"))
config[["modelling"]] <- parse_config_modelling(config)

fig_dpi <- 600
skill_measure <- "crpss"
all_skill_measures <- c("crps_fcst", "crps_ens_fcst", "crps_climat", "crpss", "aic")

output_dir <- file.path(outputroot, "fig", aggregation_period)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

aggregation_period_label_list <- list(
  "yr2" = "Year 2",
  "yr2to5_lag" = "Year 2-5",
  "yr6to9_lag" = "Year 6-9",
  "yr2to9_lag" = "Year 2-9"
)

obs_aggregation_period_list <- list(
  "yr2" = "yr2",
  "yr2to5_lag" = "yr2to5",
  "yr6to9_lag" = "yr6to9",
  "yr2to9_lag" = "yr2to9"
)

aggregation_period_label <- aggregation_period_label_list[[aggregation_period]]
obs_aggregation_period <- obs_aggregation_period_list[[aggregation_period]]
obs_aggregation_period_label <- aggregation_period_label

## ####################################################### ##
## ####################################################### ##
##
## Preamble
##
## ####################################################### ##
## ####################################################### ##

## Load model skill scores for observed and hindcast experiments
obs_skill_scores <- load_skill_scores(config, "observed", obs_aggregation_period)
skill_scores <- load_skill_scores(config, "hindcast", aggregation_period)
station_ids <- skill_scores$ID %>% unique()

## Load model fit metrics for hindcast experiment
fit <- load_model_fit(config, "hindcast", aggregation_period) %>% mutate(kurtosis = kurtosis - 3)

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
## Load model skill and compute summary statistics
##
## ####################################################### ##
## ####################################################### ##

## Find out which observed stations are positive
obs_skill_subset <-
  obs_skill_scores %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID) %>%
  filter(crps_ens_fcst==min(crps_ens_fcst))

obs_ids <- obs_skill_subset %>% filter(crpss>0) %>% `$`(ID)

## Median skill using observed predictors
obs_skill_subset %>% `$`(crpss) %>% median()
obs_skill_subset %>% filter(crpss > 0) %>% `$`(crpss) %>% median()

full_ens_ids <-
  skill_scores %>%
  filter(model %in% c("P", "PT") & subset %in% "full") %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID) %>%
  filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
  filter(crpss>0) %>% `$`(ID)

nao_matched_ids <-
  skill_scores %>%
  filter(model %in% c("P", "PT") & subset %in% "best_n") %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID) %>%
  filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
  filter(crpss>0) %>% `$`(ID)

## Percent full ensemble IDs with +ve skill in observed IDs with +ve skill
sum(obs_ids %in% full_ens_ids) / length(obs_ids)    # 32%
sum(obs_ids %in% nao_matched_ids) / length(obs_ids) # 74%

## Overall statistic
## % stations with +ve MSSS
stat <-
  skill_scores %>%
  filter(period %in% aggregation_period & model %in% c("P", "PT")) %>%
  group_by(ID, subset) %>%
  filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
  ungroup() %>%
  group_by(subset) %>%
  summarize(pct_positive = sum(crpss > 0) / n() * 100)

## ####################################################### ##
## ####################################################### ##
##
## Figure 1
##
## ####################################################### ##
## ####################################################### ##

## skill_scores_subset =
##   skill_scores %>%
##   filter(model %in% c("P", "PT", "NAOPT")) %>%
##   filter(subset %in% "full", period %in% aggregation_period) %>%
##   mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
##   mutate(skill = !!sym(skill_measure))

## skill_scores_median <-
##   skill_scores_subset %>%
##   group_by(model) %>%
##   summarize(md = median(!!sym(skill_measure)), n = n()) %>%
##   mutate(md_label = sprintf(md, fmt = "%#.3f"))

## skill1 = myfun(skill_scores_subset)
## skill2 = myfun(skill_scores_subset %>% filter(!model %in% "NAOPT"))

## ## Number of stations with positive skill across models
## skill_scores_subset %>%
##   filter(model %in% c("P", "PT")) %>%
##   group_by(ID) %>%
##   filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
##   ## filter(crpss == max(crpss)) %>%
##   ungroup() %>%
##   summarize(n = sum(crpss > 0))

## obs_skill_scores %>%
##   filter(model %in% c("P", "PT")) %>%
##   group_by(ID) %>%
##   filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
##   ## filter(crpss == max(crpss)) %>%
##   ungroup() %>%
##   summarize(n = sum(crpss > 0))

## ## Number of stations with positive skill considering only model NAOPT
## skill_scores_subset %>%
##   filter(model %in% "NAOPT") %>%
##   summarize(n = sum(crpss > 0))

## ## Number of stations with positive skill not including NAOPT
## skill_scores_subset %>%
##   filter(model %in% c("P", "PT")) %>%
##   group_by(ID) %>%
##   filter(crps_ens_fcst==min(crps_ens_fcst)) %>%
##   ## filter(crpss == max(crpss)) %>%
##   ungroup() %>%
##   summarize(n = sum(crpss > 0))

## ## Median skill of all models
## median(skill1$crpss)

## ## Median skill of models P and PT
## median(skill2$crpss)

## ## Number of stations where NAOPT/PT/P is the best model
## table(skill1$model)[["NAOPT"]]
## table(skill1$model)[["PT"]]
## table(skill1$model)[["P"]]

## p1 = myplotfun1(na.omit(skill1), legend_title = toupper(skill_measure))
## p2 = myplotfun1(na.omit(skill2), legend_title = toupper(skill_measure))
## p3 = myplotfun2(skill_scores_subset, legend_title = toupper(skill_measure))

## ## Add number of stations which perform best by model
## n_naopt <- table(skill1$model)[["NAOPT"]]
## n_pt <- table(skill1$model)[["PT"]]
## n_p <- table(skill1$model)[["P"]]
## labs <- c(paste0("italic(n)==", n_p), paste0("italic(n)==", n_pt), paste0("italic(n)==", n_naopt))
## d <- data.frame(x = c(0, 0, 0), y = c(59, 58.5, 58), lab = labs, model = c("P", "PT", "NAOPT"))
## p1 <- p1 +
##   geom_point(data = d, aes(x, y, shape = model), size = 1, lwd = 0.1, show.legend = FALSE) +
##   geom_text(data = d, aes(x, y, label = lab), parse = TRUE, hjust = 0, nudge_x = 0.3, size = 2) +
##   scale_shape_manual(values = c(21, 24, 22)) +
##   theme(axis.title = element_blank())

## n_pt <- table(skill2$model)[["PT"]]
## n_p <- table(skill2$model)[["P"]]
## labs <- c(paste0("italic(n)==", n_p), paste0("italic(n)==", n_pt))
## d <- data.frame(x = c(0, 0), y = c(59, 58.5), lab = labs, model = c("P", "PT"))
## p2 <- p2 +
##   geom_point(data = d, aes(x, y, shape = model), size = 1, lwd = 0.1, show.legend = FALSE) +
##   geom_text(data = d, aes(x, y, label = lab), parse = TRUE, hjust = 0, nudge_x = 0.3, size = 2) +
##   scale_shape_manual(values = c(21, 24)) +
##   theme(axis.title = element_blank())

## p3 <- p3 + coord_fixed(ratio = 3)
## skill_scores_subset =
##   skill_scores %>%
##   filter(model %in% c("P", "PT", "NAOPT")) %>%
##   filter(subset %in% "full", period %in% aggregation_period) %>%
##   mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
##   mutate(skill = !!sym(skill_measure))

## rdbu_pal = brewer.pal(12, "RdBu")
## mylabelfun <- function(breaks) {
##   int_breaks <- as.integer(breaks * 1000)
##   idx <- sapply(int_breaks, FUN=function(x) isTRUE(all.equal(x %% 200, 0)))
##   labs <- formatC(breaks, digits=1, format="f")
##   labs[!idx] <- ""
##   labs
## }
## p1 <- p1 +
##   scale_fill_stepsn(
##     ## colours = rev(rdbu_pal)[c(4:5, 7:11)], #[3:9],
##     ## breaks = seq(-0.4, 0.8, 0.2),
##     ## limits = c(-0.3, 0.9)
##     colours = rev(rdbu_pal)[c(4:6, 7:11)], #[3:9],
##     breaks = seq(-0.2, 0.6, 0.1),
##     limits = c(-0.3, 0.7),
##     labels=mylabelfun
##   )
## p2 <- p2 +
##   scale_fill_stepsn(
##     colours = rev(rdbu_pal)[c(4:6, 7:11)], #[3:9],
##     breaks = seq(-0.2, 0.6, 0.1),
##     limits = c(-0.3, 0.7),
##     labels=mylabelfun
##     ## colours = rev(rdbu_pal)[c(4:5, 7:11)], #[3:9],
##     ## breaks = seq(-0.4, 0.8, 0.2),
##     ## limits = c(-0.3, 0.9)
##   )
## p2 <- p2 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
## p2 <- p2 + guides(shape = "none")
## skill_scores_subset =
##   skill_scores %>%
##   filter(model %in% c("P", "PT", "NAOPT")) %>%
##   filter(subset %in% "full", period %in% aggregation_period) %>%
##   mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
##   mutate(skill = !!sym(skill_measure))

## rdbu_pal = brewer.pal(12, "RdBu")
## mylabelfun <- function(breaks) {
##   int_breaks <- as.integer(breaks * 1000)
##   idx <- sapply(int_breaks, FUN=function(x) isTRUE(all.equal(x %% 200, 0)))
##   labs <- formatC(breaks, digits=1, format="f")
##   labs[!idx] <- ""
##   labs
## }
## obs_skill_scores_subset <-
##   obs_skill_scores %>%
##   filter(model %in% c("P", "PT", "NAOPT")) %>%
##   filter(period %in% obs_aggregation_period) %>%
##   mutate(period = factor(period, levels = obs_aggregation_period, labels = obs_aggregation_period_label)) %>%
##   mutate(skill = !!sym(skill_measure))

## obs_skill_scores_median <-
##   obs_skill_scores_subset %>%
##   group_by(model) %>%
##   ## summarize(md = median(msss), n = n()) %>%
##   summarize(md = median(crpss), n = n()) %>%
##   mutate(md_label = sprintf(md, fmt = "%#.3f"))

## ## Perfect prediction case
## skill1 = myfun(obs_skill_scores_subset)
## skill2 = myfun(obs_skill_scores_subset %>% filter(!model %in% "NAOPT"))

## ## Median skill of all models
## ## median(skill1$msss)
## median(skill1$crpss)
## ## Median skill of models P and PT
## ## median(skill2$msss)
## median(skill2$crpss)
## ## Number of stations where NAOPT/PT/P is the best model
## table(skill1$model)[["NAOPT"]]
## table(skill1$model)[["PT"]]
## table(skill1$model)[["P"]]

## p4 = myplotfun1(na.omit(skill1), legend_title = toupper(skill_measure))
## p5 = myplotfun1(na.omit(skill2), legend_title = toupper(skill_measure))
## p6 = myplotfun2(obs_skill_scores_subset, legend_title = toupper(skill_measure))

## n_naopt <- table(skill1$model)[["NAOPT"]]
## n_pt <- table(skill1$model)[["PT"]]
## n_p <- table(skill1$model)[["P"]]
## labs <- c(paste0("italic(n)==", n_p), paste0("italic(n)==", n_pt), paste0("italic(n)==", n_naopt))
## d <- data.frame(x = c(0, 0, 0), y = c(59, 58.5, 58), lab = labs, model = c("P", "PT", "NAOPT"))
## p4 <- p4 +
##   geom_point(data = d, aes(x, y, shape = model), size = 1, lwd = 0.1, show.legend = FALSE) +
##   geom_text(data = d, aes(x, y, label = lab), parse = TRUE, hjust = 0, nudge_x = 0.3, size = 2) +
##   scale_shape_manual(values = c(21, 24, 22)) +
##   theme(axis.title = element_blank())

## n_pt <- table(skill2$model)[["PT"]]
## n_p <- table(skill2$model)[["P"]]
## labs <- c(paste0("italic(n)==", n_p), paste0("italic(n)==", n_pt))
## d <- data.frame(x = c(0, 0), y = c(59, 58.5), lab = labs, model = c("P", "PT"))
## p5 <- p5 +
##   geom_point(data = d, aes(x, y, shape = model), size = 1, lwd = 0.1, show.legend = FALSE) +
##   geom_text(data = d, aes(x, y, label = lab), parse = TRUE, hjust = 0, nudge_x = 0.3, size = 2) +
##   scale_shape_manual(values = c(21, 24)) +
##   theme(axis.title = element_blank())

## p6 <- p6 + coord_fixed(ratio = 3)
## p4 <- p4 +
##   scale_fill_stepsn(
##     colours = rev(rdbu_pal)[c(4:5, 7:11)], #[3:9],
##     breaks = seq(-0.4, 0.8, 0.2),
##     limits = c(-0.3, 0.9)
##   )
## p5 <- p5 +
##   scale_fill_stepsn(
##     colours = rev(rdbu_pal)[c(4:5, 7:11)], #[3:9],
##     breaks = seq(-0.4, 0.8, 0.2),
##     limits = c(-0.3, 0.9)
##   )
## p5 <- p5 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
## p5 <- p5 + guides(shape = "none")

## p = p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(widths = c(2, 2, 2), ncol = 3, nrow = 2)
## p = p +
##   plot_layout(guides = "collect") &
##   theme(legend.position = "bottom",
##         legend.box = "vertical",
##         legend.justification = "left",
##         legend.box.just = "left",
##         legend.margin = margin(0, 0, 0, 0, unit = "cm"))

## p$patches$plots[[1]] =
##   p$patches$plots[[1]] +
##   labs(tag = "a") +
##   theme(plot.tag.position = c(0.145, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p$patches$plots[[2]] =
##   p$patches$plots[[2]] +
##   labs(tag = "b") +
##   theme(plot.tag.position = c(0.05, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p$patches$plots[[3]] =
##   p$patches$plots[[3]] +
##   labs(tag = "c") +
##   theme(plot.tag.position = c(0.185, 0.98),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p$patches$plots[[4]] =
##   p$patches$plots[[4]] +
##   labs(tag = "d") +
##   theme(plot.tag.position = c(0.145, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p$patches$plots[[5]] =
##   p$patches$plots[[5]] +
##   labs(tag = "e") +
##   theme(plot.tag.position = c(0.05, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p =
##   p +
##   labs(tag = "f") +
##   theme(plot.tag.position = c(0.185, 0.965),
##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))

## ## Include both modelled and observed results
## ggsave(file.path(output_dir, "fig1.png"), plot = p, width = 6, height = 7.25, units = "in")

## ## Alternative without boxplots:
## p = p1 + p2 + p4 + p5 + plot_layout(widths = c(2, 2), ncol = 2, nrow = 2)
## p = p +
##   plot_layout(guides = "collect") &
##   theme(legend.position = "bottom",
##         legend.box = "vertical",
##         legend.justification = "left",
##         legend.box.just = "left",
##         legend.margin = margin(0, 0, 0, 0, unit = "cm"))

## p$patches$plots[[1]] =
##   p$patches$plots[[1]] +
##   labs(tag = "a") +
##   theme(plot.tag.position = c(0.145, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p$patches$plots[[2]] =
##   p$patches$plots[[2]] +
##   labs(tag = "b") +
##   theme(plot.tag.position = c(0.05, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p$patches$plots[[3]] =
##   p$patches$plots[[3]] +
##   labs(tag = "c") +
##   theme(plot.tag.position = c(0.145, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p =
##   p +
##   labs(tag = "d") +
##   theme(plot.tag.position = c(0.05, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))

## ## Include both modelled and observed results
## ggsave(file.path(output_dir, "fig1_alt1.png"), plot = p, width = 6, height = 8, units = "in")

## ## Alternative without boxplots and with all models:
## p = p2 + p5 + plot_layout(widths = c(2, 2), ncol = 2, nrow = 1)
## p = p +
##   plot_layout(guides = "collect") &
##   theme(legend.position = "bottom",
##         legend.box = "vertical",
##         legend.justification = "left",
##         legend.box.just = "left",
##         legend.margin = margin(0, 0, 0, 0, unit = "cm"))

## p$patches$plots[[1]] =
##   p$patches$plots[[1]] +
##   labs(tag = "a") +
##   theme(plot.tag.position = c(0.145, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))
## p =
##   p +
##   labs(tag = "b") +
##   theme(plot.tag.position = c(0.05, 1.02),
##         plot.tag = element_text(size = tag_label_size, face="bold"))

## ## Include both modelled and observed results
## ggsave(file.path(output_dir, "fig1_alt2.png"), plot = p, width = 6, height = 6, units = "in")

## ####################################################### ##
## ####################################################### ##
##
## Figure 1
##
## ####################################################### ##
## ####################################################### ##

skill_scores_subset =
  skill_scores %>%
  filter(model %in% c("P", "PT")) %>% #, "NAOPT")) %>%
  mutate(model = droplevels(model)) %>%
  filter(subset %in% "full", period %in% aggregation_period) %>%
  mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
  mutate(skill = !!sym(skill_measure))

obs_skill_scores_subset <-
  obs_skill_scores %>%
  filter(model %in% c("P", "PT")) %>% #, "NAOPT")) %>%
  mutate(model = droplevels(model)) %>%
  filter(period %in% obs_aggregation_period) %>%
  mutate(period = factor(period, levels = obs_aggregation_period, labels = obs_aggregation_period_label)) %>%
  mutate(skill = !!sym(skill_measure))

## Schematic
example_discharge_data <-
  open_dataset(file.path(outputroot, "nrfa-discharge-summaries")) %>%
  collect() %>%
  filter(ID %in% 21017 & clim_season %in% "DJFM" & season_year %in% 1978:1991)
pp1 <- myplotfun_schematic(example_discharge_data)

## Raw ensemble
skill1 <- myfun(skill_scores_subset) # %>% filter(!model %in% "NAOPT"))
p1 <- myplotfun1(na.omit(skill1), legend_title = toupper(skill_measure))

## Observed data (i.e. perfect predictors)
skill1 <- myfun(obs_skill_scores_subset) # %>% filter(!model %in% "NAOPT"))
p2 <- myplotfun1(na.omit(skill1), legend_title = toupper(skill_measure))
p2 <- p2 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

## North Atlantic region
p3 <- myplotfun0()
p3 <- p3 + theme(plot.margin = margin(0, 0, 0, 0, unit="cm"))

## Add labels to plots
p1 <- p1 + labs(title = "c") + theme(plot.title = element_text(size = tag_label_size, face="bold"))
p2 <- p2 + labs(title = "d") + theme(plot.title = element_text(size = tag_label_size, face="bold"))
p3 <- p3 + labs(title = "e") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

design <- "
AAAABBBBCCC
AAAABBBBCCC
AAAABBBBDDD
AAAABBBBDDD
"
pp2 = p1 + p2 + p3 + guide_area() + plot_layout(design=design, guides="collect") &
  theme(
    legend.spacing.y = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0, unit="cm"),
    plot.background = element_blank()
  )

## Now use cowplot
p = plot_grid(pp1, pp2, nrow=2, align="v", rel_heights = c(1, 1.25))
ggsave(file.path(output_dir, "fig1.png"), plot = p, width = 6, height = 6.05, units = "in", dpi=fig_dpi)

## ## ####################################################### ##
## ## ####################################################### ##
## ##
## ## Figure 1: Presentation version
## ##
## ## ####################################################### ##
## ## ####################################################### ##

## skill1 = myfun(skill_scores_subset %>% filter(!model %in% "NAOPT"), skill_measure = skill_measure)
## p1 = myplotfun1(na.omit(skill1), legend_title = toupper(skill_measure))
## ## n_naopt <- table(skill1$model)[["NAOPT"]]
## n_pt <- table(skill1$model)[["PT"]]
## n_p <- table(skill1$model)[["P"]]
## ## labs <- c(paste0("italic(n)==", n_p), paste0("italic(n)==", n_pt), paste0("italic(n)==", n_naopt))
## ## d <- data.frame(x = c(0, 0, 0), y = c(59, 58.5, 58), lab = labs, model = c("P", "PT", "NAOPT"))
## labs <- c(paste0("italic(n)==", n_p), paste0("italic(n)==", n_pt))
## d <- data.frame(x = c(0, 0), y = c(59, 58.5), lab = labs, model = c("P", "PT"))
## p1 <- p1 +
##   geom_point(
##     data = d,
##     aes(x, y, shape = model),
##     size = 2,
##     lwd = 0.1,
##     show.legend = FALSE
##   ) +
##   geom_text(
##     data = d,
##     aes(x, y, label = lab),
##     parse = TRUE,
##     hjust = 0,
##     nudge_x = 0.3,
##     size = 4
##   ) +
##   scale_shape_manual(values = c(21, 24)) +
##   theme(axis.title = element_blank(),
##         axis.text.x = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         legend.position = "bottom",
##         legend.box = "vertical",
##         legend.justification = "left",
##         legend.box.just = "left",
##         legend.title = element_text(size = 11),
##         legend.text = element_text(size = 9),
##         strip.text = element_blank(),
##         panel.grid.major = element_line(size = 0.25))

## ## Perfect prediction case
## skill1 = myfun(obs_skill_scores_subset %>% filter(!model %in% "NAOPT"), skill_measure = skill_measure)
## p2 = myplotfun1(na.omit(skill1), legend_title = toupper(skill_measure))
## n_naopt <- table(skill1$model)[["NAOPT"]]
## n_pt <- table(skill1$model)[["PT"]]
## n_p <- table(skill1$model)[["P"]]
## labs <- c(paste0("italic(n)==", n_p), paste0("italic(n)==", n_pt))
## d <- data.frame(x = c(0, 0), y = c(59, 58.5), lab = labs, model = c("P", "PT"))
## p2 <- p2 +
##   geom_point(
##     data = d,
##     aes(x, y, shape = model),
##     size = 2,
##     lwd = 0.1,
##     show.legend = FALSE
##   ) +
##   geom_text(
##     data = d,
##     aes(x, y, label = lab),
##     parse = TRUE,
##     hjust = 0,
##     nudge_x = 0.3,
##     size = 4
##   ) +
##   scale_shape_manual(values = c(21, 24)) +
##   ## scale_shape_manual(values = c(21, 24, 22)) +
##   theme(axis.title = element_blank(),
##         axis.text.x = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         legend.position = "bottom",
##         legend.box = "vertical",
##         legend.justification = "left",
##         legend.box.just = "left",
##         legend.title = element_text(size = 11),
##         legend.text = element_text(size = 9),
##         strip.text = element_blank(),
##         panel.grid.major = element_line(size = 0.25))

## p1 <- p1 +
##   scale_fill_stepsn(
##     colours = rev(rdbu_pal)[c(4:5, 7:11)], #[3:9],
##     breaks = seq(-0.4, 0.8, 0.2),
##     limits = c(-0.3, 0.9)
##   )
## p2 <- p2 +
##   scale_fill_stepsn(
##     colours = rev(rdbu_pal)[c(4:5, 7:11)], #[3:9],
##     breaks = seq(-0.4, 0.8, 0.2),
##     limits = c(-0.3, 0.9)
##   )
## p2 <- p2 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
## p2 <- p2 + guides(shape = "none")

## p = p2 + p1 + plot_layout(widths = c(2, 2), ncol = 2, nrow = 1)
## p = p +
##   plot_layout(guides = "collect") &
##   theme(legend.position = "bottom",
##         legend.box = "vertical",
##         legend.justification = "left",
##         legend.box.just = "left",
##         legend.margin = margin(0, 0, 0, 0, unit = "cm"))
## ggsave(file.path(output_dir, "fig1_ppt.png"), plot = p, width = 6, height = 6, units = "in")

## ####################################################### ##
## ####################################################### ##
##
## Figure 2
##
## ####################################################### ##
## ####################################################### ##

skill_scores_subset =
  skill_scores %>%
  filter(model %in% c("P", "PT")) %>% #, "NAOPT")) %>%
  mutate(model = droplevels(model)) %>%
  filter(subset %in% c("best_n", "full")) %>%
  filter(period %in% aggregation_period) %>%
  mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
  mutate(skill = !!sym(skill_measure)) %>%
  dplyr::select(-all_of(all_skill_measures)) %>%
  pivot_wider(names_from = subset, values_from = skill) %>%
  mutate(skill_diff = best_n - full)

p4 <- myplotfun444(
  skill_scores_subset %>% filter(period %in% aggregation_period_label),
  legend_title = toupper(skill_measure)
)

skill_scores_subset <-
  skill_scores_subset %>% ungroup() %>%
  left_join(
    (skill_scores %>%
     filter(subset %in% "best_n") %>%
     filter(period %in% aggregation_period) %>%
     mutate(skill = !!sym(skill_measure)) %>%
     dplyr::select(model, crps_ens_fcst, skill, ID)),
    by = c("model", "ID")
  )

## Select best model and show compute change in skill as a result of NAO-matching
skill <-
  skill_scores_subset %>%
  group_by(ID, period) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst)) %>%
  mutate(increase = skill_diff > 0) %>%
  mutate(skill_diff = abs(skill_diff))
skill <- gauge_stns %>% left_join(skill) %>% na.omit()

p2 = myplotfun3(skill)
p2 = p2 + theme(legend.margin = margin(0, 0, 0, 0, unit = "cm"))

skill_scores_subset =
  skill_scores %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model)) %>%
  filter(subset %in% c("best_n", "full"), period %in% aggregation_period) %>%
  mutate(
    subset = factor(
      subset,
      levels = c("full", "best_n"),
      labels = c("Full ensemble", "NAO-matched ensemble")
    )
  ) %>%
  mutate(
    period = factor(
      period,
      levels = aggregation_period,
      labels = aggregation_period_label
    )
  ) %>%
  mutate(skill = !!sym(skill_measure))

p3 <- myplotfun22(
  skill_scores_subset,
  legend_title = toupper(skill_measure)
)

## Add significance [one-sided Wilcoxon signed rank test]
x <- skill_scores_subset %>%
  filter(model %in% "P" & subset %in% "Full ensemble") %>%
  `[[`(skill_measure)
y <- skill_scores_subset %>%
  filter(model %in% "P" & subset %in% "NAO-matched ensemble") %>%
  `[[`(skill_measure)
pval1 = wilcox.test(y, x, paired = TRUE, alternative = "greater")$p.value

x <- skill_scores_subset %>%
  filter(model %in% "PT" & subset %in% "Full ensemble") %>%
  `[[`(skill_measure)
y <- skill_scores_subset %>%
  filter(model %in% "PT" & subset %in% "NAO-matched ensemble") %>%
  `[[`(skill_measure)
pval2 = wilcox.test(y, x, paired = TRUE, alternative = "greater")$p.value

p3 <- p3 +
  geom_signif(
    y_position = 0.7, xmin = c(1, 3), xmax = c(2, 4),
    annotation = c(format_p(pval1), format_p(pval2)), #, format_p(pval3)),
    step_increase = 0.06,
    tip_length = 0.02,
    size = 0.25,
    textsize = 2) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = legend_label_size),
        legend.direction = "vertical",
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.box = "vertical",
        legend.box.margin = margin(-0.5, 0, 0, 0, unit = "cm"))

skill = myfun(
  skill_scores_subset %>% filter(!model %in% "NAOPT" & subset %in% "NAO-matched ensemble")
) %>% na.omit()
p1 = myplotfun5(skill, legend_title = toupper(skill_measure))
p1 <- p1 +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(.25, 'cm')
  ) +
  guides(
    fill = guide_colorbar(
      title="CRPSS",
      title.position="top",
      frame.colour = "black",
      ticks.colour = "black",
      frame.linewidth = 0.25,
      ticks.linewidth = 0.25,
      barwidth = 0.5,
      barheight = 8,
      order = 2
    ))

p2 <- p2 +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.margin=margin(0, 0, 0, 0),
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(0.25, 'cm'),
    legend.text = element_text(margin = margin(0, 0, 0, 0))
  )

design = "
AAAABBBB
AAAABBBB
AAAABBBB
AAAABBBB
CCCCDDDD
CCCCDDDD"

p1 <- p1 + labs(title = "a") + theme(plot.title = element_text(size = tag_label_size, face="bold"))
p2 <- p2 + labs(title = "b") + theme(plot.title = element_text(size = tag_label_size, face="bold"))
p3 <- p3 + labs(title = "c") + theme(plot.title = element_text(size = tag_label_size, face="bold"))
p4 <- p4 + labs(title = "d") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

p <- p1 + p2 + p3 + p4 + plot_layout(design = design)

ggsave(file.path(output_dir, "fig2.png"), plot = p, width = 6, height = 5.6, units = "in", dpi=fig_dpi)

## ## ####################################################### ##
## ## ####################################################### ##
## ##
## ## Figure 2 (presentation)
## ##
## ## ####################################################### ##
## ## ####################################################### ##

## skill_scores_subset <-
##   skill_scores %>%
##   filter(model %in% c("P", "PT", "NAOPT")) %>%
##   filter(subset %in% c("best_n", "full")) %>%
##   filter(period %in% aggregation_period) %>%
##   mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
##   mutate(skill = !!sym(skill_measure)) %>%
##   dplyr::select(-all_of(all_skill_measures)) %>%
##   pivot_wider(names_from = subset, values_from = skill) %>%
##   mutate(skill_diff = best_n - full)

## skill_scores_subset <-
##   skill_scores_subset %>% ungroup() %>%
##   filter(model %in% c("P", "PT")) %>%
##   left_join(
##     (skill_scores %>%
##      filter(subset %in% "best_n") %>%
##      filter(period %in% aggregation_period) %>%
##      mutate(skill = !!sym(skill_measure)) %>%
##      dplyr::select(model, skill, ID)),
##     by = c("model", "ID")
##   )

## skill <-
##   skill_scores_subset %>%
##   group_by(ID, period) %>%
##   filter(skill == max(skill)) %>% #msss == max(msss)) %>%
##   mutate(increase = skill_diff > 0) %>%
##   mutate(skill_diff = abs(skill_diff))
## skill <- gauge_stns %>% left_join(skill) %>% na.omit()

## p2 = myplotfun3(skill) # %>% filter(period %in% aggregation_period_label))
## p2 = p2 +
##   guides(
##     shape = guide_legend(order = 1),
##     size = guide_legend(
##       order = 2, nrow = 3,
##       byrow = FALSE,
##       override.aes = list(shape = c(rep(21, 3), rep(24, 3))),
##       direction = "vertical"
##     ),
##     fill = guide_legend(
##       order = 3,
##       override.aes = list(shape = 21, fill = c("#FC8D62", "#66C2A5"))
##     )
##   )
## p2 = p2 +
##   theme(legend.margin = margin(0, 0, 0, 0, unit = "cm"),
##         axis.title = element_blank(),
##         axis.text.x = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         ## legend.position = "bottom",
##         legend.box = "vertical",
##         legend.justification = "left",
##         legend.box.just = "left",
##         legend.title = element_text(size = 11),
##         legend.text = element_text(size = 9),
##         strip.text = element_blank(),
##         panel.grid.major = element_line(size = 0.25))

## skill_scores_subset =
##   skill_scores %>%
##   filter(model %in% c("P", "PT")) %>%
##   mutate(model = factor(model, levels = c("P", "PT"), labels = c("P", "PT"))) %>%
##   filter(subset %in% c("best_n", "full"), period %in% aggregation_period) %>%
##   mutate(
##     subset = factor(
##       subset,
##       levels = c("full", "best_n"),
##       labels = c("Full ensemble", "NAO-matched ensemble")
##     )
##   ) %>%
##   mutate(
##     period = factor(
##       period,
##       levels = aggregation_period,
##       labels = aggregation_period_label
##     )
##   ) %>%
##   mutate(skill = !!sym(skill_measure))

## skill = myfun(
##   skill_scores_subset %>% filter(!model %in% "NAOPT" & subset %in% "NAO-matched ensemble"),
##   skill_measure = skill_measure
## ) %>%
##   na.omit()

## p1 = myplotfun5(skill, legend_title = toupper(skill_measure))
## rdbu_pal = brewer.pal(11, "RdBu")
## p1 <- p1 +
##   scale_fill_stepsn(
##     colours = rev(rdbu_pal)[c(4:5, 7:11)], #[3:9],
##     breaks = seq(-0.4, 0.8, 0.2),
##     limits = c(-0.3, 0.9)
##   ) +
##   theme(legend.margin = margin(0, 0, 0, 0, unit = "cm"),
##         axis.title = element_blank(),
##         axis.text.x = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         legend.position = "right",
##         legend.box = "vertical",
##         legend.justification = "left",
##         legend.box.just = "left",
##         legend.title = element_text(size = 11),
##         legend.text = element_text(size = 9),
##         strip.text = element_blank(),
##         panel.grid.major = element_line(size = 0.25))

## p <- p1 + p2 + plot_layout(nrow = 1, ncol = 2)
## p <- p + plot_layout(heights = 1, widths = c(2, 2))

## ggsave(file.path(output_dir, "fig2_ppt.png"), plot = p, width = 8, height = 6, units = "in")

## ####################################################### ##
## ####################################################### ##
##
## Figure 3
##
## ####################################################### ##
## ####################################################### ##

skill_scores_subset =
  skill_scores %>%
  filter(subset %in% c("best_n", "full")) %>%
  filter(period %in% aggregation_period) %>%
  mutate(skill = !!sym(skill_measure)) %>%
  dplyr::select(-all_of(all_skill_measures)) %>%
  pivot_wider(names_from = subset, values_from = skill) %>%
  mutate(skill_diff = best_n - full)

skill_scores_subset <-
  skill_scores_subset %>% ungroup() %>%
  left_join(
    (skill_scores %>%
     filter(subset %in% "best_n") %>%
     filter(period %in% aggregation_period) %>%
     mutate(skill = !!sym(skill_measure)) %>%
     dplyr::select(model, ID, period, crps_ens_fcst, skill)),
    by = c("model", "ID", "period")
  )

skill =
  skill_scores_subset %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst)) %>%
  arrange(desc(best_n))

## skill =
##   skill_scores_subset %>%
##   filter(model %in% c("P", "PT")) %>%
##   group_by(ID) %>%
##   filter(skill_diff == max(skill_diff)) %>%
##   arrange(desc(best_n))

ids_best = head(skill, n=5)$ID
models_best = as.character(head(skill, n=5)$model)

predictions <-
  load_model_predictions(config, "hindcast", aggregation_period) %>%
  mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))
predictions <- predictions %>% mutate(obs = Q_95_obs, exp = Q_95_exp)

library(viridis)
p1 <- predictions %>%
  filter(ID %in% ids_best[1] & model %in% models_best[1]) %>%
  myplotfun6()
p1 <- p1 +
  labs(title = paste0("ID=", ids_best[1])) + #, " (", models_best[1], ")")) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p2 <- predictions %>%
  filter(ID %in% ids_best[2] & model %in% models_best[2]) %>%
  myplotfun6()
p2 <- p2 +
  labs(title = paste0("ID=", ids_best[2])) + #, " (", models_best[2], ")")) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p3 <- predictions %>%
  filter(ID %in% ids_best[3] & model %in% models_best[3]) %>%
  myplotfun6()
p3 <- p3 +
  labs(title = paste0("ID=", ids_best[3])) + #, " (", models_best[3], ")")) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p4 <- predictions %>%
  filter(ID %in% ids_best[4] & model %in% models_best[4]) %>%
  myplotfun6()
p4 <- p4 +
  labs(title = paste0("ID=", ids_best[4])) + #, " (", models_best[4], ")")) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p5 <- predictions %>%
  filter(ID %in% ids_best[5] & model %in% models_best[5]) %>%
  myplotfun6()
p5 <- p5 +
  labs(title = paste0("ID=", ids_best[5])) + #, " (", models_best[5], ")")) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

format_p_value <- function(p_value) {
  if (p_value < 0.01) {
    return("(P < 0.01)")
  } else {
    return(paste0("(P = ", sprintf(p_value, fmt = "%#.2f"), ")"))
  }
}

make_annotation <- function(crpss, id, mod) {
  crpss_full <- crpss %>% filter(ID %in% id & subset %in% "full" & model %in% mod) %>% `$`(crpss)
  crpss_matched <- crpss %>% filter(ID %in% id & subset %in% "best_n" & model %in% mod) %>% `$`(crpss)
  annotation = paste0(
    paste0(
      "CRPSS (Full) = ", sprintf(crpss_full, fmt = "%#.2f"), "\n"#, " ", format_p_value(acc_full$acc_p), "\n"
    ),
    paste0(
      "CRPSS (NAO-matched) = ", sprintf(crpss_matched, fmt = "%#.2f")#, "\n"#, " ", format_p_value(acc_matched$acc_p), "\n"
    )
  )
  annotation
}

get_y_range <- function(p) {
  yrange <- layer_scales(p)$y$range$range
  return(yrange)
}
get_y_position <- function(p, rel_pos) {
  yrange <- get_y_range(p)
  return(yrange[1] + diff(yrange) * rel_pos)
}

annotate_plot <- function(p, id, mod) {
  annotation_size = 4
  annotation_rel_pos = 1
  yrange <- get_y_range(p)
  yrange[2] <- yrange[2] * 1.025
  p <- p +
    ylim(yrange) +
    annotate(
      geom = "text",
      x = 1960,
      y = yrange[1] + diff(yrange) * annotation_rel_pos,
      ## label = make_annotation(acc, ids_best[1], models_best[1]),
      label = make_annotation(skill_scores, id, mod), #ids_best[1], models_best[1]),
      hjust=0,
      vjust=1,
      size = annotation_size / ggplot2::.pt
    ) +
  ## p <- p +
  ##   ylim(range) +
    annotate(
      geom = "text",
      x = 2005,
      y = yrange[1] + diff(yrange) * annotation_rel_pos,
      label = mod,
      hjust = 1,
      vjust = 1,
      size = annotation_size / ggplot2::.pt,
      fontface = "bold"
    )
  p
}
p1 <- p1 %>% annotate_plot(ids_best[1], models_best[1])
p2 <- p2 %>% annotate_plot(ids_best[2], models_best[2])
p3 <- p3 %>% annotate_plot(ids_best[3], models_best[3])
p4 <- p4 %>% annotate_plot(ids_best[4], models_best[4])
p5 <- p5 %>% annotate_plot(ids_best[5], models_best[5])

## Plot location of the selected gauge stations
gauge_stns_subset = gauge_stns %>% filter(ID %in% ids_best)
p6 = myplotfun777(gauge_stns_subset)

p1 = p1 + theme(plot.margin = margin(0, 0, 0, 0))
p2 = p2 + theme(plot.margin = margin(0, 0, 0, 0.25, unit="cm"),
                axis.title.y = element_blank())
p3 = p3 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_blank())
p4 = p4 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_text(size = axis_title_size))
p5 = p5 + theme(plot.margin = margin(0, 0, 0, 0.25, unit="cm"),
                axis.title.y = element_blank())
p6 = p6 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title = element_blank(),
                axis.text = element_text(size = axis_label_size_small))

p = p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3, nrow = 2, widths = c(2, 2, 2)) & theme(legend.position = "bottom")
p = p + plot_layout(guides = "collect") & theme(legend.margin = margin(0, 0, 0, 0), legend.key.size = unit(1, 'lines'), legend.box.spacing = unit(0, 'cm'))
ggsave(file.path(output_dir, "fig3.png"), plot = p, width = 6, height = 4.25, units = "in", dpi=fig_dpi)

## ## ####################################################### ##
## ## ####################################################### ##
## ##
## ## Figure 3 presentation
## ##
## ## ####################################################### ##
## ## ####################################################### ##

## p1 = predictions %>% filter(ID %in% ids_best[1] & model %in% models_best[1]) %>% myplotfun6()
## p2 = predictions %>% filter(ID %in% ids_best[2] & model %in% models_best[2]) %>% myplotfun6()
## p3 = predictions %>% filter(ID %in% ids_best[3] & model %in% models_best[3]) %>% myplotfun6()
## p4 = predictions %>% filter(ID %in% ids_best[4] & model %in% models_best[4]) %>% myplotfun6()
## p5 = predictions %>% filter(ID %in% ids_best[5] & model %in% models_best[5]) %>% myplotfun6()

## annotation_size = 6
## annotation_rel_pos = 1
## yrange <- get_y_range(p1)
## yrange[2] <- yrange[2] * 1.025
## p1 <- p1 +
##   ylim(yrange) +
##   annotate(
##     geom = "text",
##     x = 1960,
##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
##     label = make_annotation(acc, ids_best[1], models_best[1]),
##     hjust=0,
##     vjust=1,
##     size = annotation_size / ggplot2::.pt
##   )

## yrange <- get_y_range(p2)
## yrange[2] <- yrange[2] * 1.025
## p2 <- p2 +
##   annotate(
##     geom = "text",
##     x = 1960,
##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
##     label = make_annotation(acc, ids_best[2], models_best[2]),
##     hjust=0,
##     vjust=1,
##     size = annotation_size / ggplot2::.pt
##   )

## yrange <- get_y_range(p3)
## yrange[2] <- yrange[2] * 1.025
## p3 <- p3 +
##   annotate(
##     geom = "text",
##     x = 1960,
##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
##     label = make_annotation(acc, ids_best[3], models_best[3]),
##     hjust=0,
##     vjust=1,
##     size = annotation_size / ggplot2::.pt
##   )

## yrange <- get_y_range(p4)
## yrange[2] <- yrange[2] * 1.025
## p4 <- p4 +
##   annotate(
##     geom = "text",
##     x = 1960,
##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
##     label = make_annotation(acc, ids_best[4], models_best[4]),
##     hjust=0,
##     vjust=1,
##     size = annotation_size / ggplot2::.pt
##   )

## yrange <- get_y_range(p5)
## yrange[2] <- yrange[2] * 1.025
## p5 <- p5 +
##   annotate(
##     geom = "text",
##     x = 1960,
##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
##     label = make_annotation(acc, ids_best[5], models_best[5]),
##     hjust=0,
##     vjust=1,
##     size = annotation_size / ggplot2::.pt
##   )

## gauge_stns_subset = gauge_stns %>% filter(ID %in% ids_best)
## p6 = myplotfun777(gauge_stns_subset)
## p6 =
##   p6 +
##   ggrepel::geom_label_repel(
##              data=gauge_stns_subset,
##              aes(label=ID, geometry=geometry),
##              stat="sf_coordinates",
##              min.segment.length = 0,
##              size = 1.75)

## p1 = p1 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = 9),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = 9),
##                 axis.title.y = element_text(size = 9),
##                 axis.text.y = element_text(size = 9),
##                 axis.text.x = element_text(size = 9))
## p2 = p2 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = 9),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = 9),
##                 axis.title.y = element_blank(),
##                 axis.text.y = element_text(size = 9),
##                 axis.text.x = element_text(size = 9))
## p3 = p3 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = 9),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = 9),
##                 axis.title.y = element_blank(),
##                 axis.text.y = element_text(size = 9),
##                 axis.text.x = element_text(size = 9))
## p4 = p4 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = 9),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = 9),
##                 axis.title.y = element_text(size = 9),
##                 axis.text.y = element_text(size = 9),
##                 axis.text.x = element_text(size = 9))
## p5 = p5 + theme(panel.grid = element_blank(),
##                 strip.text = element_text(size = 9),
##                 ## legend.title = element_text(size = legend_title_size),
##                 legend.text = element_text(size = 9),
##                 axis.title.y = element_blank(),
##                 axis.text.y = element_text(size = 9),
##                 axis.text.x = element_text(size = 9))
## p6 = p6 + theme(panel.grid.major = element_line(size = 0.25),
##                 axis.title = element_blank(),
##                 axis.text = element_text(size = 9))

## p = p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, nrow = 2, widths = c(2, 2, 2)) & theme(legend.position = "bottom")
## p = p + plot_layout(guides = "collect")

## ggsave(file.path(output_dir, "fig3_ppt.png"), plot = p, width = 7, height = 6, units = "in")

## ####################################################### ##
## ####################################################### ##
##
## Figure 4
##
## ####################################################### ##
## ####################################################### ##

r2 <-
  open_dataset("results/exp1/analysis/yr2to9_lag/input") %>%
  collect() %>%
  filter(subset %in% "observed") %>%
  group_by(ID) %>%
  summarize(R2 = cor.test(Q_95, nao)$estimate ** 2)

skill_scores_subset <-
  skill_scores %>%
  filter(subset %in% c("best_n", "full")) %>%
  filter(period %in% aggregation_period) %>%
  mutate(skill = !!sym(skill_measure)) %>%
  dplyr::select(-all_of(all_skill_measures)) %>%
  pivot_wider(names_from = subset, values_from = skill) %>%
  mutate(skill_diff = best_n - full)

x <- skill_scores_subset %>%
  left_join(r2) %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model))

p <- myplotfun_scatter(x)
## p <- p +
##   guides(
##     color = guide_legend(
##       ## title = "Model",
##       ## title.position = "top",
##       ## title.hjust = 0,
##       order = 1,
##       ## direction = "horizontal"
##     ),
##     fill = guide_legend(
##       ## title = "Model",
##       ## title.position = "top",
##       ## title.hjust = 0,
##       order = 1,
##       ## direction = "horizontal"
##     )
##   )
ggsave(file.path(output_dir, "fig4.png"), plot = p, width = 5, height = 5, units = "in", dpi=fig_dpi)

## ## ####################################################### ##
## ## ####################################################### ##
## ##
## ## Figure 4 presentation
## ##
## ## ####################################################### ##
## ## ####################################################### ##

## p = myplotfun9(x, legend_title = toupper(skill_measure))
## p =
##   p + #xlim(-0.2, 0.5) +
##   scale_x_continuous(
##     name=paste0(toupper(skill_measure), " (Observations)"),
##     breaks=c(-0.2, 0, 0.2, 0.4, 0.6, 0.8),
##     limits=c(-0.2, 0.4)
##   ) +
##   geom_vline(xintercept = 0, size = 0.5) +
##   geom_hline(yintercept = 0, size = 0.5) +
##   scale_color_manual(name = "Model", values = cbbPalette) +
##   geom_smooth(method = "lm", se = FALSE) +
##   theme(
##     ## legend.position = "right",
##     legend.position = "bottom",
##     panel.grid = element_blank(),
##     axis.title.y = element_text(size = 9),
##     axis.text.y = element_text(size = 9),
##     axis.text.x = element_text(size = 9),
##     axis.title.x = element_text(size = 9),
##     legend.title = element_text(size = 9),
##     legend.text = element_text(size = 9))

## ggsave(file.path(output_dir, "fig4_ppt.png"), plot = p, width = 5, height = 5, units = "in")

## ####################################################### ##
## ####################################################### ##
##
## Figure S1
##
## ####################################################### ##
## ####################################################### ##

fit <- load_model_fit(config, "hindcast", aggregation_period) %>% mutate(kurtosis = kurtosis - 3) #%>%
  ## filter(model %in% model_levels) %>%
  ## mutate(model = factor(model, levels = model_levels, labels = model_labels))
skill_scores_subset <- skill_scores %>% dplyr::select(-aic) %>% left_join(fit)

skill_scores_subset <-
  skill_scores %>%
  dplyr::select(-aic) %>%
  filter(model %in% c("P", "PT") & subset %in% "best_n") %>%
  mutate(model = droplevels(model)) %>%
  myfun() # Selects the best performing model

skill_scores_subset <-
  skill_scores_subset %>%
  left_join(fit, by = c("ID", "model", "subset", "period"))

## Squish range to fit color scales
skill_scores_subset <-
  skill_scores_subset %>%
  mutate(
    mean = oob_squish(mean, c(-0.005, 0.005)),
    variance = oob_squish(variance, c(0.95, 1.05)),
    kurtosis = oob_squish(kurtosis, c(-2, 2)),
    skewness = oob_squish(skewness, c(-1, 1))
  )

p1 <- myplotfun_residuals(skill_scores_subset, "mean")
p1 <- p1 +
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    direction = 1,
    limits = c(-0.0050, 0.0050)
  ) +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0.25, 0, 0, unit="cm"),
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(.25, 'cm')
  ) +
  guides(
    fill = guide_colorbar(
      title="Mean",
      title.position="top",
      frame.colour = "black",
      ticks.colour = "black",
      frame.linewidth = 0.25,
      ticks.linewidth = 0.25,
      barwidth = 0.4,
      barheight = 6
    ))

p2 <- myplotfun_residuals(skill_scores_subset, "variance")
p2 <- p2 +
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    direction = 1,
    limits = c(0.95, 1.05)
  ) +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(.25, 'cm')
  ) +
  guides(
    fill = guide_colorbar(
      title="Variance",
      title.position="top",
      frame.colour = "black",
      ticks.colour = "black",
      frame.linewidth = 0.25,
      ticks.linewidth = 0.25,
      barwidth = 0.4,
      barheight = 6
    ))

p3 <- myplotfun_residuals(skill_scores_subset, "kurtosis")
p3 <- p3 +
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    direction = 1,
    limits = c(-2, 2)
    ## limits = c(-1.5, 1.5)
  ) +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0.25, 0, 0, unit="cm"),
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(.25, 'cm')
  ) +
  guides(
    fill = guide_colorbar(
      title="Kurtosis",
      title.position="top",
      frame.colour = "black",
      ticks.colour = "black",
      frame.linewidth = 0.25,
      ticks.linewidth = 0.25,
      barwidth = 0.4,
      barheight = 6
    ))

p4 <- myplotfun_residuals(skill_scores_subset, "skewness")
p4 <- p4 +
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    direction = 1,
    limits = c(-1, 1)
  ) +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(.25, 'cm')
  ) +
  guides(
    fill = guide_colorbar(
      title="Skewness",
      title.position="top",
      frame.colour = "black",
      ticks.colour = "black",
      frame.linewidth = 0.25,
      ticks.linewidth = 0.25,
      barwidth = 0.4,
      barheight = 6
    ))

p1 <- p1 + labs(title = "a") + theme(plot.title = element_text(size = tag_label_size, face="bold"))
p2 <- p2 + labs(title = "b") + theme(plot.title = element_text(size = tag_label_size, face="bold"))
p3 <- p3 + labs(title = "c") + theme(plot.title = element_text(size = tag_label_size, face="bold"))
p4 <- p4 + labs(title = "d") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

p <- p1 + p2 + p3 + p4 + plot_layout(ncol=2, nrow=2)

ggsave(file.path(output_dir, "figS1.png"), plot = p, width = 6, height = 6, units = "in", dpi=fig_dpi)

## ## ####################################################### ##
## ## ####################################################### ##
## ##
## ## Figure S1
## ##
## ## ####################################################### ##
## ## ####################################################### ##

## d <- tibble(
##   Q = c(3.5, 5, 8, 10.5),
##   Lead = 1,
##   ## Period = 8,
##   Period = 4,
##   Init = c(1979, 1980, 1981, 1982)
## ) %>%
##   mutate(Start = Init + Lead, End = Start + Period - 1)

## cbbPalette <- RColorBrewer::brewer.pal(3, "Set2")
## library(ggnewscale)
## p1 <- ggplot(d) +
##   geom_segment(
##     data = d,
##     aes(x = Start, y = Q, xend = End, yend = Q, colour = "Forecast period"),
##     size = 2,
##     alpha = .5
##   ) +
##   geom_segment(
##     data = d,
##     aes(x = Init, y = Q, xend = Start, yend = Q, colour = "Lead time"),
##     size = 2,
##     alpha = .5
##   ) +
##   scale_color_discrete(
##     name = "",
##     limits = c("Lead time", "Forecast period"),
##     guide = guide_legend(ncol = 1),
##     type = cbbPalette
##   ) +
##   new_scale_color() +
##   geom_point(
##     data = d %>%
##       gather(-Q, -Lead, -Period, key = key, value = value) %>%
##       mutate(
##         key = factor(
##           key,
##           levels = c("Init", "Start", "End"),
##           labels = c("Initialization", "Period start", "Period end")
##         )
##       ),
##     aes(x = value, y = Q, colour= key),
##     size = 3.5
##   ) +
##   scale_color_discrete(
##     name = "",
##     guide = guide_legend(ncol = 1),
##     type = cbbPalette
##   ) +
##   scale_y_continuous(
##     name = "X", #expression(bar(Y)),
##     limits = c(0, 17)
##   ) +
##   scale_x_continuous(
##     name="",
##     breaks=seq(1978, 1991, by = 2),
##     ## limits=c(1978.5, 1991.5),
##     limits=c(1978.5, 1988.5),
##     expand = c(0, 0)
##   ) +
##   theme_bw() +
##   theme(
##     panel.grid.major.x = element_line(colour = "lightgrey", size = 0.25),
##     panel.grid.minor.x = element_line(colour = "lightgrey", size = 0.25),
##     panel.grid.major.y = element_blank(),
##     panel.grid.minor.y = element_blank(),
##     axis.title.y = element_text(size = axis_title_size),
##     axis.title.x = element_blank(), #element_text(size = axis_title_size_large),
##     axis.text.y = element_blank(),
##     axis.ticks.x = element_blank(),
##     axis.ticks.y = element_blank(),
##     legend.position = "bottom",
##     plot.margin = margin(0.5, 0, 0.5, 0, unit = "cm")
##   )

## g <- ggplot_build(p)
## blue <- "#8DA0CB"
## green <- "#FC8D62"
## turquoise <- "#FC8D62"

## ## Discharge
## period_length = 4
## d <-
##   open_dataset(file.path(outputroot, "nrfa-discharge-summaries")) %>%
##   collect() %>%
##   filter(ID %in% 21017 & clim_season %in% "DJFM" & season_year %in% 1978:1991) %>%
##   mutate(Q_95_multiyear = rollapply(Q_95, period_length, mean, align = "left", fill = NA)) %>%
##   mutate(Start = season_year, End = Start + period_length - 1) %>%
##   dplyr::select(clim_season, season_year, Q_95, Q_95_multiyear, Start, End) %>%
##   mutate(group = ifelse(season_year %in% seq(1983, 1983 + period_length - 1), "a", "b")) %>%
##   mutate(Q_95_multiyear = ifelse(season_year %in% 1983, Q_95_multiyear, NA)) %>%
##   mutate(Init = Start - 1)

## p2 <- ggplot(d) +
##   geom_segment(
##     data = d,
##     aes(x = Start, y = Q_95_multiyear, xend = End, yend = Q_95_multiyear, colour = "Aggregation period"),
##     size = 2,
##     alpha = .5
##   ) +
##   scale_color_discrete(
##     name = "",
##     limits = c("Lead time", "Forecast period"),
##     guide = guide_legend(ncol = 1),
##     type = cbbPalette
##   ) +
##   scale_color_manual(
##     name = "",
##     values = turquoise,
##     limits = c("Aggregation period"),
##     guide = "none"
##   ) +
##   new_scale_color() +
##   geom_point(
##     data = d %>%
##       dplyr::select(-group, -Init) %>%
##       na.omit() %>%
##       gather(-(clim_season:Q_95_multiyear), key = key, value = value) %>%
##       mutate(
##         key = factor(
##           key,
##           levels = c("Start", "End"),
##           labels = c("Period start", "Period end")
##         )
##       ),
##     aes(x = value, y = Q_95_multiyear, colour= key),
##     size = 3.5
##   ) +
##   scale_color_manual(
##     name = "",
##     values = c(green, blue),
##     guide = "none"
##   ) +
##   new_scale_color() +
##   geom_point(
##     data = d,
##     aes(x = season_year, y = Q_95, colour = group)
##   ) +
##   scale_color_manual(
##     name = "",
##     values = c(turquoise, "darkgrey"),
##     guide = "none"
##   ) +
##   scale_x_continuous(
##     name="",
##     breaks=seq(1978, 1991, by = 2),
##     ## limits=c(1978.5, 1991.5),
##     limits=c(1978.5, 1988.5),
##     expand = c(0, 0)
##   ) +
##   scale_y_continuous(name = "Q") +
##   theme_bw() +
##   theme(
##     panel.grid.major.x = element_line(colour = "lightgrey", size = 0.25),
##     panel.grid.minor.x = element_line(colour = "lightgrey", size = 0.25),
##     panel.grid.major.y = element_blank(),
##     panel.grid.minor.y = element_blank(),
##     axis.title.y = element_text(size = axis_title_size),
##     axis.title.x = element_blank(), #element_text(size = axis_title_size_large),
##     axis.text.y = element_blank(),
##     axis.ticks.x = element_blank(),
##     axis.ticks.y = element_blank(),
##     plot.margin = margin(0.5, 0, 0.5, 0, unit = "cm")
##   )

## p <-
##   p2 + p1 +
##   plot_layout(nrow = 2, ncol = 1, guides = "collect") &
##   theme(legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 8)) #legend_label_size))
## p$patches$plots[[1]] =
##   p$patches$plots[[1]] +
##   labs(tag = "a") +
##   theme(plot.tag.position = c(0.04, 1.01),
##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))
## p = p +
##   labs(tag = "b") +
##   theme(plot.tag.position = c(0.04, 1.01),
##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))

## ggsave(file.path(output_dir, "figS1.png"), plot = p, width = 6, height = 6, units = "in")

## period_length = 4
## d <-
##   open_dataset(file.path(outputroot, "nrfa-discharge-summaries")) %>%
##   collect() %>%
##   filter(ID %in% 21017 & clim_season %in% "DJFM" & season_year %in% 1978:1991) %>%
##   mutate(Q_95_multiyear = rollapply(Q_95, period_length, mean, align = "left", fill = NA)) %>%
##   mutate(Start = season_year, End = Start + period_length - 1) %>%
##   dplyr::select(clim_season, season_year, Q_95, Q_95_multiyear, Start, End) %>%
##   mutate(group = ifelse(season_year %in% 1983:1991, "a", "b")) %>%
##   mutate(Q_95_multiyear = ifelse(season_year %in% 1983:1988, Q_95_multiyear, NA)) %>%
##   mutate(Init = Start - 1)

## p3 <- ggplot(d) +
##   geom_segment(
##     data = d,
##     aes(x = Start, y = Q_95_multiyear, xend = End, yend = Q_95_multiyear, colour = "Forecast period"),
##     size = 2,
##     alpha = .5
##   ) +
##   geom_segment(
##     data = d,
##     aes(x = Init, y = Q_95_multiyear, xend = Start, yend = Q_95_multiyear, colour = "Lead time"),
##     size = 2,
##     alpha = .5
##   ) +
##   scale_color_discrete(
##     name = "",
##     limits = c("Lead time", "Forecast period"),
##     guide = guide_legend(ncol = 1),
##     type = cbbPalette
##   ) +
##   new_scale_color() +
##   geom_point(
##     data = d %>%
##       gather(-(clim_season:Q_95_multiyear), -group, key = key, value = value) %>%
##       mutate(
##         key = factor(
##           key,
##           levels = c("Init", "Start", "End"),
##           labels = c("Initialization", "Period start", "Period end")
##         )
##       ),
##     aes(x = value, y = Q_95_multiyear, colour= key),
##     size = 3.5
##   ) +
##   scale_color_discrete(
##     name = "",
##     guide = guide_legend(ncol = 1),
##     type = cbbPalette
##   ) +
##   ## scale_color_discrete(
##   ##   name = "",
##   ##   limits = c("Lead time", "Forecast period"),
##   ##   guide = guide_legend(ncol = 1),
##   ##   type = cbbPalette
##   ## ) +
##   ## scale_color_manual(
##   ##   name = "",
##   ##   values = turquoise,
##   ##   limits = c("Forecast period"),
##   ##   guide = "none"
##   ## ) #+
##   ## new_scale_color() +
##   ## geom_point(
##   ##   data = d %>%
##   ##     dplyr::select(-group) %>%
##   ##     na.omit() %>%
##   ##     gather(-(clim_season:Q_95_multiyear), key = key, value = value) %>%
##   ##     mutate(
##   ##       key = factor(
##   ##         key,
##   ##         levels = c("Start", "End"),
##   ##         labels = c("Period start", "Period end")
##   ##       )
##   ##     ),
##   ##   aes(x = value, y = Q_95_multiyear, colour= key),
##   ##   size = 3.5
##   ## ) +
##   ## scale_color_manual(
##   ##   name = "",
##   ##   values = c(green, blue),
##   ##   guide = "none"
##   ## ) +
##   new_scale_color() +
##   geom_point(
##     data = d,
##     aes(x = season_year, y = Q_95, colour = group)
##   ) +
##   scale_color_manual(
##     name = "",
##     values = c(turquoise, "darkgrey"),
##     guide = "none"
##   ) +
##   scale_x_continuous(
##     name="",
##     breaks=seq(1978, 1991, by = 2),
##     limits=c(1978.5, 1991.5),
##     expand = c(0, 0)
##   ) +
##   scale_y_continuous(name = "Q") +
##   theme_bw() +
##   theme(legend.title = element_blank(),
##         legend.position = "bottom")

## ## p3 <- p3 +
## ##   theme(
## ##     panel.grid.major.x = element_line(colour = "lightgrey", size = 0.25),
## ##     panel.grid.minor.x = element_line(colour = "lightgrey", size = 0.25),
## ##     panel.grid.major.y = element_blank(),
## ##     panel.grid.minor.y = element_blank(),
## ##     axis.title.y = element_text(size = axis_title_size),
## ##     axis.text.x = element_text(size = axis_label_size),
## ##     axis.title.x = element_blank(), #element_text(size = axis_title_size_large),
## ##     axis.text.y = element_blank(),
## ##     axis.ticks.x = element_blank(),
## ##     axis.ticks.y = element_blank(),
## ##     legend.text = element_text(size = legend_label_size)
## ##   )

## ## ggsave(file.path(output_dir, "figS1.png"), plot = p3, width = 6, height = 4.5, units = "in")

## p3 <- p3 +
##   theme(
##     panel.grid.major.x = element_line(colour = "lightgrey", size = 0.25),
##     panel.grid.minor.x = element_line(colour = "lightgrey", size = 0.25),
##     panel.grid.major.y = element_blank(),
##     panel.grid.minor.y = element_blank(),
##     axis.title.y = element_text(size = 12),
##     axis.text.x = element_text(size = 12),
##     axis.title.x = element_blank(), #element_text(size = axis_title_size_large),
##     axis.text.y = element_blank(),
##     axis.ticks.x = element_blank(),
##     axis.ticks.y = element_blank(),
##     legend.text = element_text(size = 12)
##   )

## ggsave(file.path(output_dir, "figS1_alt1.png"), plot = p3, width = 6, height = 4.5, units = "in")

## ####################################################### ##
## ####################################################### ##
##
## Figure S2
##
## ####################################################### ##
## ####################################################### ##

skill_scores_subset =
  skill_scores %>%
  filter(
    subset %in% "full",
    period %in% aggregation_period,
    model %in% c("P", "PT")
  )

skill =
  skill_scores_subset %>%
  group_by(ID, subset, period) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst)) %>%
  arrange(desc(!!sym(skill_measure)))
ids_best = head(skill$ID, n=5)
ids_worst = tail(skill$ID, n=5)

dataset_dir = config$modelling[["hindcast"]]$input_dataset
predictions <- load_model_predictions(config, "hindcast", aggregation_period)
predictions <- predictions %>% filter(subset %in% "full")
## predictions = open_dataset(
##   file.path(outputroot, "analysis", "hindcast", "gamlss", aggregation_period, "prediction")
## ) %>%
##   collect() %>%
##   filter(model %in% c("P", "P_T", "NAO_P_T"))
## predictions$model = factor(predictions$model, levels = model_levels, labels = model_labels)
predictions <- predictions %>% mutate(obs = Q_95_obs, exp = Q_95_exp)
predictions <-
  predictions %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model))

p1 <- predictions %>%
  filter(ID %in% ids_best[1]) %>%
  myplotfun11()
p1 <- p1 +
  labs(title = paste0("ID=", ids_best[1])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p2 <- predictions %>%
  filter(ID %in% ids_best[2]) %>%
  myplotfun11()
p2 <- p2 +
  labs(title = paste0("ID=", ids_best[2])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p3 <- predictions %>%
  filter(ID %in% ids_best[3]) %>%
  myplotfun11()
p3 <- p3 +
  labs(title = paste0("ID=", ids_best[3])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p4 <- predictions %>%
  filter(ID %in% ids_best[4]) %>%
  myplotfun11()
p4 <- p4 +
  labs(title = paste0("ID=", ids_best[4])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p5 <- predictions %>%
  filter(ID %in% ids_best[5]) %>%
  myplotfun11()
p5 <- p5 +
  labs(title = paste0("ID=", ids_best[5])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p6 <- predictions %>%
  filter(ID %in% ids_worst[1]) %>%
  myplotfun11()
p6 <- p6 +
  labs(title = paste0("ID=", ids_worst[1])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p7 <- predictions %>%
  filter(ID %in% ids_worst[2]) %>%
  myplotfun11()
p7 <- p7 +
  labs(title = paste0("ID=", ids_worst[2])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p8 <- predictions %>%
  filter(ID %in% ids_worst[3]) %>%
  myplotfun11()
p8 <- p8 +
  labs(title = paste0("ID=", ids_worst[3])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p9 <- predictions %>%
  filter(ID %in% ids_worst[4]) %>%
  myplotfun11()
p9 <- p9 +
  labs(title = paste0("ID=", ids_worst[4])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p10 <- predictions %>%
  filter(ID %in% ids_worst[5]) %>%
  myplotfun11()
p10 <- p10 +
  labs(title = paste0("ID=", ids_worst[15])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

make_annotation <- function(skill_scores, id) {
  crpss_p <- skill_scores %>% filter(ID %in% id & subset %in% "full" & model %in% "P") %>% `$`(crpss)
  crpss_pt <- skill_scores %>% filter(ID %in% id & subset %in% "best_n" & model %in% "PT") %>% `$`(crpss)
  annotation = paste0(
    paste0(
      "CRPSS (P) = ", sprintf(crpss_p, fmt = "%#.2f"), "\n"#, " ", format_p_value(acc_full$acc_p), "\n"
    ),
    paste0(
      "CRPSS (PT) = ", sprintf(crpss_pt, fmt = "%#.2f")#, "\n"#, " ", format_p_value(acc_matched$acc_p), "\n"
    )
  )
  annotation
}

annotate_plot <- function(p, id) {
  annotation_size = 4
  annotation_rel_pos = 1
  yrange <- get_y_range(p)
  yrange[2] <- yrange[2] * 1.025
  p <- p +
    ylim(yrange) +
    annotate(
      geom = "text",
      x = 1960,
      y = yrange[1] + diff(yrange) * annotation_rel_pos,
      label = make_annotation(skill_scores, id),
      hjust=0,
      vjust=1,
      size = annotation_size / ggplot2::.pt
    )
  p
}

p1 <- p1 %>% annotate_plot(ids_best[1])
p2 <- p2 %>% annotate_plot(ids_best[2])
p3 <- p3 %>% annotate_plot(ids_best[3])
p4 <- p4 %>% annotate_plot(ids_best[4])
p5 <- p5 %>% annotate_plot(ids_best[5])
p6 <- p6 %>% annotate_plot(ids_worst[1])
p7 <- p7 %>% annotate_plot(ids_worst[2])
p8 <- p8 %>% annotate_plot(ids_worst[3])
p9 <- p9 %>% annotate_plot(ids_worst[4])
p10 <- p10 %>% annotate_plot(ids_worst[5])

p1 = p1 + theme(plot.margin = margin(0, 0, 0, 0))
p2 = p2 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_blank())
p3 = p3 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_blank())
p4 = p4 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_blank())
p5 = p5 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_blank())
p6 = p6 + theme(plot.margin = margin(0, 0, 0, 0))
p7 = p7 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_blank())
p8 = p8 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_blank())
p9 = p9 + theme(plot.margin = margin(0, 0, 0, 0),
                axis.title.y = element_blank())
p10 = p10 + theme(plot.margin = margin(0, 0, 0, 0),
                  axis.title.y = element_blank())

p = p1 + p2 + p3 + p4 + p6 + p7 + p8 + p9 +
  plot_layout(ncol = 4, nrow = 2, guides = "collect") &
  theme(
    legend.margin = margin(0, 0, 0, 0),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(0, 'cm')
  )

ggsave(file.path(output_dir, "figS2.png"), plot = p, width = 6, height = 4, units = "in", dpi=fig_dpi)

## ####################################################### ##
## ####################################################### ##
##
## Figure S3
##
## ####################################################### ##
## ####################################################### ##

load_fcst <- function(aggregation_period) {
  fcst = read_parquet(
    file.path(outputroot, "analysis", aggregation_period, "ensemble_mean_fcst.parquet")
  )
  fcst <-
    fcst %>%
    group_by(variable) %>%
    arrange(init_year) %>%
    mutate(
      ens_q95_lag = zoo::rollmean(ens_q95, 4, na.pad = TRUE, align = "right"),
      ens_q05_lag = zoo::rollmean(ens_q05, 4, na.pad = TRUE, align = "right")
    )

  ## Add full ensemble mean
  full_fcst <-
    fcst %>%
    dplyr::select(init_year, variable, starts_with("obs"), starts_with("ens")) %>%
    rename(full_ens_mean = ens_mean)
  full_fcst
}

compute_acc <- function(fcst, var_name, obs_name, model_name) {
  obs = fcst %>% filter(variable %in% var_name) %>% `[[`(obs_name)
  mod = fcst %>% filter(variable %in% var_name) %>% `[[`(model_name)
  na_ix = is.na(obs) | is.na(mod)
  obs = obs[!na_ix]
  mod = mod[!na_ix]
  acc = cor.test(obs, mod, method = "pearson")
  acc
}

compute_predictable_sd <- function(fcst, var_name, model_name) {
  mod = fcst %>% filter(variable %in% var_name) %>% `[[`(model_name)
  pred_sd = sd(mod, na.rm = T)
  pred_sd
}

compute_total_sd <- function(ensemble_fcst, var_name) {
  sd =
    ensemble_fcst %>%
    group_by(project, mip, source_id, member) %>%
    summarize(across(all_of(var_name), list(sd = ~sd(.x, na.rm = T)), .names = "sd"))
  tot_sd = mean(sd$sd, na.rm = T)
  tot_sd
}

compute_rpc <- function(acc, pred_sd, tot_sd) {
  rpc = acc / (pred_sd / tot_sd)
  rpc
}

format_p_value <- function(p_value) {
  if (p_value < 0.01) {
    return("(P < 0.01)")
  } else {
    return(paste0("(P = ", sprintf(p_value, fmt = "%#.2f"), ")"))
  }
}

make_annotation <- function(acc, rpc) {
  annotation = paste0(
    "ACC = ", sprintf(acc$estimate, fmt = "%#.2f"), " ", format_p_value(acc$p.value), ", ",
    "RPC = ", sprintf(rpc, fmt = "%#.1f")
  )
  annotation
}

full_fcst = load_fcst("yr2to9_lag")
obs = read_parquet(
  file.path(outputroot, "analysis", "yr2to9_lag", "obs_study_period.parquet")
)

ensemble_fcst = read_parquet(
  file.path(outputroot, "analysis", "yr2to9_lag", "ensemble_fcst.parquet")
)

nao_matched_ensemble_fcst = read_parquet(
  file.path(outputroot, "analysis", "yr2to9_lag", "matched_ensemble.parquet")
)

## Select n best performing members
nao_matched_ensemble_fcst_best =
  nao_matched_ensemble_fcst %>%
  group_by(source_id, member, init_year) %>%
  group_by(init_year, variable) %>%
  slice_min(error, n = 20)

nao_matched_fcst =
  nao_matched_ensemble_fcst_best %>%
  group_by(init_year, variable) %>%
  summarize(ens_mean = mean(value, na.rm=TRUE)) %>%
  ungroup()

## Join with observed data
nao_matched_fcst =
  nao_matched_fcst %>%
  left_join(obs, by = c("init_year", "variable"))

## Smooth
nao_matched_fcst =
  nao_matched_fcst %>%
  group_by(variable) %>%
  mutate(
    ens_mean_lag = rollmean(
      ens_mean,
      4,
      na.pad=TRUE,
      align="right"
    )
  )

## Adjust variance to match that of observed
nao_matched_fcst =
  nao_matched_fcst %>%
  group_by(variable) %>%
  mutate(ens_mean_var_adj = ens_mean * sd(obs) / sd(ens_mean, na.rm=T)) %>%
  mutate(ens_mean_lag_std = ens_mean_lag / sd(ens_mean_lag, na.rm=T)) %>%
  mutate(ens_mean_lag_var_adj = ens_mean_lag_std * sd(obs))

## Plot 1 [analog of Fig 2a from Smith et al. 2020]
p1 <- myplotfun_nao_raw(full_fcst, ensemble_fcst)
p1 <- p1 + labs(title = "a") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 2 [analog of Fig 2b from Smith et al. 2020]
p2 <- myplotfun_nao_matched(full_fcst, ensemble_fcst)
p2 <- p2 + labs(title = "b") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 3 [analog of Fig 2c from Smith et al. 2020]
p3 <- myplotfun_amv_raw(full_fcst, ensemble_fcst)
p3 <- p3 + labs(title = "c") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 4 [analog of Fig 2d from Smith et al. 2020]
p4 <- myplotfun_amv_matched(nao_matched_fcst, ensemble_fcst)
p4 <- p4 + labs(title = "d") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 5 [analog of Fig 2e from Smith et al. 2020]
p5 <- myplotfun_precip_raw(full_fcst, ensemble_fcst)
p5 <- p5 + labs(title = "e") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 6 [analog of Fig 2f from Smith et al. 2020]
p6 <- myplotfun_precip_matched(nao_matched_fcst, ensemble_fcst)
p6 <- p6 + labs(title = "f") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

## ## Plot 7
## p7 <- myplotfun_temp_raw(full_fcst, ensemble_fcst)
## p7 <- p7 + labs(title = "e") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

## ## Plot 8
## p8 <- myplotfun_temp_matched(nao_matched_fcst, ensemble_fcst)
## p8 <- p8 + labs(title = "f") + theme(plot.title = element_text(size = tag_label_size, face="bold"))

p1 = p1 + theme(axis.title.y = element_text(size = axis_title_size_small))
p3 = p3 + theme(axis.title.y = element_text(size = axis_title_size_small))
p5 = p5 + theme(axis.title.y = element_text(size = axis_title_size_small))
## p7 = p7 + theme(axis.title.y = element_text(size = axis_title_size_small))

p2 = p2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
p4 = p4 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
p6 = p6 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
## p8 = p8 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())

p = p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2, nrow = 3) & theme(legend.position = "bottom")
## p = p1 + p2 + p3 + p4 + p7 + p8 + p5 + p6 + plot_layout(ncol = 2, nrow = 4) & theme(legend.position = "bottom")
p = p + plot_layout(guides = "collect")

ggsave(file.path(output_dir, "figS3.png"), p, width = 6, height = 6, units = "in", dpi=fig_dpi)

## ## ####################################################### ##
## ## ####################################################### ##
## ##
## ## Figure S3 presentation
## ##
## ## ####################################################### ##
## ## ####################################################### ##

## full_fcst = load_fcst("yr2to9_lag")
## obs = read_parquet(
##   file.path(outputroot, "analysis", "yr2to9_lag", "obs_study_period.parquet")
## )

## ensemble_fcst = read_parquet(
##   file.path(outputroot, "analysis", "yr2to9_lag", "ensemble_fcst.parquet")
## )

## nao_matched_ensemble_fcst = read_parquet(
##   file.path(outputroot, "analysis", "yr2to9_lag", "matched_ensemble.parquet")
## )

## ## Select n best performing members
## nao_matched_ensemble_fcst_best =
##   nao_matched_ensemble_fcst %>%
##   group_by(source_id, member, init_year) %>%
##   group_by(init_year, variable) %>%
##   slice_min(error, n = 20)

## nao_matched_fcst =
##   nao_matched_ensemble_fcst_best %>%
##   group_by(init_year, variable) %>%
##   summarize(ens_mean = mean(value, na.rm=TRUE)) %>%
##   ungroup()

## ## Join with observed data
## nao_matched_fcst =
##   nao_matched_fcst %>%
##   left_join(obs, by = c("init_year", "variable"))

## ## Smooth
## nao_matched_fcst =
##   nao_matched_fcst %>%
##   group_by(variable) %>%
##   mutate(
##     ens_mean_lag = rollmean(
##       ens_mean,
##       4,
##       na.pad=TRUE,
##       align="right"
##     )
##   )

## ## Adjust variance to match that of observed
## nao_matched_fcst =
##   nao_matched_fcst %>%
##   group_by(variable) %>%
##   mutate(ens_mean_var_adj = ens_mean * sd(obs) / sd(ens_mean, na.rm=T)) %>%
##   mutate(ens_mean_lag_std = ens_mean_lag / sd(ens_mean_lag, na.rm=T)) %>%
##   mutate(ens_mean_lag_var_adj = ens_mean_lag_std * sd(obs))

## ## Plot 1 [analog of Fig 2a from Smith et al. 2020]
## acc = compute_acc(full_fcst, "nao", "obs", "full_ens_mean")
## pred_sd = compute_predictable_sd(full_fcst, "nao", "full_ens_mean")
## tot_sd = compute_total_sd(ensemble_fcst, "nao")
## rpc = compute_rpc(acc$estimate, pred_sd, tot_sd)

## plotdata =
##   full_fcst %>%
##   filter(variable %in% "nao") %>%
##   pivot_longer(c(-init_year, -variable), names_to = "statistic", values_to = "value") %>%
##   filter(statistic %in% c("obs", "full_ens_mean", "ens_q95", "ens_q05")) %>%
##   mutate(statistic = factor(
##            statistic,
##            levels = c("obs", "full_ens_mean", "ens_q95", "ens_q05"),
##            labels = c("Observed", "Modelled", "ens_q95", "ens_q05")))

## p1 = ggplot() +
##   geom_ribbon(
##     data = plotdata %>% filter(statistic %in% c("ens_q95", "ens_q05")) %>% pivot_wider(names_from = statistic, values_from = value),
##     aes(x = init_year, ymin = ens_q05, ymax = ens_q95), fill = "red", alpha=0.15
##   ) +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
##     aes(x = init_year, y = value, color = statistic),
##     size = 1
##   ) +
##   geom_hline(yintercept=0, size=0.25) +
##   scale_y_continuous(
##     name="NAO anomaly (hPa)",
##     breaks=seq(-7.5, 7.5, by=2.5),
##     limits=c(-7.5, 7.5)
##   ) +
##   scale_x_continuous(
##     name = "",
##     breaks = seq(1960, 2000, 10),
##     limits = c(1960, 2005)
##   ) +
##   scale_color_discrete(
##     name = "",
##     labels = c("Observed", "Modelled"),
##     type = cbbPalette[2:1]
##   ) +
##   theme_bw() +
##   theme(panel.grid = element_blank(),
##         strip.text = element_text(size = 9),
##         legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 9),
##         axis.title.y = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         axis.text.x = element_text(size = 9))

## ## p1 =
## ##   p1 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = 7.5,
## ##     label = make_annotation(acc, rpc),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = 9 / ggplot2::.pt
## ##   ) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = -7.5,
## ##     label = "Raw ensemble",
## ##     hjust=0,
## ##     vjust=0,
## ##     size = 9 / ggplot2::.pt
## ##   )

## ## Plot 2 [analog of Fig 2b from Smith et al. 2020]
## acc = compute_acc(full_fcst, "nao", "obs", "ens_mean_lag_var_adj")
## pred_sd = compute_predictable_sd(full_fcst, "nao", "ens_mean_lag")
## tot_sd = compute_total_sd(ensemble_fcst, "nao")
## rpc = compute_rpc(acc$estimate, pred_sd, tot_sd)

## plotdata =
##   full_fcst %>%
##   filter(variable %in% "nao") %>%
##   pivot_longer(c(-init_year, -variable), names_to = "statistic", values_to = "value") %>%
##   filter(statistic %in% c("obs", "ens_mean_lag_var_adj", "ens_mean_var_adj")) %>%
##   mutate(statistic = factor(
##            statistic,
##            levels = c("obs", "ens_mean_lag_var_adj", "ens_mean_var_adj") ,
##            labels = c("Observed", "Modelled", "ens_mean_var_adj"))) %>%
##   mutate(value = ifelse(init_year < 1964, NA, value))

## p2 = ggplot() +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("ens_mean_var_adj")),
##     aes(x = init_year, y = value), color = "#F8766D", size = 0.25
##   ) +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
##     aes(x = init_year, y = value, color = statistic),
##     size = 1
##   ) +
##   geom_hline(yintercept=0, size=0.25) +
##   scale_y_continuous(
##     name="NAO anomaly (hPa)",
##     breaks=seq(-7.5, 7.5, by=2.5),
##     limits=c(-7.5, 7.5)
##   ) +
##   scale_x_continuous(
##     name = "",
##     breaks = seq(1960, 2000, 10),
##     limits = c(1960, 2005)
##   ) +
##   scale_color_discrete(
##     name = "",
##     labels = c("Observed", "Modelled"),
##     type = cbbPalette[2:1]
##   ) +
##   theme_bw() +
##   theme(panel.grid = element_blank(),
##         strip.text = element_text(size = 9),
##         legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 9),
##         axis.title.y = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         axis.text.x = element_text(size = 9))

## ## p2 =
## ##   p2 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = 7.5,
## ##     label = make_annotation(acc, rpc),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = 9 / ggplot2::.pt
## ##   ) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = -7.5,
## ##     label = "Variance-adjusted and lagged",
## ##     hjust=0,
## ##     vjust=0,
## ##     size = 9 / ggplot2::.pt
## ##   )

## ## Plot 3 [analog of Fig 2c from Smith et al. 2020]
## acc = compute_acc(full_fcst, "amv", "obs", "ens_mean_lag")
## pred_sd = compute_predictable_sd(full_fcst, "amv", "ens_mean_lag")
## tot_sd = compute_total_sd(ensemble_fcst, "amv")
## rpc = compute_rpc(acc$estimate, pred_sd, tot_sd)

## plotdata =
##   full_fcst %>%
##   filter(variable %in% "amv") %>%
##   pivot_longer(c(-init_year, -variable), names_to = "statistic", values_to = "value") %>%
##   filter(statistic %in% c("obs", "ens_mean_lag", "ens_q95_lag", "ens_q05_lag")) %>%
##   mutate(statistic = factor(
##            statistic,
##            levels = c("obs", "ens_mean_lag", "ens_q95_lag", "ens_q05_lag"),
##            labels = c("Observed", "Modelled", "ens_q95", "ens_q05"))) %>%
##   mutate(value = ifelse(init_year < 1964, NA, value))

## p3 = ggplot() +
##   geom_ribbon(
##     data = plotdata %>% filter(statistic %in% c("ens_q95", "ens_q05")) %>% pivot_wider(names_from = statistic, values_from = value),
##     aes(x = init_year, ymin = ens_q05, ymax = ens_q95), fill = "red", alpha=0.15
##   ) +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
##     aes(x = init_year, y = value, color = statistic),
##     size = 1
##   ) +
##   geom_hline(yintercept=0, size=0.25) +
##   scale_y_continuous(
##     name="AMV anomaly (K)",
##     breaks=c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3),
##     limits=c(-0.3, 0.3)
##   ) +
##   scale_x_continuous(
##     name = "",
##     breaks = seq(1960, 2000, 10),
##     limits = c(1960, 2005)
##   ) +
##   scale_color_discrete(
##     name = "",
##     labels = c("Observed", "Modelled"),
##     type = cbbPalette[2:1]
##   ) +
##   theme_bw() +
##   theme(panel.grid = element_blank(),
##         strip.text = element_text(size = 9),
##         legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 9),
##         axis.title.y = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         axis.text.x = element_text(size = 9))

## ## p3 =
## ##   p3 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = 0.3,
## ##     label = make_annotation(acc, rpc),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = 9 / ggplot2::.pt
## ##   ) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = -0.3,
## ##     label = "Raw lagged ensemble",
## ##     hjust=0,
## ##     vjust=0,
## ##     size = 9 / ggplot2::.pt
## ##   )

## ## Plot 4 [analog of Fig 2d from Smith et al. 2020]
## acc = compute_acc(nao_matched_fcst, "amv", "obs", "ens_mean")
## pred_sd = compute_predictable_sd(nao_matched_fcst, "amv", "ens_mean")
## tot_sd = compute_total_sd(ensemble_fcst, "amv")
## rpc = compute_rpc(acc$estimate, pred_sd, tot_sd)

## plotdata =
##   nao_matched_fcst %>%
##   filter(variable %in% "amv") %>%
##   pivot_longer(c(-init_year, -variable), names_to = "statistic", values_to = "value") %>%
##   filter(statistic %in% c("obs", "ens_mean_var_adj")) %>%
##   mutate(statistic = factor(statistic, levels = c("obs", "ens_mean_var_adj"), labels = c("Observed", "Modelled"))) %>%
##   mutate(value = ifelse(init_year < 1964, NA, value))

## p4 = ggplot() +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
##     aes(x = init_year, y = value, color = statistic),
##     size = 1
##   ) +
##   geom_hline(yintercept=0, size=0.25) +
##   scale_y_continuous(
##     name="AMV anomaly (K)",
##     breaks=c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3),
##     limits=c(-0.3, 0.3)
##   ) +
##   scale_x_continuous(
##     name = "",
##     breaks = seq(1960, 2000, 10),
##     limits = c(1960, 2005)
##   ) +
##   scale_color_discrete(
##     name = "",
##     labels = c("Observed", "Modelled"),
##     type = cbbPalette[2:1]
##   ) +
##   theme_bw() +
##   theme(panel.grid = element_blank(),
##         strip.text = element_text(size = 9),
##         legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 9),
##         axis.title.y = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         axis.text.x = element_text(size = 9))

## ## p4 =
## ##   p4 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = 0.3,
## ##     label = make_annotation(acc, rpc),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = 9 / ggplot2::.pt
## ##   ) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = -0.3,
## ##     label = "NAO-matched",
## ##     hjust=0,
## ##     vjust=0,
## ##     size = 9 / ggplot2::.pt
## ##   )

## ## Plot 5 [analog of Fig 2e from Smith et al. 2020]
## acc = compute_acc(full_fcst, "european_precip", "obs", "ens_mean_lag")
## pred_sd = compute_predictable_sd(full_fcst, "european_precip", "ens_mean_lag")
## tot_sd = compute_total_sd(ensemble_fcst, "european_precip")
## rpc = compute_rpc(acc$estimate, pred_sd, tot_sd)

## plotdata =
##   full_fcst %>%
##   filter(variable %in% "european_precip") %>%
##   pivot_longer(c(-init_year, -variable), names_to = "statistic", values_to = "value") %>%
##   filter(statistic %in% c("obs", "ens_mean_lag", "ens_q95_lag", "ens_q05_lag")) %>%
##   mutate(statistic = factor(
##            statistic,
##            levels = c("obs", "ens_mean_lag", "ens_q95_lag", "ens_q05_lag"),
##            labels = c("Observed", "Modelled", "ens_q95", "ens_q05"))) %>%
##   mutate(value = ifelse(init_year < 1964, NA, value))

## p5 = ggplot() +
##   geom_ribbon(
##     data = plotdata %>% filter(statistic %in% c("ens_q95", "ens_q05")) %>% pivot_wider(names_from = statistic, values_from = value),
##     aes(x = init_year, ymin = ens_q05, ymax = ens_q95), fill = "red", alpha=0.15
##   ) +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
##     aes(x = init_year, y = value, color = statistic),
##     size = 1
##   ) +
##   geom_hline(yintercept=0, size=0.25) +
##   scale_y_continuous(
##     name=expression(atop("Northern European", paste("precipitation anomaly " (mm~day^{-1})))),
##     breaks=c(-0.5, -0.25, 0.0, 0.25, 0.5),
##     limits=c(-0.5, 0.5)
##   ) +
##   scale_x_continuous(
##     name = "",
##     breaks = seq(1960, 2000, 10),
##     limits = c(1960, 2005)
##   ) +
##   scale_color_discrete(
##     name = "",
##     labels = c("Observed", "Modelled"),
##     type = cbbPalette[2:1]
##   ) +
##   theme_bw() +
##   theme(panel.grid = element_blank(),
##         strip.text = element_text(size = 9),
##         legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 9),
##         axis.title.y = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         axis.text.x = element_text(size = 9))

## ## p5 =
## ##   p5 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = 0.5,
## ##     label = make_annotation(acc, rpc),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = 9 / ggplot2::.pt
## ##   ) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = -0.5,
## ##     label = "Raw lagged ensemble",
## ##     hjust=0,
## ##     vjust=0,
## ##     size = 9 / ggplot2::.pt
## ##   )

## ## Plot 6 [analog of Fig 2f from Smith et al. 2020]
## acc = compute_acc(nao_matched_fcst, "european_precip", "obs", "ens_mean")
## pred_sd = compute_predictable_sd(nao_matched_fcst, "european_precip", "ens_mean")
## tot_sd = compute_total_sd(ensemble_fcst, "european_precip")
## rpc = compute_rpc(acc$estimate, pred_sd, tot_sd)

## plotdata =
##   nao_matched_fcst %>%
##   filter(variable %in% "european_precip") %>%
##   pivot_longer(c(-init_year, -variable), names_to = "statistic", values_to = "value") %>%
##   filter(statistic %in% c("obs", "ens_mean_var_adj")) %>%
##   mutate(statistic = factor(statistic, levels = c("obs", "ens_mean_var_adj"), labels = c("Observed", "Modelled"))) %>%
##   mutate(value = ifelse(init_year < 1964, NA, value))

## p6 = ggplot() +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
##     aes(x = init_year, y = value, color = statistic),
##     size = 1
##   ) +
##   geom_hline(yintercept=0, size=0.25) +
##   scale_y_continuous(
##     name=expression(atop("Northern European", paste("precipitation anomaly " (mm~day^{-1})))),
##     breaks=c(-0.5, -0.25, 0.0, 0.25, 0.5),
##     limits=c(-0.5, 0.5)
##   ) +
##   scale_x_continuous(
##     name = "",
##     breaks = seq(1960, 2000, 10),
##     limits = c(1960, 2005)
##   ) +
##   scale_color_discrete(
##     name = "",
##     labels = c("Observed", "Modelled"),
##     type = cbbPalette[2:1]
##   ) +
##   theme_bw() +
##   theme(panel.grid = element_blank(),
##         strip.text = element_text(size = 9),
##         legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 9),
##         axis.title.y = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         axis.text.x = element_text(size = 9))

## ## p6 =
## ##   p6 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = 0.5,
## ##     label = make_annotation(acc, rpc),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = 9 / ggplot2::.pt
## ##   ) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = -0.5,
## ##     label = "NAO-matched",
## ##     hjust=0,
## ##     vjust=0,
## ##     size = 9 / ggplot2::.pt
## ##   )

## ## Plot 7
## acc = compute_acc(full_fcst, "uk_temp", "obs", "ens_mean_lag")
## pred_sd = compute_predictable_sd(full_fcst, "uk_temp", "ens_mean_lag")
## tot_sd = compute_total_sd(ensemble_fcst, "uk_temp")
## rpc = compute_rpc(acc$estimate, pred_sd, tot_sd)

## plotdata =
##   full_fcst %>%
##   filter(variable %in% "uk_temp") %>%
##   pivot_longer(c(-init_year, -variable), names_to = "statistic", values_to = "value") %>%
##   filter(statistic %in% c("obs", "ens_mean_lag", "ens_q95_lag", "ens_q05_lag")) %>%
##   mutate(statistic = factor(
##            statistic,
##            levels = c("obs", "ens_mean_lag", "ens_q95_lag", "ens_q05_lag"),
##            labels = c("Observed", "Modelled", "ens_q95", "ens_q05"))) %>%
##   mutate(value = ifelse(init_year < 1964, NA, value))

## p7 = ggplot() +
##   geom_ribbon(
##     data = plotdata %>% filter(statistic %in% c("ens_q95", "ens_q05")) %>% pivot_wider(names_from = statistic, values_from = value),
##     aes(x = init_year, ymin = ens_q05, ymax = ens_q95), fill = "red", alpha=0.15
##   ) +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
##     aes(x = init_year, y = value, color = statistic),
##     size = 1
##   ) +
##   geom_hline(yintercept=0, size=0.25) +
##   scale_y_continuous(
##     name=expression(atop("Northern European", "temperature anomaly (K)")),
##     breaks=c(-1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2),
##     limits=c(-1.2, 1.2)
##   ) +
##   scale_x_continuous(
##     name = "",
##     breaks = seq(1960, 2000, 10),
##     limits = c(1960, 2005)
##   ) +
##   scale_color_discrete(
##     name = "",
##     labels = c("Observed", "Modelled"),
##     type = cbbPalette[2:1]
##   ) +
##   theme_bw() +
##   theme(panel.grid = element_blank(),
##         strip.text = element_text(size = 9),
##         legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 9),
##         axis.title.y = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         axis.text.x = element_text(size = 9))

## ## p7 =
## ##   p7 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = 1.2,
## ##     label = make_annotation(acc, rpc),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = 9 / ggplot2::.pt
## ##   ) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = -1.2,
## ##     label = "Raw lagged ensemble",
## ##     hjust=0,
## ##     vjust=0,
## ##     size = 9 / ggplot2::.pt
## ##   )

## ## Plot 8
## acc = compute_acc(nao_matched_fcst, "uk_temp", "obs", "ens_mean")
## pred_sd = compute_predictable_sd(nao_matched_fcst, "uk_temp", "ens_mean")
## tot_sd = compute_total_sd(ensemble_fcst, "uk_temp")
## rpc = compute_rpc(acc$estimate, pred_sd, tot_sd)

## plotdata =
##   nao_matched_fcst %>%
##   filter(variable %in% "uk_temp") %>%
##   pivot_longer(c(-init_year, -variable), names_to = "statistic", values_to = "value") %>%
##   filter(statistic %in% c("obs", "ens_mean_var_adj")) %>%
##   mutate(statistic = factor(statistic, levels = c("obs", "ens_mean_var_adj"), labels = c("Observed", "Modelled"))) %>%
##   mutate(value = ifelse(init_year < 1964, NA, value))

## p8 = ggplot() +
##   geom_line(
##     data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
##     aes(x = init_year, y = value, color = statistic),
##     size = 1
##   ) +
##   geom_hline(yintercept=0, size=0.25) +
##   scale_y_continuous(
##     name=expression(atop("Northern European", "temperature anomaly (K)")),
##     breaks=c(-1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2),
##     limits=c(-1.2, 1.2)
##   ) +
##   scale_x_continuous(
##     name = "",
##     breaks = seq(1960, 2000, 10),
##     limits = c(1960, 2005)
##   ) +
##   scale_color_discrete(
##     name = "",
##     labels = c("Observed", "Modelled"),
##     type = cbbPalette[2:1]
##   ) +
##   theme_bw() +
##   theme(panel.grid = element_blank(),
##         strip.text = element_text(size = 9),
##         legend.title = element_blank(),
##         legend.position = "bottom",
##         legend.text = element_text(size = 9),
##         axis.title.y = element_text(size = 9),
##         axis.text.y = element_text(size = 9),
##         axis.text.x = element_text(size = 9))

## ## p8 =
## ##   p8 +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = 1.2,
## ##     label = make_annotation(acc, rpc),
## ##     hjust=0,
## ##     vjust=1,
## ##     size = 9 / ggplot2::.pt
## ##   ) +
## ##   annotate(
## ##     geom = "text",
## ##     x = 1960, y = -1.2,
## ##     label = "NAO-matched",
## ##     hjust=0,
## ##     vjust=0,
## ##     size = 9 / ggplot2::.pt
## ##   )

## p1 = p1 + theme(axis.title.y = element_text(size = 9))
## p3 = p3 + theme(axis.title.y = element_text(size = 9))
## p5 = p5 + theme(axis.title.y = element_text(size = 9))
## p7 = p7 + theme(axis.title.y = element_text(size = 9))

## p2 = p2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
## p4 = p4 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
## p6 = p6 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
## p8 = p8 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())

## ## p = p1 + p2 + p3 + p4 + p7 + p8 + p5 + p6 + plot_layout(ncol = 2, nrow = 4) & theme(legend.position = "bottom")
## p = p1 + p2 + p7 + p8 + p5 + p6 + plot_layout(ncol = 2, nrow = 3) & theme(legend.position = "bottom")
## p = p +
##   plot_layout(guides = "collect")
##   theme(legend.position = "bottom",
##         ## legend.box = "vertical",
##         ## legend.justification = "left",
##         ## legend.box.just = "left",
##         legend.margin = margin(0, 0, 0, 0, unit = "cm"))

## ## p$patches$plots[[1]] =
## ##   p$patches$plots[[1]] +
## ##   labs(tag = "a") +
## ##   theme(plot.tag.position = c(0.215, 1),
## ##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))
## ## p$patches$plots[[2]] =
## ##   p$patches$plots[[2]] +
## ##   labs(tag = "b") +
## ##   theme(plot.tag.position = c(0.03, 1),
## ##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))
## ## p$patches$plots[[3]] =
## ##   p$patches$plots[[3]] +
## ##   labs(tag = "c") +
## ##   theme(plot.tag.position = c(0.215, 1),
## ##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))
## ## p$patches$plots[[4]] =
## ##   p$patches$plots[[4]] +
## ##   labs(tag = "d") +
## ##   theme(plot.tag.position = c(0.03, 1),
## ##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))
## ## p$patches$plots[[5]] =
## ##   p$patches$plots[[5]] +
## ##   labs(tag = "e") +
## ##   theme(plot.tag.position = c(0.215, 1),
## ##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))
## ## p$patches$plots[[6]] =
## ##   p$patches$plots[[6]] +
## ##   labs(tag = "f") +
## ##   theme(plot.tag.position = c(0.03, 1),
## ##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))
## ## p$patches$plots[[7]] =
## ##   p$patches$plots[[7]] +
## ##   labs(tag = "g") +
## ##   theme(plot.tag.position = c(0.215, 1),
## ##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))
## ## p =
## ##   p +
## ##   labs(tag = "h") +
## ##   theme(plot.tag.position = c(0.03, 1),
## ##         plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))

## ggsave(file.path(output_dir, "figS3_ppt.png"), p, width = 6, height = 7, units = "in")

## ####################################################### ##
## ####################################################### ##
##
## Figure S4
##
## ####################################################### ##
## ####################################################### ##

dataset_dir = config$modelling[["hindcast"]]$input_dataset

skill_scores_subset <-
  skill_scores %>%
  filter(model %in% c("P", "PT") & subset %in% c("full", "best_n")) %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID, subset, period) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst)) %>%
  dplyr::select(-any_of(all_skill_measures)) %>% ungroup()

predictions <-
  load_model_predictions(config, "hindcast", aggregation_period) %>%
  mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))

predictions <- predictions %>%
  dplyr::select(-matches("Q[0-9]{2}")) %>%
  dplyr::select(-mu, -sigma, -shape, -scale, -date) %>%
  dplyr::select(-any_of(all_skill_measures)) #%>%
  ## mutate(model = factor(model, levels = model_levels, labels = model_labels))
## predictions = open_dataset(
##   file.path(outputroot, "analysis", "hindcast", "gamlss", aggregation_period, "prediction")
## ) %>%
##   collect() %>%
##   dplyr::select(-matches("Q[0-9]{2}")) %>%
##   dplyr::select(-mu, -sigma, -shape, -scale, -date) %>%
##   dplyr::select(-any_of(all_skill_measures)) %>%
##   mutate(model = factor(model, levels = model_levels, labels = model_labels))

## Take the best model for each ID/subset combination
id_year <- predictions %>% dplyr::select(ID, year) %>% distinct()

skill_scores_subset <-
  id_year %>%
  full_join(skill_scores_subset) %>% dplyr::select(-period)

skill_scores_subset <-
  skill_scores_subset %>%
  mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))

predictions <- skill_scores_subset %>% left_join(predictions, by = c("ID", "subset", "model", "year"))
## predictions <-
##   predictions %>%
##   ## mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))

## predictions$model = factor(predictions$model, levels = model_levels, labels = model_labels)
predictions <- predictions %>% mutate(obs = Q_95_obs, exp = Q_95_exp)

## Compute anomalies
predictions0 <-
  predictions %>%
  group_by(ID, subset) %>% #period, ID, predictand, subset) %>%
  ## group_by(period, ID, predictand, subset) %>%
  summarize(
    ens_mean = mean(exp, na.rm = TRUE),
    obs_mean = mean(obs, na.rm = TRUE),
    ens_sd = sd(exp, na.rm = TRUE),
    obs_sd = sd(obs, na.rm = TRUE)
  )

predictions1 <-
  predictions %>%
  filter(year %in% 1985:1993) %>%
  group_by(ID, subset) %>% #period, ID, predictand, subset) %>%
  ## group_by(period, ID, predictand, subset) %>%
  summarize(
    ens_subset_mean = mean(exp, na.rm = TRUE),
    obs_subset_mean = mean(obs, na.rm = TRUE),
    ens_subset_sd = sd(exp, na.rm = TRUE),
    obs_subset_sd = sd(obs, na.rm = TRUE)
  )

predictions0 <- predictions0 %>% left_join(predictions1)

predictions0 <-
  predictions0 %>%
  mutate(
    obs_anom = (obs_subset_mean - obs_mean) / obs_sd,
    ens_anom = (ens_subset_mean - ens_mean) / ens_sd
  )

dat1 <- predictions0 %>% filter(subset %in% "Full ensemble") %>% mutate(anom = obs_anom)
dat1 <- gauge_stns %>% left_join(dat1, by = "ID")
p1 <- myplotfun888(dat1)
p1 <- p1 +
  theme(
    axis.title.y = element_text(size = axis_title_size, margin = margin(0, -0.5, 0, 0, unit="cm")),
    plot.margin = margin(0, 0.25, 0, 0, unit = "cm")
  )
p1 <- p1 +
  labs(title = "a. Observed") +
  theme(plot.title = element_text(size = tag_label_size, face = "bold", margin = margin(0,0,1,0, unit="pt")))

dat2 <- predictions0 %>% filter(subset %in% "Full ensemble") %>% mutate(anom = ens_anom)
dat2 <- gauge_stns %>% left_join(dat2, by = "ID")
p2 <- myplotfun888(dat2)
p2 <- p2 +
  theme(
    axis.title.y = element_text(size = axis_title_size, margin = margin(0, 0, 0, 0, unit="cm")),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 0.25, 0, 0, unit = "cm")
  )
p2 <- p2 +
  labs(title = "b. Raw ensemble") +
  theme(plot.title = element_text(size = tag_label_size, face = "bold", margin = margin(0,0,1,0, unit="pt")))

dat3 <- predictions0 %>% filter(subset %in% "NAO-matched ensemble") %>% mutate(anom = ens_anom)
dat3 <- gauge_stns %>% left_join(dat3, by = "ID")
p3 <- myplotfun888(dat3)
p3 <- p3 +
  theme(
    axis.title.y = element_text(size = axis_title_size, margin = margin(0, 0, 0, 0, unit="cm")),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "cm")
  )
p3 <- p3 +
  labs(title = "c. NAO-matched ensemble") +
  theme(plot.title = element_text(size = tag_label_size, face = "bold", margin = margin(0,0,1,0, unit="pt")))

p <- p1 + p2 + p3 + plot_layout(ncol = 3, nrow = 1) & theme(legend.position = "bottom")
p = p + plot_layout(guides = "collect")

ggsave(file.path(output_dir, "figS4.png"), plot = p, width = 6, height = 4.5, units = "in", dpi=fig_dpi)

## ## ####################################################### ##
## ## ####################################################### ##
## ##
## ## Figure S4 presentation
## ##
## ## ####################################################### ##
## ## ####################################################### ##

## dat1 <- predictions0 %>% filter(subset %in% "Full ensemble") %>% mutate(anom = obs_anom)
## dat1 <- gauge_stns %>% left_join(dat1, by = "ID")
## p1 <- myplotfun888(dat1)

## dat2 <- predictions0 %>% filter(subset %in% "Full ensemble") %>% mutate(anom = ens_anom)
## dat2 <- gauge_stns %>% left_join(dat2, by = "ID")
## p2 <- myplotfun888(dat2)

## dat3 <- predictions0 %>% filter(subset %in% "NAO-matched ensemble") %>% mutate(anom = ens_anom)
## dat3 <- gauge_stns %>% left_join(dat3, by = "ID")
## p3 <- myplotfun888(dat3)

## p1 = p1 + theme(legend.text = element_text(size = 9),
##                 axis.title.y = element_text(size = 9),
##                 axis.text.y = element_text(size = 9),
##                 axis.text.x = element_text(size = 9))
## p2 = p2 + theme(legend.text = element_text(size = 9),
##                 axis.title.y = element_blank(),
##                 axis.text.y = element_blank(),
##                 axis.ticks.y = element_blank(),
##                 axis.text.x = element_text(size = 9))
## p3 = p3 + theme(legend.text = element_text(size = 9),
##                 axis.title.y = element_blank(),
##                 axis.text.y = element_blank(),
##                 axis.ticks.y = element_blank(),
##                 axis.text.x = element_text(size = 9))

## p <- p1 + p2 + p3 + plot_layout(ncol = 3, nrow = 1) & theme(legend.position = "bottom")
## p = p + plot_layout(guides = "collect")

## ## p$patches$plots[[1]] =
## ##   p$patches$plots[[1]] +
## ##   labs(tag = "a. Observed") +
## ##   theme(plot.tag.position = c(0.15, 1.03),
## ##         plot.tag = element_text(hjust = 0, size = tag_label_size, face="bold"))
## ## p$patches$plots[[2]] =
## ##   p$patches$plots[[2]] +
## ##   labs(tag = "b. Full ensemble") +
## ##   theme(plot.tag.position = c(0.025, 1.03),
## ##         plot.tag = element_text(hjust = 0, size = tag_label_size, face="bold"))
## ## p =
## ##   p +
## ##   labs(tag = "c. NAO-matched ensemble") +
## ##   theme(plot.tag.position = c(0.025, 1.03),
## ##         plot.tag = element_text(#vjust = -0.7,
## ##                                 hjust = 0, size = tag_label_size, face="bold"))

## ggsave(file.path(output_dir, "figS4_ppt.png"), plot = p, width = 6, height = 4.5, units = "in")
