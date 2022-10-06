#!/usr/bin/env Rscript

## Author : Simon Moulds
## Date   : Jan 2022

library(tidyverse)
library(ggsignif)
library(lubridate)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(scales)
library(arrow)
library(sf)
library(yaml)
library(rnrfa)

options(dplyr.summarise.inform = FALSE)
options(bitmapType = 'cairo') # For server

## source("workflow/scripts/external/R/utils.R")

## Extract configuration info
if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=TRUE)
  inputdir <- args[1]
  aggregation_period = args[2]
  outputroot <- args[3]
  config <- read_yaml(args[4])
  args = commandArgs()
  m <- regexpr("(?<=^--file=).+", args, perl=TRUE)
  cwd <- dirname(regmatches(args, m))
}
## TODO put these functions in an R package
source(file.path(cwd, "utils.R"))

## #################################### ##
## Set some variables needed for plotting
## #################################### ##

## Spatial data for plotting
uk_boundary = st_read(
  file.path(
    inputdir,
    "CNTR_RG_01M_2020_4326.shp"
  )) %>%
  filter(CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

europe_boundary = st_read(
  file.path(
    inputdir,
    "CNTR_RG_01M_2020_4326.shp"
  )) %>%
  filter(!CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

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

plot_outputdir = file.path(outputroot, "fig", aggregation_period)
if (dir.exists(plot_outputdir)) {
  unlink(plot_outputdir, recursive = TRUE)
}
dir.create(plot_outputdir, recursive = TRUE)

## #################################### ##
## Load data
## #################################### ##

outputroot <- file.path(outputroot, "analysis", "hindcast")
gamlss_datadir <- file.path(outputroot, "gamlss", aggregation_period, "prediction")
## lstm_datadir <- file.path(outputroot, "lstm", aggregation_period, "prediction")
xgboost_datadir <- file.path(outputroot, "xgboost", aggregation_period, "prediction")
tabnet_datadir <- file.path(outputroot, "tabnet", aggregation_period, "prediction")
## metadata <- read_parquet(file.path(outputroot, "nrfa-metadata.parquet"))

## Join data frames and reshape
gamlss_output <- open_dataset(gamlss_datadir) %>%
  collect() %>%
  rename(Q95_exp = Q_95_exp, Q95_obs = Q_95_obs) %>%
  dplyr::select(Q95_obs, Q95_exp, year, ID, date, model)

## lstm_output <- open_dataset(lstm_datadir) %>%
##   collect() %>%
##   rename(Q95_exp = Q95_sim) %>%
##   mutate(model = "LSTM")

xgboost_output <- open_dataset(xgboost_datadir) %>%
  collect() %>%
  mutate(model = "XGBoost")

tabnet_output <- open_dataset(tabnet_datadir) %>%
  collect() %>%
  mutate(model = "TabNet")

## #################################### ##
## Get skill scores for each model
## #################################### ##

compute_skill_scores <- function(x, ...) {
  models <- unique(x$model) %>% sort()
  station_ids <- unique(x$ID) %>% sort()
  skill_scores <- list()
  for (i in 1:length(models)) {
    xx <- x %>% filter(model %in% models[i])
    model_skill_scores <- list()
    print(paste0("Computing skill scores for model ", models[i]))
    pb <- txtProgressBar(min = 0, max = length(station_ids), initial = 0)
    for (j in 1:length(station_ids)) {
      id <- station_ids[j]
      yy <- xx %>% filter(ID %in% id) %>% arrange(year)
      skill <- mean_square_error_skill_score(yy$Q95_obs, yy$Q95_exp) %>% as_tibble()
      skill <- skill %>% mutate(ID = id, .before = names(.)[1])
      model_skill_scores[[j]] <- skill
      setTxtProgressBar(pb, j)
    }
    close(pb)
    model_skill_scores <- do.call("rbind", model_skill_scores) %>% mutate(model = models[i], .before = ID)
    skill_scores[[i]] <- model_skill_scores
  }
  skill_scores <- do.call("rbind", skill_scores)
  skill_scores
}

gamlss_skill_scores <-
  compute_skill_scores(gamlss_output) %>%
  arrange(ID)

gamlss_skill_scores_best <-
  gamlss_skill_scores %>%
  group_by(ID) %>%
  filter(acc == max(acc))

gamlss_output <-
  gamlss_skill_scores_best %>%
  dplyr::select(model, ID) %>%
  left_join(gamlss_output)

output <- rbind(
  xgboost_output,
  tabnet_output,
  ## lstm_output,
  gamlss_output
) %>%
  mutate(date = as.Date(date)) %>%
  arrange(ID, date, model) %>%
  mutate(submodel = model) %>% # So we don't lose the information completely
  mutate(model = ifelse(str_detect(model, "GAMLSS_(.*)_(full|best_n)"), "GAMLSS", model))

## lstm_skill_scores <- compute_skill_scores(lstm_output)
xgboost_skill_scores <- compute_skill_scores(xgboost_output)
tabnet_skill_scores <- compute_skill_scores(tabnet_output)

## Combine skill scores in a single data frame
skill_scores <- rbind(
  xgboost_skill_scores,
  tabnet_skill_scores,
  ## lstm_skill_scores,
  gamlss_skill_scores_best
)

## Provide a catch-all model name for GAMLSS
skill_scores <-
  skill_scores %>%
  mutate(submodel = model) %>% # So we don't lose the information completely
  mutate(model = ifelse(str_detect(model, "GAMLSS_(.*)_(full|best_n)"), "GAMLSS", model))

## Only include stations where all models are available
models <- c("GAMLSS", "TabNet", "XGBoost")
n_models <- length(models)
common_ids <-
  skill_scores %>%
  group_by(ID) %>%
  summarize(n = n()) %>%
  filter(n == 3) %>% `$`(ID)

skill_scores <- skill_scores %>% filter(ID %in% common_ids)

## Add the best model
skill_scores_best <- skill_scores %>% group_by(ID) %>% filter(acc == max(acc)) %>% mutate(model = "Best")
skill_scores <- rbind(skill_scores, skill_scores_best)

station_ids <- unique(skill_scores$ID) %>% sort()
skill_scores <- skill_scores %>% mutate(model = factor(model, levels = c("GAMLSS", "XGBoost", "TabNet", "Best")))

## ################################### ##
## Plot 1: Boxplot comparison of models
## ################################### ##

## Need gauge station locations
gauge_stns =
  catalogue() %>%
  rename(ID = id, area = "catchment-area") %>%
  filter(ID %in% station_ids) %>%
  dplyr::select(ID, name, area, latitude, longitude) %>%
  st_as_sf(coords=c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_transform(27700)

plotfun1 <- function(x, legend_title = "ACC") {
  rdbu_pal = RColorBrewer::brewer.pal(9, "RdBu")
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
      aes(fill = skill, shape = model),
      size = 1.5,
      lwd = 0.1,
      alpha = 0.8
    ) +
    ## facet_wrap(. ~ period, ncol = 1) +
    coord_sf(
      xlim = c(-8, 2),
      ylim = c(50, 59)#,
      ## default_crs = st_crs(4326)
    ) +
    scale_shape_manual(values = c(21, 24, 22)) +
    scale_fill_stepsn(
      colours = rev(rdbu_pal),
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
      legend.position = "bottom",
      legend.box = "vertical",
      legend.justification = "left",
      legend.box.just = "left",
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_label_size),
      strip.text = element_blank(),
      panel.grid.major = element_line(size = 0.25),
      axis.text = element_text(size = axis_label_size_small),
    ) +
    guides(
      shape = guide_legend(
        title = "Model",
        title.position = "top",
        order = 1
      ),
      fill = guide_colorbar(
        title=legend_title,
        ## title="MSSS",
        title.position="top",
        frame.colour = "black",
        ticks.colour = "black",
        frame.linewidth = 0.25,
        ticks.linewidth = 0.25,
        barwidth = 12,
        barheight = 0.75,
        order = 2
      )
    )
  p
}

skill_scores <- skill_scores %>% mutate(skill = acc)
skill_scores_subset <- skill_scores %>% filter(!model %in% "Best") %>% group_by(ID) %>% filter(skill == max(skill))
skill_scores_subset <- gauge_stns %>% left_join(skill_scores_subset, by = "ID")
p1 <- plotfun1(skill_scores_subset)

## Add number of stations which perform best by model
n_best_model <- table(skill_scores_subset$model)
n_gamlss <- n_best_model[["GAMLSS"]]
n_xgboost <- n_best_model[["XGBoost"]]
n_tabnet <- n_best_model[["TabNet"]]
labs <- c(paste0("italic(n)==", n_gamlss), paste0("italic(n)==", n_xgboost), paste0("italic(n)==", n_tabnet))
d <- data.frame(x = c(0, 0, 0), y = c(59, 58.5, 58), lab = labs, model = c("GAMLSS", "XGBoost", "TabNet"))
p1 <- p1 +
  geom_point(data = d, aes(x, y, shape = model), size = 1, lwd = 0.1, show.legend = FALSE) +
  geom_text(data = d, aes(x, y, label = lab), parse = TRUE, hjust = 0, nudge_x = 0.3, size = 2) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(axis.title = element_blank())

rdbu_pal = brewer.pal(9, "RdBu")
p1 <- p1 +
  scale_fill_stepsn(
    colours = rev(rdbu_pal)[3:9],
    breaks = seq(-0.4, 0.8, 0.2),
    limits = c(-0.3, 0.9)
  )
## p1

## skill_scores_i = skill_scores %>% filter(model %in% "GAMLSS") %>% mutate(skill = acc)
## p1 <- plotfun1(skill_scores_i)
## skill_scores_i = skill_scores %>% filter(model %in% "TabNet") %>% mutate(skill = acc)
## p2 <- plotfun1(skill_scores_i)
## skill_scores_i = skill_scores %>% filter(model %in% "XGBoost") %>% mutate(skill = acc)
## p3 <- plotfun1(skill_scores_i)

## plotfun1 <- function(x, ...) {
##   ## Boxplot comparison
##   p = ggplot(x, aes(x = model, y=acc)) +
##     ## stat_boxplot(coef = NULL) +
##     geom_boxplot(
##       lwd = 0.25,
##       outlier.size = 0.25
##     ) +
##     scale_x_discrete(name = "") +
##     ## scale_y_continuous(
##     ##   name=legend_title,
##     ##   ## breaks=seq(-0.2, 1, by=0.2),
##     ##   ## limits=c(-0.25, 1.1)
##     ##   breaks=seq(-1, 1, by=0.2),
##     ##   limits=c(-1, 1)
##     ## ) +
##     theme_bw() +
##     theme(
##       panel.grid = element_blank(),
##       axis.title = element_text(size = axis_title_size_small),
##       axis.text = element_text(size = axis_label_size_small)
##     )
##   p
## }
##

plotfun2 <- function(x, legend_title = "ACC") {
  p <- ggplot(x, aes(x = model, y = skill)) + #interaction(subset, model), y = skill)) +
    ## geom_boxplot(aes(colour = subset), lwd = 0.25, outlier.size = 0.25) +
    geom_boxplot(lwd = 0.25, outlier.size = 0.25) +
    ## geom_line(aes(group = ID), lwd = 0.15, alpha = 0.35, colour = "darkgrey") +
    ## scale_colour_discrete(
    ##   name = "",
    ##   labels = c("Full ensemble", "NAO-matched ensemble"),
    ##   type = c("#FC8D62", "#66C2A5")
    ##   ## type = rev(cbbPalette)
    ## ) +
    ## geom_signif(
    ##   ## comparisons = list(c("P", "PT"), c("PT", "NAOPT"), c("P", "NAOPT")),
    ##   comparisons = list(c("best_n", "full")),
    ##   map_signif_level = TRUE,
    ##   step_increase = 0.06,
    ##   tip_length = 0.01,
    ##   size = 0.25,
    ##   textsize = 1.5
    ## ) +
    ## scale_x_discrete(name = "") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = axis_title_size_small),
      axis.text = element_text(size = axis_label_size_small)
    )
  p
}

## Signif [one-sided Wilcoxon signed rank test]
format_p <- function(x) {
  if (x <= 0.001) return("***")
  if (x <= 0.01) return("**")
  if (x <= 0.05) return("*")
  if (x > 0.05) return("NS.")
}

print("e")
skill_measure = "acc"
x <- skill_scores %>%
  filter(model %in% "GAMLSS") %>%
  `[[`(skill_measure)
y <- skill_scores %>%
  filter(model %in% "XGBoost") %>%
  `[[`(skill_measure)
pval1 = wilcox.test(x, y, paired = TRUE, alternative = "two.sided")$p.value %>% format_p

x <- skill_scores %>%
  filter(model %in% "GAMLSS") %>%
  `[[`(skill_measure)
y <- skill_scores %>%
  filter(model %in% "TabNet") %>%
  `[[`(skill_measure)
pval2 = wilcox.test(y, x, paired = TRUE, alternative = "two.sided")$p.value %>% format_p

x <- skill_scores %>%
  filter(model %in% "XGBoost") %>%
  `[[`(skill_measure)
y <- skill_scores %>%
  filter(model %in% "TabNet") %>%
  `[[`(skill_measure)
pval3 = wilcox.test(y, x, paired = TRUE, alternative = "two.sided")$p.value %>% format_p

print("f")
skill_scores_median <-
  skill_scores %>%
  group_by(model) %>%
  summarize(md = median(skill), n = n()) %>%
  mutate(md_label = sprintf(md, fmt = "%#.3f"))

print("g")
p2 <- plotfun2(skill_scores)
print("h")
p2 <- p2 +
  geom_signif(
    y_position = c(0.84, 0.92, 0.76), xmin = c(1, 1, 2), xmax = c(2, 3, 3),
    annotation = c(pval1, pval2, pval3), #, format_p(pval3)),
    step_increase = 0.06,
    tip_length = 0.02,
    size = 0.25,
    textsize = 2) +
  geom_text(data = skill_scores_median, aes(x = model, y = md, label = md_label), nudge_y = 0.05, size = 2) +
  geom_text(data = skill_scores_median, aes(x = model, y = -0.6, label = n), nudge_y = -0.05, size = 2, fontface = "italic") +
  scale_x_discrete(name = "", labels = levels(skill_scores$model)) +
  scale_y_continuous(
    name="ACC",
    breaks=seq(-0.8, 1, by=0.2),
    limits=c(-0.8, 1)
  ) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = legend_label_size),
        legend.direction = "vertical",
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.box = "vertical",
        legend.box.margin = margin(-0.5, 0, 0, 0, unit = "cm"),
        aspect.ratio = 1
        )

print("i")
p <- p1 + p2 + plot_layout(nrow = 1, ncol = 2)
p <- p + plot_layout(widths = c(2, 2))
p$patches$plots[[1]] =
  p$patches$plots[[1]] +
  labs(tag = "a") +
  theme(plot.tag.position = c(0.1, 1.01),
        plot.tag = element_text(size = tag_label_size, face="bold"))
p =
  p +
  labs(tag = "b") +
  theme(plot.tag.position = c(0.125, 0.87),
        plot.tag = element_text(vjust = -0.7, size = tag_label_size, face="bold"))

ggsave(file.path(plot_outputdir, "fig1.png"), plot = p, width = 6, height = 6, units = "in")

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
      size = 0.6 #0.2
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
skill_scores_subset <- skill_scores %>% filter(!model %in% "Best") %>% arrange(desc(acc))
ids_best = head(skill_scores_subset$ID, n=5)

p1 <- plotfun2(output, ids_best[1])
p2 <- plotfun2(output, ids_best[2])
p3 <- plotfun2(output, ids_best[3])
p4 <- plotfun2(output, ids_best[4])
p5 <- plotfun2(output, ids_best[5])
## p6 <- plotfun2(output, 84019)

format_p_value = function(p_value) {
  if (p_value < 0.01) {
    return("(P < 0.01)")
  } else {
    return(paste0("(P = ", sprintf(p_value, fmt = "%#.2f"), ")"))
  }
}

## Regionalization:
## ================
##
## See: https://stackoverflow.com/q/67579831 for extracting contours from raster
##      [raster::rasterToContour / grDevices::contourLines]
## See: https://swilke-geoscience.net/post/2020-09-10-kriging_with_r/kriging/ for kriging using gstat
##
## Basic idea:
## ===========
## * Use ordinary kriking to interpolate skill
## * Extract zero skill contour
## * Remove section of raster with zero skill
## * Interpolate over the resulting grid

## make_annotation = function(acc, id) {
##   acc_full <- acc %>% filter(ID %in% id & subset %in% "Full ensemble")
##   acc_matched <- acc %>% filter(ID %in% id & subset %in% "NAO-matched ensemble")
##   annotation = paste0(
##     paste0(
##       "ACC (Full) = ", sprintf(acc_full$acc, fmt = "%#.2f"), " ", format_p_value(acc_full$acc_p), "\n"
##     ),
##     paste0(
##       "ACC (NAO-matched) = ", sprintf(acc_matched$acc, fmt = "%#.2f"), " ", format_p_value(acc_matched$acc_p), "\n"
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

## annotation_size = 4
## annotation_rel_pos = 1
## yrange <- get_y_range(p1)
## yrange[2] <- yrange[2] * 1.025
## p1 <- p1 +
##   ylim(yrange) +
##   annotate(
##     geom = "text",
##     x = 1960,
##     y = yrange[1] + diff(yrange) * annotation_rel_pos,
##     label = make_annotation(acc, ids_best[1]),
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
##     label = make_annotation(acc, ids_best[2]),
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
##     label = make_annotation(acc, ids_best[3]),
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
##     label = make_annotation(acc, ids_best[4]),
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
##     label = make_annotation(acc, ids_best[5]),
##     hjust=0,
##     vjust=1,
##     size = annotation_size / ggplot2::.pt
##   )
myplotfun3 <- function(x) {
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
      ## shape = 21,
      size = 0,
      ## lwd = 0.1
      ## alpha = 0.8
    ) +
    coord_sf(
      xlim = c(-8, 2),
      ylim = c(50, 60)#,
      ## default_crs = st_crs(4326)
    ) +
    theme_bw() #+
  p
}

gauge_stns_subset = gauge_stns %>% filter(ID %in% ids_best)
p6 = myplotfun3(gauge_stns_subset)
p6 =
  p6 +
  ggrepel::geom_label_repel(
             data=gauge_stns_subset,
             aes(label=ID, geometry=geometry),
             stat="sf_coordinates",
             min.segment.length = 0,
             size = 1.75)

p1 = p1 + theme(panel.grid = element_blank(),
                strip.text = element_text(size = strip_label_size),
                ## legend.title = element_text(size = legend_title_size),
                legend.text = element_text(size = legend_label_size),
                axis.title.y = element_text(size = axis_title_size, angle = 90),
                axis.text.y = element_text(size = axis_label_size_small),
                axis.text.x = element_text(size = axis_label_size_small))
p2 = p2 + theme(panel.grid = element_blank(),
                strip.text = element_text(size = strip_label_size),
                ## legend.title = element_text(size = legend_title_size),
                legend.text = element_text(size = legend_label_size),
                ## axis.title.y = element_text(size = axis_title_size, angle = 90),
                axis.title.y = element_blank(),
                axis.text.y = element_text(size = axis_label_size_small),
                axis.text.x = element_text(size = axis_label_size_small))
p3 = p3 + theme(panel.grid = element_blank(),
                strip.text = element_text(size = strip_label_size),
                ## legend.title = element_text(size = legend_title_size),
                legend.text = element_text(size = legend_label_size),
                ## axis.title.y = element_text(size = axis_title_size, angle = 90),
                axis.title.y = element_blank(),
                axis.text.y = element_text(size = axis_label_size_small),
                axis.text.x = element_text(size = axis_label_size_small))
p4 = p4 + theme(panel.grid = element_blank(),
                strip.text = element_text(size = strip_label_size),
                ## legend.title = element_text(size = legend_title_size),
                legend.text = element_text(size = legend_label_size),
                axis.title.y = element_text(size = axis_title_size, angle = 90),
                ## axis.title.y = element_text(size = axis_title_size),
                axis.text.y = element_text(size = axis_label_size_small),
                axis.text.x = element_text(size = axis_label_size_small))
p5 = p5 + theme(panel.grid = element_blank(),
                strip.text = element_text(size = strip_label_size),
                ## legend.title = element_text(size = legend_title_size),
                legend.text = element_text(size = legend_label_size),
                axis.title.y = element_blank(),
                axis.text.y = element_text(size = axis_label_size_small),
                axis.text.x = element_text(size = axis_label_size_small))
p6 = p6 + theme(panel.grid.major = element_line(size = 0.25),
                axis.title = element_blank(),
                axis.text = element_text(size = axis_label_size_small))

p = p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, nrow = 2, widths = c(2, 2, 2)) & theme(legend.position = "bottom")
p = p + plot_layout(guides = "collect")

ggsave(file.path(plot_outputdir, "fig2.png"), plot = p, width = 5, height = 5, units = "in")
