#!/usr/bin/env Rscript

## Author : Simon Moulds
## Date   : Jan 2022

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(ggsignif)
library(scales)
library(arrow)
library(sf)
library(gamlss)
library(rnrfa)
library(yaml)

options(dplyr.summarise.inform = FALSE)

myplotfun0 <- function() {
  north_atlantic_box <-
    st_sf(a=1:2, geom=st_sfc(st_point(c(-10, 55)), st_point(c(25, 70))), crs=4326) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    densify(n=100) %>%
    st_transform(crs = 27700)

  ## TODO reproject
  p =
    ggplot() +
    geom_sf(
      data = europe_boundary,
      color=NA,
      fill="lightgrey"
    ) +
    geom_sf(
      data = uk_boundary,
      lwd = 0.1
    ) +
    geom_sf(
      data = north_atlantic_box,
      lwd = 0.25,
      ## color="black",
      fill=NA
    ) +
    coord_sf(
      xlim = c(-15, 40),
      ylim = c(50, 75),
      default_crs = st_crs(4326)
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
    )
  p
}

myplotfun1 <- function(x, legend_title = "MSSS") {
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
      ## size = 1.5,
      size = 2,
      lwd = 0.1,
      alpha = 0.8
    ) +
    facet_wrap(. ~ period, ncol = 1) +
    coord_sf(
      xlim = c(-8, 2),
      ylim = c(50, 59),
      default_crs = st_crs(4326)
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
        title=legend_title, #"MSSS",
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

## myplotfun2 <- function(x) {
##   p = ggplot(x, aes(x = model, y=msss)) +
##     geom_boxplot(
##       lwd = 0.25,
##       outlier.size = 0.25
##     ) +
##     scale_x_discrete(name = "") +
##     scale_y_continuous(
##       name="MSSS",
##       breaks=seq(-0.2, 0.8, by=0.2),
##       limits=c(-0.25, 0.85)
##     ) +
##     theme_bw() +
##     theme(
##       panel.grid = element_blank(),
##       axis.title = element_text(size = axis_title_size_small),
##       axis.text = element_text(size = axis_label_size_small)
##     )
##   p
## }

format_p <- function(x) {
  if (x <= 0.001) return("***")
  if (x <= 0.01) return("**")
  if (x <= 0.05) return("*")
  if (x > 0.05) return("NS.")
}

myplotfun2 <- function(x, legend_title = "MSSS") {
  x_median <-
    x %>%
    group_by(model) %>%
    summarize(md = median(skill), n = n()) %>%
    mutate(md_label = sprintf(md, fmt = "%#.3f"))

  x1 <- x %>% filter(model %in% "P") %>% `$`(skill)
  y1 <- x %>% filter(model %in% "PT") %>% `$`(skill)
  pval1 <- wilcox.test(y1, x1, alternative = "two.sided", paired = TRUE)$p.value

  x2 <- x %>% filter(model %in% "PT") %>% `$`(skill)
  y2 <- x %>% filter(model %in% "NAOPT") %>% `$`(skill)
  pval2 <- wilcox.test(y2, x2, alternative = "two.sided", paired = TRUE)$p.value

  x3 <- x %>% filter(model %in% "P") %>% `$`(skill)
  y3 <- x %>% filter(model %in% "NAOPT") %>% `$`(skill)
  pval3 <- wilcox.test(y3, x3, alternative = "two.sided", paired = TRUE)$p.value

  p = ggplot(x, aes(x = model, y=skill)) +
    ## stat_boxplot(coef = NULL) +
    geom_boxplot(
      lwd = 0.25,
      outlier.size = 0.25
    ) +
    geom_signif(
      comparisons = list(c("P", "PT"), c("PT", "NAOPT"), c("P", "NAOPT")),
      ## test = wilcox.test,
      ## test.args = list(paired = TRUE),
      ## comparisons = list(c("PT", "NAOPT_full")),
      annotation = c(format_p(pval1), format_p(pval2), format_p(pval3)), #, format_p(pval3)),
      map_signif_level = TRUE,
      step_increase = 0.06,
      tip_length = 0.02,
      size = 0.25,
      textsize = 2
    ) +
    geom_text(data = x_median, aes(x = model, y = md, label = md_label), nudge_y = 0.05, size = 2) +
    ## geom_text(data = skill_scores_median, aes(x = model, y = -0.2, label = n), nudge_y = 0.03, size = 1.5) +
    geom_text(data = x_median, aes(x = model, y = -0.2, label = n), nudge_y = -0.05, size = 2, fontface = "italic") +
    scale_x_discrete(name = "") +
    scale_y_continuous(
      name=legend_title,
      ## breaks=seq(-0.2, 1, by=0.2),
      ## limits=c(-0.25, 1.1)
      breaks=seq(-1, 1, by=0.2),
      limits=c(-1, 1)
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = axis_title_size_small),
      axis.text = element_text(size = axis_label_size_small)
    )
  p
}

myplotfun22 <- function(x, legend_title = "MSSS") {
  x_median <-
    x %>%
    group_by(subset, model) %>%
    summarize(md = median(skill), n = n()) %>%
    mutate(md_label = sprintf(md, fmt = "%#.3f"))
  p <- ggplot(x, aes(x = interaction(subset, model), y = skill)) +
    ## stat_boxplot(coef = NULL) +
    geom_boxplot(aes(colour = subset), lwd = 0.25, outlier.size = 0.25) +
    geom_line(aes(group = interaction(ID, model)), lwd = 0.15, alpha = 0.35, colour = "darkgrey") +
    scale_colour_discrete(
      name = "",
      labels = c("Full ensemble", "NAO-matched ensemble"),
      type = c("#FC8D62", "#66C2A5")
      ## type = rev(cbbPalette)
    ) +
    ## geom_signif(
    ##   ## comparisons = list(c("P", "PT"), c("PT", "NAOPT"), c("P", "NAOPT")),
    ##   comparisons = list(c("best_n", "full")),
    ##   map_signif_level = TRUE,
    ##   step_increase = 0.06,
    ##   tip_length = 0.01,
    ##   size = 0.25,
    ##   textsize = 1.5
    ## ) +
    geom_text(data = x_median, aes(x = interaction(subset, model), y = md, label = md_label), nudge_y = 0.05, size = 2) +
    ## geom_text(data = skill_scores_median, aes(x = model, y = -0.2, label = n), nudge_y = 0.03, size = 1.5) +
    ## geom_text(data = x_median, aes(x = model, y = -0.7, label = n), nudge_y = -0.05, size = 2, fontface = "italic") +
    ## geom_text(data = x_median, aes(x = interaction(subset, model), y = -0.7, label = n), nudge_y = -0.05, size = 2, fontface = "italic") +
    ## scale_x_discrete(name = "") +
    scale_x_discrete(name = "", labels = rep(levels(x$model), each = 2)) +
    scale_y_continuous(
      name=legend_title,
      ## breaks=seq(-0.2, 1, by=0.2),
      ## limits=c(-0.3, 1.1)
      breaks=seq(-1, 1, by=0.2),
      ## limits=c(-1, 1)
      limits=c(-0.8, 0.8)
      ## limits=c(-0.7, 0.7)
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = axis_title_size_small),
      axis.text = element_text(size = axis_label_size_small)
    )
  p
}

myplotfun3 <- function(x) {
  rdbu_pal = RColorBrewer::brewer.pal(9, "RdBu")
  cbbPalette <- RColorBrewer::brewer.pal(3, "Set2")
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
      aes(size = skill_diff, shape = model, fill = increase),
      alpha = 0.6
    ) +
    facet_wrap(. ~ period, ncol = 1) +
    coord_sf(
      xlim=c(-8, 2),
      ylim=c(50, 59),
      default_crs = st_crs(4326)
    ) +
    scale_shape_manual(
      name = "Model",
      ## values = c(21, 24, 22)
      values = c(21, 24)
    ) +
    scale_size_continuous(
      name = "Difference",
      ## breaks = c(0.2, 0.4, 0.6),
      ## labels = c("±0.2", "±0.4", "±0.6"),
      breaks = rep(c(0.2, 0.4, 0.6), 2),
      labels = c(c("±0.2", "±0.4", "±0.6"), rep("", 3)),
      limits = c(0, 0.75),
      range = c(0.1, 4)
    ) +
    ## scale_fill_manual(name = "Direction", labels = c("Decrease", "Increase"), values = rev(cbbPalette)) +
    scale_fill_discrete(
      name = "Direction",
      labels = c("Decrease", "Increase"),
      ## type = c("#FC8D62", "#66C2A5")
      ## type = c("#E1BE6A", "#40B0A6")
      type = cbbPalette[2:1]
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
    )
  p
}

myplotfun444 <- function(x, legend_title = "MSSS") {
  x_median <-
    x %>%
    group_by(model) %>%
    summarize(md = median(skill_diff), n = n()) %>%
    mutate(md_label = sprintf(md, fmt = "%#.3f"))
  p = ggplot(x, aes(x = model, y=skill_diff)) +
    ## stat_boxplot(coef = NULL) +
    geom_boxplot(
      lwd = 0.25,
      outlier.size = 0.25,
      position = position_dodge(.85)
    ) +
    ## geom_signif(
    ##   ## comparisons = list(c("P", "PT")), #, c("PT", "NAOPT"), c("P", "NAOPT")),
    ##   comparisons = list(c("P", "PT"), c("PT", "NAOPT"), c("P", "NAOPT")),
    ##   map_signif_level = TRUE,
    ##   step_increase = 0.06,
    ##   tip_length = 0.01,
    ##   size = 0.25,
    ##   textsize = 1.5
    ## ) +
    geom_text(data = x_median, aes(x = model, y = md, label = md_label), nudge_y = 0.05, size = 2) +
    ## geom_text(data = skill_scores_median, aes(x = model, y = -0.2, label = n), nudge_y = 0.03, size = 1.5) +
    ## geom_text(data = x_median, aes(x = model, y = -0.7, label = n), nudge_y = -0.05, size = 2, fontface = "italic") +
    scale_x_discrete(name = "") +
    scale_y_continuous(
      name=paste0("\u0394 ", legend_title),
      ## breaks=c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6),
      breaks=seq(-1, 1, by=0.2),
      ## limits=c(-0.7, 0.7)
      limits=c(-0.8, 0.8)
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = axis_title_size_small),
      axis.text = element_text(size = axis_label_size_small)
    )
  p
}

myplotfun5 <- function(x, legend_title = "MSSS") {
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
      ## size = 1.5,
      size = 1.5,
      lwd = 0.1,
      alpha = 0.8
    ) +
    ## geom_sf(
    ##   data = x,
    ##   aes(fill = skill),
    ##   shape = 21,
    ##   ## size = 1.5,
    ##   size = 2,
    ##   lwd = 0.1,
    ##   alpha = 0.8
    ## ) +
    facet_wrap(. ~ period, ncol = 1) +
    coord_sf(
      xlim = c(-8, 2),
      ylim = c(50, 59),
      default_crs = st_crs(4326)
    ) +
    scale_shape_manual(values = c(21, 24, 22)) +
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
    ) +
    guides(
      ## shape = "none",
      shape = guide_legend(
        title = "Model",
        title.position = "top",
        order = 1
      ),
      fill = guide_colorbar(
        title=legend_title,
        title.position="top",
        frame.colour = "black",
        ticks.colour = "black",
        frame.linewidth = 0.25,
        ticks.linewidth = 0.25,
        barwidth = 0.75,
        barheight = 10,
        order = 2
      )
    )
  p
}


myplotfun6 <- function(x) {
  ## We construct the legend using the model, so here we convert
  ## it to a factor and change the labels to those we want to
  ## display (defined in preamble)
  ## x = x %>% mutate(model = as.factor(model))
  ## x$model = do.call("recode_factor", c(list(x$model), model_display_names))
  cbbPalette <- RColorBrewer::brewer.pal(2, "Set2")
  obs =
    x %>%
    dplyr::select(year, obs) %>%
    ## dplyr::select(clim_season, year, obs) %>%
    distinct(year, .keep_all = TRUE) %>%
    mutate(type = "Observed")
  x = x %>% mutate(ID = paste0("ID=", ID))
  p = ggplot() +
    theme_bw() +
    geom_ribbon(
      aes(ymin=Q25, ymax=Q75, x=year, fill=subset), #model),
      alpha=0.5, data=x
    ) +
    geom_line(
      aes(y=Q50, x=year, colour=subset), data=x #model), data=x
    ) +
    scale_fill_manual(values = cbbPalette) +
    scale_color_manual(values = cbbPalette) +
    ## scale_fill_discrete(values = cbbPalette) +
    ## scale_colour_discrete(values = cbbPalette) +
    ## facet_wrap(. ~ ID, ncol = 1) + #, labeller = label_parsed) +
    ## labs(title = paste0("ID=", id)) +
    ylab(expression(Streamflow~(m^{3}~s^{-1}))) +
    xlab("") +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(expression(Streamflow~(m^{3}~s^{-1})), breaks = pretty_breaks()) +
    ## N.B. use alpha to create another legend
    ## https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
    ## geom_point(
    geom_line(
      aes(y=obs, x=year), #, alpha="Observed"),
      color = "black",
      data=obs,
      size = 1 #0.2
    ) +
    ## scale_alpha_manual(name=NULL, values=1, breaks="Observed") +
    ## ggtitle(sprintf("ID = %d", stn_id)) +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          legend.title = element_blank(),
          strip.background = element_blank())
  p = p + guides(#alpha = guide_legend(order = 1),
                 fill = guide_legend(order = 2, direction = "horizontal"),
                 color = guide_legend(order = 2, direction = "horizontal"))
  p
}

myplotfun777 <- function(x) {
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
      ylim = c(50, 60),
      default_crs = st_crs(4326)
    ) +
    theme_bw() #+
  p
}

myplotfun888 <- function(x) {
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
      aes(fill = anom), #, shape = model),
      shape = 21,
      size = 1.5,
      ## size = 2,
      lwd = 0.1,
      alpha = 0.8
    ) +
    coord_sf(
      xlim = c(-8, 2),
      ylim = c(50, 60),
      default_crs = st_crs(4326)
    ) +
    ## scale_shape_manual(values = c(21, 24, 22)) +
    scale_fill_stepsn(
      ## colours = rdbu_pal,
      ## breaks = seq(-0.8, 0.8, 0.2),
      ## values = scales::rescale(c(-0.8, 0, 0.8)),
      ## limits = c(-0.3, 0.9)
      ## breaks = seq(-0.2, 0.8, 0.2),
      ## values = scales::rescale(c(-0.2, 0, 0.8)),
      ## limits = c(-0.1, 0.9)
      colours = rev(RColorBrewer::brewer.pal(9, "RdBu")),
      breaks = seq(-2, 2, by=0.2),
      ## values = scales::rescale(-2, 2),
      limits = c(-2, 2),
      labels = c("", "-1.8", "", "", "-1.2", "", "", "-0.6", "", "", "0.0", "", "", "0.6", "", "", "1.2", "", "", "1.8", "")
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
      fill = guide_colorbar(
        title="DJFM Q95 anomaly",
        title.position="top",
        frame.colour = "black",
        ticks.colour = "black",
        frame.linewidth = 0.25,
        ticks = FALSE,
        ## ticks.linewidth = 0.25,
        barwidth = 16,
        ## barwidth = 12,
        barheight = 0.75,
        order = 2
      )
    )
  p
}

myplotfun9 <- function(x, legend_title = "MSSS") {
  ## p = ggplot(x, aes(x = msss_diff, y=obs_msss, color = model)) +
  ## p = ggplot(x, aes(x = obs_skill, y = skill_diff, color = model)) +
  p = ggplot(x, aes(x = R2, y = skill_diff, color = model)) +
    geom_point()  +
    scale_y_continuous(
      name=paste0("\u0394 ", legend_title, " (Predictions)"),
      breaks=c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8),
      limits=c(-0.8, 0.8)
    ) +
    scale_x_continuous(
      name=paste0(legend_title, " (Observations)"),
      ## breaks=c(-0.2, 0, 0.2, 0.4, 0.6, 0.8),
      ## limits=c(-0.2, 0.8)
      breaks=seq(0, 1, by = 0.2),
      limits=c(0, 1)
    ) +
    theme_bw()
  p
}

myplotfun1010 <- function(full_fcst, varname) {

  ## Plot 3 [analog of Fig 2c from Smith et al. 2020]
  plotdata =
    full_fcst %>%
    filter(variable %in% varname) %>%
    pivot_longer(c(-init_year, -variable, -period), names_to = "statistic", values_to = "value") %>%
    filter(statistic %in% c("obs", "ens_mean_lag", "ens_q95_lag", "ens_q05_lag")) %>%
    mutate(statistic = factor(
             statistic,
             levels = c("obs", "ens_mean_lag", "ens_q95_lag", "ens_q05_lag"),
             labels = c("Observed", "Modelled", "ens_q95", "ens_q05"))) %>%
    mutate(value = ifelse(init_year < 1964, NA, value))

  p3 = ggplot() +
    geom_ribbon(
      data = plotdata %>% filter(statistic %in% c("ens_q95", "ens_q05")) %>% pivot_wider(names_from = statistic, values_from = value),
      aes(x = init_year, ymin = ens_q05, ymax = ens_q95), fill = "red", alpha=0.25
    ) +
    geom_line(
      data = plotdata %>% filter(statistic %in% c("Observed", "Modelled")),
      aes(x = init_year, y = value, color = statistic)
    ) +
    geom_hline(yintercept=0, size=0.25) +
    facet_wrap(. ~ period, ncol = 1) + #, labeller = label_parsed) +
    ## scale_y_continuous(
    ##   name="AMV anomaly (K)",
    ##   breaks=c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3),
    ##   limits=c(-0.3, 0.3)
    ## ) +
    scale_x_continuous(
      name = "",
      breaks = seq(1960, 2000, 10),
      limits = c(1960, 2005)
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(size = strip_label_size),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = legend_label_size),
          axis.title.y = element_text(size = axis_title_size),
          axis.text.y = element_text(size = axis_label_size_small),
          axis.text.x = element_text(size = axis_label_size_small))


  p3
}

myplotfun11 <- function(x) {
  ## We construct the legend using the model, so here we convert
  ## it to a factor and change the labels to those we want to
  ## display (defined in preamble)
  ## cbbPalette <- c(
  ##   "#000000", "#E69F00", "#56B4E9", "#009E73",
  ##   "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  ## )
  ## cbbPalette <- cbbPalette[c(2, 4, 8)]
  cbbPalette <- RColorBrewer::brewer.pal(3, "Set2")
  obs =
    x %>%
    dplyr::select(year, obs) %>%
    distinct(year, .keep_all = TRUE) %>%
    mutate(type = "Observed")
  x = x %>% mutate(ID = paste0("ID=", ID))
  p = ggplot() +
    theme_bw() +
    geom_ribbon(
      aes(ymin=Q25, ymax=Q75, x=year, fill=model),
      alpha=0.5, data=x
    ) +
    geom_line(
      aes(y=Q50, x=year, colour=model), data=x
    ) +
    ## facet_wrap(. ~ ID, ncol = 1) +
    ylab(expression(Streamflow~(m^{3}~s^{-1}))) +
    xlab("") +
    scale_x_continuous(breaks = pretty_breaks(), limits = c(1960, 2005)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_fill_manual(values = cbbPalette) +
    scale_color_manual(values = cbbPalette) +
    ## N.B. use alpha to create another legend
    ## https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
    geom_line(
      aes(y=obs, x=year), #, alpha="Observed"),
      color = "black",
      data=obs,
      size = 1 #0.2
    ) +
    ## geom_point(
    ##   aes(y=obs, x=year), #, alpha="Observed"),
    ##   color = "black",
    ##   data=obs,
    ##   size = 0.2
    ## ) +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          legend.title = element_blank(),
          strip.background = element_blank())
  p = p +
    guides(
      fill = guide_legend(
        title = "Model",
        title.position = "top",
        title.hjust = 0,
        order = 2,
        direction = "horizontal"
      ),
      color = guide_legend(
        title = "Model",
        title.position = "top",
        title.hjust = 0,
        order = 2,
        direction = "horizontal"
      )
    )
  p
}

myplotfun12 <- function(x) {
  p = ggplot(x, aes(x = model, y=msss, fill=subset)) +
    ## stat_boxplot(coef = NULL) +
    geom_boxplot(
      lwd = 0.25,
      outlier.size = 0.25,
      position = position_dodge(.85)
    ) + scale_x_discrete(name = "") +
    facet_wrap(. ~ period, ncol = 1) + #, labeller = label_parsed) +
    scale_y_continuous(
      name="MSSS",
      breaks=seq(-0.2, 0.8, by=0.2),
      limits=c(-0.25, 0.85)
    ) +
    theme_bw() +
    theme(strip.background = element_blank())
  p
}

## myplotfun7 = function(x) {
##   p = ggplot(x, aes(x = subset, y=msss, fill=lag)) +
##     geom_boxplot(
##       lwd = 0.25,
##       outlier.size = 0.25,
##       position = position_dodge(.85)
##     ) +
##     scale_x_discrete(
##       name="",
##       labels=c("CMIP5", "CMIP6", "CMIP5 & 6")
##       ## labels=c("CMIP5", "CMIP6", "CMIP5 & CMIP6")
##     ) +
##     facet_wrap(. ~ period, ncol = 1) + #, labeller = label_parsed) +
##     scale_y_continuous(
##       name="MSSS",
##       breaks=seq(-0.2, 0.8, by=0.2),
##       limits=c(-0.25, 0.85)
##     ) +
##     theme_bw() +
##     theme(strip.background = element_blank())
##   p
## }

## myplotfun8 = function(x) {
##   p = ggplot(x, aes(x = model, y=msss, fill=lag)) +
##     geom_boxplot(
##       lwd = 0.25,
##       outlier.size = 0.25,
##       position = position_dodge(.85)
##     ) +
##     ## scale_x_discrete(
##     ##   name="",
##     ##   labels=c("CMIP5", "CMIP6", "CMIP5 & 6")
##     ##   ## labels=c("CMIP5", "CMIP6", "CMIP5 & CMIP6")
##     ## ) +
##     facet_wrap(. ~ period, ncol = 1) + #, labeller = label_parsed) +
##     scale_y_continuous(
##       name="MSSS",
##       breaks=seq(-0.2, 0.8, by=0.2),
##       limits=c(-0.25, 0.85)
##     ) +
##     theme_bw() +
##     theme(strip.background = element_blank())
##   p
## }

## myplotfun66 = function(x) {
##   p = ggplot(x, aes(x = model, y=msss)) +
##     geom_boxplot(
##       lwd = 0.25,
##       outlier.size = 0.25
##     ) +
##     scale_x_discrete(name = "") +
##     scale_y_continuous(
##       name="MSSS",
##       breaks=seq(-0.2, 0.8, by=0.2),
##       limits=c(-0.25, 0.85)
##     ) +
##     theme_bw() +
##     theme(
##       panel.grid = element_blank(),
##       axis.title = element_text(size = axis_title_size_small),
##       axis.text = element_text(size = axis_label_size_small)
##     )
##   p
## }
