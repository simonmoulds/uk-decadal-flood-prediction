#!/usr/bin/env Rscript

library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(scales)
library(arrow)
library(sf)
library(smoothr)
library(yaml)
library(viridis)

options(dplyr.summarise.inform = FALSE)
options(bitmapType = 'cairo')

if (exists("snakemake")) {
  config <- snakemake@config
  aggregation_period <- snakemake@wildcards[["aggr"]]
  outputroot <- snakemake@params[["outputroot"]]
  snakemake@source("utils.R")
  snakemake@source("plotting.R")
} else {
  ## FOR TESTING:
  config = read_yaml('config/config_2.yml')
  aggregation_period = "yr2"
  outputroot = 'results'
  cwd = 'workflow/decadal-prediction-scripts/R'
  source(file.path(cwd, "utils.R"))
  source(file.path(cwd, "plotting.R"))
}

config[["modelling"]] <- parse_config_modelling(config)

fig_dpi <- 600
skill_measure <- "crpss"
all_skill_measures <- c(
  "crps_fcst", "crps_ens_fcst",
  "crps_climat", "crpss", "aic", "r"
)
input_dir <- file.path(outputroot, "input")
output_dir <- file.path(outputroot, "fig", aggregation_period)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

aggregation_period_label_list <- list(
  "yr2" = "Year 2",
  "yr2to5_lag" = "Year 2-5",
  "yr6to9_lag" = "Year 6-9",
  "yr2to9_lag" = "Year 2-9"
)

obs_aggregation_period_list <- list(
  "yr2" = "yr2",
  "yr2to5_lag" = "yr2to5_lag",
  "yr6to9_lag" = "yr6to9_lag",
  "yr2to9_lag" = "yr2to9_lag"
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
obs_skill_scores <- load_skill_scores(
  config, "observed_Q95", obs_aggregation_period, season="DJFM"
)
skill_scores <- load_skill_scores(
  config, "hindcast_Q95", aggregation_period, season="DJFM"
)
station_ids <- skill_scores$ID %>% unique()

## Load model fit metrics for hindcast experiment
fit <- load_model_fit(
  config, "hindcast_Q95", aggregation_period, season="DJFM"
) %>% mutate(kurtosis = kurtosis - 3)

## For spatial plots:
uk_boundary =
  st_read("resources/CNTR_RG_01M_2020_4326.shp") %>%
  filter(CNTR_NAME %in% "United Kingdom") %>%
  st_transform(crs = 27700)

europe_boundary =
  st_read("resources/CNTR_RG_01M_2020_4326.shp") %>%
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

skill_scores_subset =
  skill_scores %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model)) %>%
  filter(subset %in% "full", period %in% aggregation_period) %>%
  mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
  mutate(skill = !!sym(skill_measure))

obs_skill_scores_subset <-
  obs_skill_scores %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model)) %>%
  filter(period %in% obs_aggregation_period) %>%
  mutate(period = factor(period, levels = obs_aggregation_period, labels = obs_aggregation_period_label)) %>%
  mutate(skill = !!sym(skill_measure))

## Schematic
example_discharge_data <-
  open_dataset(
    file.path(input_dir, "discharge/NRFA"),
    partitioning=c("clim_season", "ID")
  ) %>%
  collect() %>%
  filter(ID %in% 21017 & clim_season %in% "DJFM" & season_year %in% 1978:1991)
pp1 <- plot_grl2023_fig1a_b(example_discharge_data)

## Raw ensemble
skill1 <- myfun(skill_scores_subset)
p1 <- plot_grl2023_fig1c_d(na.omit(skill1), legend_title = toupper(skill_measure))

## Observed data (i.e. perfect predictors)
skill1 <- myfun(obs_skill_scores_subset)
p2 <- plot_grl2023_fig1c_d(na.omit(skill1), legend_title = toupper(skill_measure))
p2 <- p2 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

## North Atlantic region
p3 <- plot_grl2023_fig1e()
p3 <- p3 + theme(plot.margin = margin(0, 0, 0, 0, unit="cm"))

## Add labels to plots
p1 <- p1 +
  labs(title = "c") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

p2 <- p2 +
  labs(title = "d") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

p3 <- p3 +
  labs(title = "e") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

design <- "
AAAABBBBCCC
AAAABBBBCCC
AAAABBBBDDD
AAAABBBBDDD
"
pp2 <-
  p1 + p2 + p3 +
  guide_area() +
  plot_layout(design=design, guides="collect") &
  theme(
    legend.spacing.y = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0, unit="cm"),
    plot.background = element_blank()
  )

## Now use cowplot to join the two parts
p = plot_grid(pp1, pp2, nrow=2, align="v", rel_heights = c(1, 1.25))
ggsave(
  file.path(output_dir, "fig1.png"),
  plot = p, width = 6, height = 6.05, units = "in", dpi=fig_dpi
)

## ####################################################### ##
## ####################################################### ##
##
## Figure 2
##
## ####################################################### ##
## ####################################################### ##

skill_scores_subset =
  skill_scores %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model)) %>%
  filter(subset %in% c("best_n", "full")) %>%
  filter(period %in% aggregation_period) %>%
  mutate(period = factor(period, levels = aggregation_period, labels = aggregation_period_label)) %>%
  mutate(skill = !!sym(skill_measure)) %>%
  dplyr::select(-all_of(all_skill_measures)) %>%
  pivot_wider(names_from = subset, values_from = skill) %>%
  mutate(skill_diff = best_n - full)

p4 <- plot_grl2023_fig2d(
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

p2 <- plot_grl2023_fig2b(skill)
p2 <- p2 +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.margin=margin(0, 0, 0, 0),
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(0.25, 'cm'),
    legend.text = element_text(margin = margin(0, 0, 0, 0))
  )

skill_scores_subset <-
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

p3 <- plot_grl2023_fig2c(
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
    annotation = c(format_p(pval1), format_p(pval2)),
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

skill <- myfun(
  skill_scores_subset %>%
  filter(!model %in% "NAOPT" & subset %in% "NAO-matched ensemble")
) %>% na.omit()
p1 <- plot_grl2023_fig2a(skill, legend_title = toupper(skill_measure))
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

design = "
AAAABBBB
AAAABBBB
AAAABBBB
AAAABBBB
CCCCDDDD
CCCCDDDD"

p1 <- p1 +
  labs(title = "a") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

p2 <- p2 +
  labs(title = "b") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

p3 <- p3 +
  labs(title = "c") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

p4 <- p4 +
  labs(title = "d") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

p <- p1 + p2 + p3 + p4 + plot_layout(design = design)

ggsave(
  file.path(output_dir, "fig2.png"),
  plot = p, width = 6, height = 5.6, units = "in", dpi=fig_dpi
)

## ####################################################### ##
## ####################################################### ##
##
## Figure 3
##
## ####################################################### ##
## ####################################################### ##

skill_scores_subset <-
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

skill <-
  skill_scores_subset %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst)) %>%
  arrange(desc(best_n))

ids_best <- head(skill, n=5)$ID
models_best <- as.character(head(skill, n=5)$model)

predictions <-
  load_model_predictions(
    config, "hindcast_Q95", aggregation_period, season="DJFM"
  ) %>%
  mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))
predictions <- predictions %>% mutate(obs = Q_95_obs, exp = Q_95_exp)

p1 <- predictions %>%
  filter(ID %in% ids_best[1] & model %in% models_best[1]) %>%
  plot_grl2023_fig3()

p1 <- p1 +
  labs(title = paste0("ID=", ids_best[1])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p2 <- predictions %>%
  filter(ID %in% ids_best[2] & model %in% models_best[2]) %>%
  plot_grl2023_fig3()

p2 <- p2 +
  labs(title = paste0("ID=", ids_best[2])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p3 <- predictions %>%
  filter(ID %in% ids_best[3] & model %in% models_best[3]) %>%
  plot_grl2023_fig3()

p3 <- p3 +
  labs(title = paste0("ID=", ids_best[3])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p4 <- predictions %>%
  filter(ID %in% ids_best[4] & model %in% models_best[4]) %>%
  plot_grl2023_fig3()

p4 <- p4 +
  labs(title = paste0("ID=", ids_best[4])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p5 <- predictions %>%
  filter(ID %in% ids_best[5] & model %in% models_best[5]) %>%
  plot_grl2023_fig3()

p5 <- p5 +
  labs(title = paste0("ID=", ids_best[5])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

format_p_value <- function(p_value) {
  if (p_value < 0.01) {
    return("(P < 0.01)")
  } else {
    return(paste0("(P = ", sprintf(p_value, fmt = "%#.2f"), ")"))
  }
}

make_annotation <- function(crpss, id, mod) {
  crpss_full <-
    crpss %>%
    filter(ID %in% id & subset %in% "full" & model %in% mod) %>%
    `$`(crpss)
  crpss_matched <-
    crpss %>%
    filter(ID %in% id & subset %in% "best_n" & model %in% mod) %>%
    `$`(crpss)
  annotation <- paste0(
    paste0("CRPSS (Full) = ", sprintf(crpss_full, fmt = "%#.2f"), "\n"),
    paste0("CRPSS (NAO-matched) = ", sprintf(crpss_matched, fmt = "%#.2f"))
  )
  annotation
}

get_y_range <- function(p) {
  yrange <- layer_scales(p)$y$range$range
  return(yrange)
}

get_y_position <- function(p, rel_pos) {
  yrange <- get_y_range(p)
  return(yrange[1] + base::diff(yrange) * rel_pos)
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
      y = yrange[1] + base::diff(yrange) * annotation_rel_pos,
      label = make_annotation(skill_scores, id, mod),
      hjust=0,
      vjust=1,
      size = annotation_size / ggplot2::.pt
    ) +
    annotate(
      geom = "text",
      x = 2005,
      y = yrange[1] + base::diff(yrange) * annotation_rel_pos,
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
gauge_stns_subset <- gauge_stns %>% filter(ID %in% ids_best)
p6 <- plot_grl2023_fig3_map(gauge_stns_subset)

p1 <- p1 + theme(plot.margin = margin(0, 0, 0, 0))
p2 <- p2 + theme(plot.margin = margin(0, 0, 0, 0.25, unit="cm"),
                 axis.title.y = element_blank())
p3 <- p3 + theme(plot.margin = margin(0, 0, 0, 0),
                 axis.title.y = element_blank())
p4 <- p4 + theme(plot.margin = margin(0, 0, 0, 0),
                 axis.title.y = element_text(size = axis_title_size))
p5 <- p5 + theme(plot.margin = margin(0, 0, 0, 0.25, unit="cm"),
                 axis.title.y = element_blank())
p6 <- p6 + theme(plot.margin = margin(0, 0, 0, 0),
                 axis.title = element_blank(),
                 axis.text = element_text(size = axis_label_size_small))

p <-
  p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3, nrow = 2, widths = c(2, 2, 2)) &
  theme(legend.position = "bottom")

p <-
  p + plot_layout(guides = "collect") &
  theme(
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(1, 'lines'),
    legend.box.spacing = unit(0, 'cm')
  )

ggsave(
  file.path(output_dir, "fig3.png"),
  plot = p, width = 6, height = 4.25, units = "in", dpi=fig_dpi
)

## ####################################################### ##
## ####################################################### ##
##
## Figure 4
##
## ####################################################### ##
## ####################################################### ##

r2 <-
  open_dataset(
    file.path(input_dir, "combined/yr2to9_lag/DJFM"),
    partitioning = c("ID", "subset")
  ) %>%
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

p <- plot_grl2023_fig4(x)
ggsave(
  file.path(output_dir, "fig4.png"),
  plot = p, width = 5, height = 5, units = "in", dpi=fig_dpi
)

## ####################################################### ##
## ####################################################### ##
##
## Figure S1
##
## ####################################################### ##
## ####################################################### ##

fit <-
  load_model_fit(config, "hindcast_Q95", aggregation_period, "DJFM") %>%
  mutate(kurtosis = kurtosis - 3)
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

p1 <- plot_grl2023_figS1a_d(skill_scores_subset, "mean")
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

p2 <- plot_grl2023_figS1a_d(skill_scores_subset, "variance")
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

p3 <- plot_grl2023_figS1a_d(skill_scores_subset, "kurtosis")
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

p4 <- plot_grl2023_figS1a_d(skill_scores_subset, "skewness")
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

p1 <- p1 +
  labs(title = "a") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))
p2 <- p2 +
  labs(title = "b") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))
p3 <- p3 +
  labs(title = "c") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))
p4 <- p4 +
  labs(title = "d") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

p <- p1 + p2 + p3 + p4 + plot_layout(ncol=2, nrow=2)

ggsave(
  file.path(output_dir, "figS1.png"),
  plot = p, width = 6, height = 6, units = "in", dpi=fig_dpi
)

## ####################################################### ##
## ####################################################### ##
##
## Figure S2
##
## ####################################################### ##
## ####################################################### ##

skill_scores_subset <-
  skill_scores %>%
  filter(
    subset %in% "full",
    period %in% aggregation_period,
    model %in% c("P", "PT")
  )

skill <-
  skill_scores_subset %>%
  group_by(ID, subset, period) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst)) %>%
  arrange(desc(!!sym(skill_measure)))
ids_best = head(skill$ID, n=5)
ids_worst = tail(skill$ID, n=5)

dataset_dir <- config$modelling[["hindcast"]]$input_dataset
predictions <- load_model_predictions(
  config, "hindcast_Q95", aggregation_period, season="DJFM"
)
predictions <-
  predictions %>%
  filter(subset %in% "full") %>%
  mutate(obs = Q_95_obs, exp = Q_95_exp) %>%
  filter(model %in% c("P", "PT")) %>%
  mutate(model = droplevels(model))

p1 <- predictions %>%
  filter(ID %in% ids_best[1]) %>%
  plot_grl2023_figS2()
p1 <- p1 +
  labs(title = paste0("ID=", ids_best[1])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p2 <- predictions %>%
  filter(ID %in% ids_best[2]) %>%
  plot_grl2023_figS2()
p2 <- p2 +
  labs(title = paste0("ID=", ids_best[2])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p3 <- predictions %>%
  filter(ID %in% ids_best[3]) %>%
  plot_grl2023_figS2()
p3 <- p3 +
  labs(title = paste0("ID=", ids_best[3])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p4 <- predictions %>%
  filter(ID %in% ids_best[4]) %>%
  plot_grl2023_figS2()
p4 <- p4 +
  labs(title = paste0("ID=", ids_best[4])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p5 <- predictions %>%
  filter(ID %in% ids_best[5]) %>%
  plot_grl2023_figS2()
p5 <- p5 +
  labs(title = paste0("ID=", ids_best[5])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p6 <- predictions %>%
  filter(ID %in% ids_worst[1]) %>%
  plot_grl2023_figS2()
p6 <- p6 +
  labs(title = paste0("ID=", ids_worst[1])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p7 <- predictions %>%
  filter(ID %in% ids_worst[2]) %>%
  plot_grl2023_figS2()
p7 <- p7 +
  labs(title = paste0("ID=", ids_worst[2])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p8 <- predictions %>%
  filter(ID %in% ids_worst[3]) %>%
  plot_grl2023_figS2()
p8 <- p8 +
  labs(title = paste0("ID=", ids_worst[3])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p9 <- predictions %>%
  filter(ID %in% ids_worst[4]) %>%
  plot_grl2023_figS2()
p9 <- p9 +
  labs(title = paste0("ID=", ids_worst[4])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

p10 <- predictions %>%
  filter(ID %in% ids_worst[5]) %>%
  plot_grl2023_figS2()
p10 <- p10 +
  labs(title = paste0("ID=", ids_worst[15])) +
  theme(plot.title = element_text(size = tag_label_size, margin = margin(0,0,1,0, unit="pt")))

make_annotation <- function(skill_scores, id) {
  crpss_p <-
    skill_scores %>%
    filter(ID %in% id & subset %in% "full" & model %in% "P") %>%
    `$`(crpss)
  crpss_pt <-
    skill_scores %>%
    filter(ID %in% id & subset %in% "best_n" & model %in% "PT") %>%
    `$`(crpss)
  annotation = paste0(
    paste0("CRPSS (P) = ", sprintf(crpss_p, fmt = "%#.2f"), "\n"),
    paste0("CRPSS (PT) = ", sprintf(crpss_pt, fmt = "%#.2f"))
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
      y = yrange[1] + base::diff(yrange) * annotation_rel_pos,
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

ggsave(
  file.path(output_dir, "figS2.png"),
  plot = p, width = 6, height = 4, units = "in", dpi=fig_dpi
)

## ####################################################### ##
## ####################################################### ##
##
## Figure S3
##
## ####################################################### ##
## ####################################################### ##

load_fcst <- function(aggregation_period) {
  fcst = read_parquet(
    file.path(
      input_dir, "nao_matching", aggregation_period,
      "DJFM/ensemble_mean_fcst.parquet"
    )
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

full_fcst <- load_fcst("yr2to9_lag")
obs <- read_parquet(
  file.path(input_dir, "meteo/yr2to9_lag/DJFM/observed.parquet")
)

ensemble_fcst <- read_parquet(
  file.path(input_dir, "meteo/yr2to9_lag/DJFM/ensemble_forecast.parquet")
)
ensemble_fcst <-
  ensemble_fcst %>%
  pivot_wider(names_from="variable", values_from="value")

nao_matched_ensemble_fcst <- read_parquet(
  file.path(input_dir, "nao_matching/yr2to9_lag/DJFM/matched_ensemble.parquet")
)

## Select n best performing members
nao_matched_ensemble_fcst_best <-
  nao_matched_ensemble_fcst %>%
  group_by(source_id, member, init_year) %>%
  group_by(init_year, variable) %>%
  slice_min(error, n = 20)

nao_matched_fcst <-
  nao_matched_ensemble_fcst_best %>%
  group_by(init_year, variable) %>%
  summarize(ens_mean = mean(value, na.rm=TRUE)) %>%
  ungroup()

## Join with observed data
nao_matched_fcst <-
  nao_matched_fcst %>%
  left_join(obs, by = c("init_year", "variable"))

## Smooth
nao_matched_fcst <-
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
nao_matched_fcst <-
  nao_matched_fcst %>%
  group_by(variable) %>%
  mutate(ens_mean_var_adj = ens_mean * sd(obs) / sd(ens_mean, na.rm=T)) %>%
  mutate(ens_mean_lag_std = ens_mean_lag / sd(ens_mean_lag, na.rm=T)) %>%
  mutate(ens_mean_lag_var_adj = ens_mean_lag_std * sd(obs))

## Plot 1 [analog of Fig 2a from Smith et al. 2020]
p1 <- plot_grl2023_figS3a(full_fcst, ensemble_fcst)
p1 <- p1 +
  labs(title = "a") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 2 [analog of Fig 2b from Smith et al. 2020]
p2 <- plot_grl2023_figS3b(full_fcst, ensemble_fcst)
p2 <- p2 +
  labs(title = "b") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 3 [analog of Fig 2c from Smith et al. 2020]
p3 <- plot_grl2023_figS3c(full_fcst, ensemble_fcst)
p3 <- p3 +
  labs(title = "c") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 4 [analog of Fig 2d from Smith et al. 2020]
p4 <- plot_grl2023_figS3d(nao_matched_fcst, ensemble_fcst)
p4 <- p4 +
  labs(title = "d") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 5 [analog of Fig 2e from Smith et al. 2020]
p5 <- plot_grl2023_figS3e(full_fcst, ensemble_fcst)
p5 <- p5 +
  labs(title = "e") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

## Plot 6 [analog of Fig 2f from Smith et al. 2020]
p6 <- plot_grl2023_figS3f(nao_matched_fcst, ensemble_fcst)
p6 <- p6 +
  labs(title = "f") +
  theme(plot.title = element_text(size = tag_label_size, face="bold"))

p1 <- p1 + theme(axis.title.y = element_text(size = axis_title_size_small))
p3 <- p3 + theme(axis.title.y = element_text(size = axis_title_size_small))
p5 <- p5 + theme(axis.title.y = element_text(size = axis_title_size_small))

p2 <- p2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
p4 <- p4 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
p6 <- p6 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())

p <- p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 2, nrow = 3) &
  theme(legend.position = "bottom")
p <- p + plot_layout(guides = "collect")

ggsave(
  file.path(output_dir, "figS3.png"),
  plot = p, width = 6, height = 6, units = "in", dpi=fig_dpi
)

## ####################################################### ##
## ####################################################### ##
##
## Figure S4
##
## ####################################################### ##
## ####################################################### ##

dataset_dir <- config$modelling[["hindcast"]]$input_dataset

skill_scores_subset <-
  skill_scores %>%
  filter(model %in% c("P", "PT") & subset %in% c("full", "best_n")) %>%
  mutate(model = droplevels(model)) %>%
  group_by(ID, subset, period) %>%
  filter(crps_ens_fcst == min(crps_ens_fcst)) %>%
  dplyr::select(-any_of(all_skill_measures)) %>% ungroup()

predictions <-
  load_model_predictions(config, "hindcast_Q95", aggregation_period, "DJFM") %>%
  mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))

predictions <- predictions %>%
  dplyr::select(-matches("Q[0-9]{2}")) %>%
  dplyr::select(-mu, -sigma, -shape, -scale, -date) %>%
  dplyr::select(-any_of(all_skill_measures))

## Take the best model for each ID/subset combination
id_year <- predictions %>% dplyr::select(ID, year) %>% distinct()

skill_scores_subset <-
  id_year %>%
  full_join(skill_scores_subset) %>% dplyr::select(-period)

skill_scores_subset <-
  skill_scores_subset %>%
  mutate(subset = ifelse(subset == "best_n", "NAO-matched ensemble", "Full ensemble"))

predictions <-
  skill_scores_subset %>%
  left_join(predictions, by = c("ID", "subset", "model", "year")) %>%
  mutate(obs = Q_95_obs, exp = Q_95_exp)

## Compute anomalies
predictions0 <-
  predictions %>%
  group_by(ID, subset) %>%
  summarize(
    ens_mean = mean(exp, na.rm = TRUE),
    obs_mean = mean(obs, na.rm = TRUE),
    ens_sd = sd(exp, na.rm = TRUE),
    obs_sd = sd(obs, na.rm = TRUE)
  )

predictions1 <-
  predictions %>%
  filter(year %in% 1985:1993) %>%
  group_by(ID, subset) %>%
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

dat1 <- predictions0 %>%
  filter(subset %in% "Full ensemble") %>%
  mutate(anom = obs_anom)
dat1 <- gauge_stns %>% left_join(dat1, by = "ID")
p1 <- plot_grl2023_figS4a_c(dat1)
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
p2 <- plot_grl2023_figS4a_c(dat2)
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
p3 <- plot_grl2023_figS4a_c(dat3)
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

p <- p1 + p2 + p3 +
  plot_layout(ncol = 3, nrow = 1) &
  theme(legend.position = "bottom")
p <- p + plot_layout(guides = "collect")

ggsave(
  file.path(output_dir, "figS4.png"),
  plot = p, width = 6, height = 4.5, units = "in", dpi=fig_dpi
)

