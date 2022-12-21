
## Output of prepare_input_data

## x <- read_parquet('results_wrong/analysis/yr2to5_lag/observed.parquet')
## y <- read_parquet('results/analysis/yr2to5_lag/observed.parquet')
## all.equal(x, y) # TRUE

x <- read_parquet('results_wrong/analysis/yr2to5_lag/ensemble_forecast.parquet')
y <- read_parquet('results/analysis/yr2to5_lag/ensemble_forecast.parquet')
all.equal(x, y) # FALSE - this is the reason for the discrepancy

dim(x)
dim(y)

table(x$source_id)
table(y$source_id)

table(x$init_year)
table(y$init_year)

table(x$variable)
table(y$variable)
## Difference in european_precip_antecedent
## This is concerning because european_precip_antecedent is not used in exp 1

## The variables we use in the analysis are european_precip, amv, nao
x_nao <- x %>% filter(variable %in% "nao")
y_nao <- y %>% filter(variable %in% "nao")
all.equal(x_nao, y_nao) # TRUE

x_european_precip <- x %>% filter(variable %in% "european_precip")
y_european_precip <- y %>% filter(variable %in% "european_precip")
all.equal(x_european_precip, y_european_precip)

x_amv <- x %>% filter(variable %in% "amv")
y_amv <- y %>% filter(variable %in% "amv")
all.equal(x_amv, y_amv)

x_ea <- x %>% filter(variable %in% "ea")
y_ea <- y %>% filter(variable %in% "ea")
all.equal(x_ea, y_ea)

x_uk_precip <- x %>% filter(variable %in% "uk_precip")
y_uk_precip <- y %>% filter(variable %in% "uk_precip")
all.equal(x_uk_precip, y_uk_precip)

x_uk_temp <- x %>% filter(variable %in% "uk_temp")
y_uk_temp <- y %>% filter(variable %in% "uk_temp")
all.equal(x_uk_temp, y_uk_temp)

## x_european_precip_antecedent <- x %>% filter(variable %in% "european_precip_antecedent")
## y_european_precip_antecedent <- y %>% filter(variable %in% "european_precip_antecedent")
## all.equal(x_european_precip_antecedent, y_european_precip_antecedent)

x <- read_parquet('results_wrong/analysis/yr2to9_lag/ensemble_mean_fcst.parquet')
y <- read_parquet('results/analysis/yr2to9_lag/ensemble_mean_fcst.parquet')
all.equal(x, y)

## The variables we use in the analysis are european_precip, amv, nao
x_nao <- x %>% filter(variable %in% "nao")
y_nao <- y %>% filter(variable %in% "nao")
all.equal(x_nao, y_nao) # TRUE

x_european_precip <- x %>% filter(variable %in% "european_precip")
y_european_precip <- y %>% filter(variable %in% "european_precip")
all.equal(x_european_precip, y_european_precip)

x_amv <- x %>% filter(variable %in% "amv")
y_amv <- y %>% filter(variable %in% "amv")
all.equal(x_amv, y_amv)

x <- read_parquet('results_wrong/analysis/yr2to9_lag/matched_ensemble.parquet')
y <- read_parquet('results/analysis/yr2to9_lag/matched_ensemble.parquet')
all.equal(x, y)

## The variables we use in the analysis are european_precip, amv, nao
x_nao <- x %>% filter(variable %in% "nao")
y_nao <- y %>% filter(variable %in% "nao")
all.equal(x_nao, y_nao) # TRUE

x_european_precip <- x %>% filter(variable %in% "european_precip")
y_european_precip <- y %>% filter(variable %in% "european_precip")
all.equal(x_european_precip, y_european_precip)

x_amv <- x %>% filter(variable %in% "amv")
y_amv <- y %>% filter(variable %in% "amv")
all.equal(x_amv, y_amv)

x <- read_parquet('results_wrong/analysis/yr2to9_lag/matched_ensemble_error.parquet')
y <- read_parquet('results/analysis/yr2to9_lag/matched_ensemble_error.parquet')
all.equal(x, y)

## build-catchment-dataset.R

x <- open_dataset("results/analysis/yr2to5_lag/input") %>% collect()
y <- open_dataset("results_wrong/analysis/yr2to5_lag/input") %>% collect()
all.equal(x, y)

x_7001 <- x %>% filter(ID %in% 7001)
y_7001 <- y %>% filter(ID %in% 7001)
all.equal(x_7001$european_precip, y_7001$european_precip)
all.equal(x_7001$nao, y_7001$nao)
all.equal(x_7001$amv, y_7001$amv)

x <- x %>% dplyr::select(-european_precip_antecedent)
y <- y %>% dplyr::select(-european_precip_antecedent)
all.equal(x, y) # missing_pct - big difference because of european_precip_antecedent
