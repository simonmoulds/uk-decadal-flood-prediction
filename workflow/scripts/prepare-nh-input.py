#!/usr/bin/env python3

import os
import sys
import shutil
import datetime
from pathlib import Path
from tqdm import tqdm

import pyarrow
import pandas as pd
import xarray as xr
from ruamel.yaml import YAML

# Get command line arguments
inputdir = sys.argv[1]
aggr_period = sys.argv[2]
outputroot = sys.argv[3]
config = sys.argv[4]

# # For testing:
# config = 'config/config.yml'
# outputroot = 'results'
# aggr_period = 'yr2to9_lag'

outputdir = os.path.join(os.getcwd(), outputroot, 'analysis', aggr_period, 'nh-input')
try:
    os.makedirs(outputdir)
except FileExistsError:
    pass

try:
    os.makedirs(os.path.join(os.getcwd(), outputdir, 'time_series'))
except FileExistsError:
    pass

try:
    os.makedirs(os.path.join(os.getcwd(), outputdir, 'attributes'))
except FileExistsError:
    pass

rundir = os.path.join(os.getcwd(), outputroot, 'analysis', aggr_period, 'nh-output', 'runs')
try:
    os.makedirs(rundir)
except FileExistsError:
    pass

# ##################################### #
# 1 - Catchment attributes              #
# ##################################### #

yaml = YAML() #typ = 'safe')
cfg = yaml.load(Path(config))
camels_datadir = os.path.join(
    inputdir,
    # cfg['input_data_root'],
    cfg['aux_data']['camels']['subdirectory'],
    'data'
)
destdir = os.path.join(outputdir, 'attributes')

camels_attr_fs = ['CAMELS_GB_climatic_attributes.csv',
                  'CAMELS_GB_humaninfluence_attributes.csv',
                  'CAMELS_GB_hydrogeology_attributes.csv',
                  'CAMELS_GB_hydrologic_attributes.csv',
                  'CAMELS_GB_hydrometry_attributes.csv',
                  'CAMELS_GB_landcover_attributes.csv',
                  'CAMELS_GB_soil_attributes.csv',
                  'CAMELS_GB_topographic_attributes.csv']
camels_attr_fs = [os.path.join(camels_datadir, f) for f in camels_attr_fs]

for f in camels_attr_fs:
    srcfile = os.path.join(camels_datadir, f)
    shutil.copy(srcfile, destdir)

camels_meta = pd.read_csv(os.path.join(camels_datadir, 'CAMELS_GB_hydrometry_attributes.csv'))
camels_station_ids = camels_meta['gauge_id'].values

# ##################################### #
# 2 - Time series data                  #
# ##################################### #

nrfa_meta = pd.read_parquet(os.path.join(outputroot, 'nrfa-metadata.parquet'))
nrfa_meta = nrfa_meta.astype({'id' : int})
nrfa_meta = nrfa_meta.set_index('id')

# This generated from rule build_catchment_dataset
x = pd.read_parquet(os.path.join(outputroot, 'analysis', aggr_period, 'input'))
x = x.astype({'year' : int, 'lead_time' : int, 'ID' : int})
station_ids = list(x['ID'].unique())
# Ensure the selected stations are included in the CAMELS dataset,
# because we currently take the metadata from this data source
station_ids = [stn for stn in station_ids if stn in camels_station_ids]
n_stations = len(station_ids)

for i in tqdm(range(n_stations)):
    stn = station_ids[i]
    xx = x.loc[(x['ID'] == stn) & (x['subset'] == 'best_n')].copy(deep = True)
    fcst_start_years = list(xx['year'] + xx['lead_time'] - 1)
    fcst_start_dates = pd.to_datetime([datetime.datetime(yr, 12, 1, 0, 0) for yr in fcst_start_years])
    xx.loc[:,'date'] = fcst_start_dates
    xx = xx.set_index('date')
    xx = xx[['Q_max', 'Q_mean', 'Q_05', 'Q_50', 'Q_95', 'nao', 'ea', 'amv', 'european_precip', 'uk_temp']]
    xx = xx.rename(columns = {
        'Q_max' : 'Qmax',
        'Q_mean' : 'Qmean',
        'Q_05' : 'Q05',
        'Q_50' : 'Q50',
        'Q_95' : 'Q95',
        'nao' : 'NAO',
        'ea' : 'EA',
        'amv' : 'AMV',
        'european_precip' : 'P',
        'uk_temp' : 'T'
    })
    # Normalize by catchment area (take step by step for clarity)
    area = nrfa_meta.loc[stn, 'catchment_area'] * 1000 * 1000 # km2 to m2
    xx['Q95'] *= (24 * 60 * 60)                              # m3 s-1 to m3 d-1
    xx['Q95'] /= area                                        # m3 d-1 to m d-1
    xx['Q95'] *= 1000                                        # m d-1 to mm d-1
    xarr = xr.Dataset.from_dataframe(xx)
    nc_filename = os.path.join(outputdir, 'time_series', str(stn) + '.nc')
    xarr.to_netcdf(nc_filename)

# ##################################### #
# 3 - NH configuration                  #
# ##################################### #
climatic_attributes = ['p_mean', 'pet_mean', 'aridity', 'p_seasonality', 'frac_snow', 'high_prec_freq', 'high_prec_dur', 'low_prec_freq', 'low_prec_dur']
human_influence_attributes = []
hydrogeology_attributes = ['inter_high_perc', 'inter_mod_perc', 'inter_low_perc', 'frac_high_perc', 'frac_mod_perc', 'frac_low_perc', 'no_gw_perc', 'low_nsig_perc', 'nsig_low_perc']
hydrologic_attributes = []
landcover_attributes = ['dwood_perc', 'ewood_perc', 'urban_perc']
soil_attributes = ['sand_perc', 'silt_perc', 'clay_perc', 'porosity_hypres', 'conductivity_hypres', 'soil_depth_pelletier_50']
topographic_attributes = ['gauge_lat', 'gauge_lon', 'gauge_elev', 'area', 'elev_10', 'elev_50', 'elev_90']

static_attributes = climatic_attributes + \
    human_influence_attributes + \
    hydrogeology_attributes + \
    hydrologic_attributes + \
    landcover_attributes + \
    soil_attributes + \
    topographic_attributes

dynamic_inputs = ['P', 'T', 'EA', 'AMV']
# static_attributes = ['gauge_lat', 'gauge_lon', 'p_mean', 'pet_mean', 'area'] #, 'q_mean']#, 'frac_snow', 'high_prec_freq', 'high_prec_dur', 'low_prec_freq', 'low_prec_dur']

yaml = YAML() #typ = 'safe')
cfg = yaml.load(Path('resources/nh-config-template.yml'))
cfg['experiment_name'] = 'cudalstm_' + str(n_stations) + '_basins_' + str(aggr_period)
cfg['run_dir'] = rundir
cfg['train_basin_file'] = os.path.join(outputdir, "basins.txt")
cfg['validation_basin_file'] = os.path.join(outputdir, "basins.txt")
cfg['test_basin_file'] = os.path.join(outputdir, "basins.txt")
cfg['train_start_date'] = ["01/12/1961", "01/12/1992"] #"01/12/1961"
cfg['train_end_date'] = ["01/12/1990", "01/12/2006"] #"01/12/1990", #"01/12/2006"
cfg['validation_start_date'] = "01/12/1991" #"01/12/1961"
cfg['validation_end_date'] = "01/12/1991"
cfg['test_start_date'] = "01/12/1961"
cfg['test_end_date']= "01/12/2006"
cfg['device'] = "cpu" # "cuda:0"
cfg['validate_every'] = int(3)
cfg['validate_n_random_basins'] = int(1)
cfg['metrics'] = ["MSE"]
cfg['save_validation_results'] = True
cfg['model'] = "cudalstm" #"cudalstm" # "ealstm"
cfg['head'] = "regression"
cfg['output_activation'] = "linear"
cfg['hidden_size'] = int(64) #int(128) #int(20) #int(64) #int(20)
cfg['initial_forget_bias'] = int(3)
cfg['output_dropout'] = 0.4
cfg['optimizer'] = "Adam"
cfg['loss'] = "MSE"
cfg['learning_rate'] = {0: 1e-2, 30: 5e-3, 40: 1e-3}
cfg['batch_size'] = int(128) #int(256)
cfg['epochs'] = int(50) #int(50) #int(500) #int(150) #int(75)
cfg['clip_gradient_norm'] = int(1)
cfg['predict_last_n'] = int(1)
cfg['seq_length'] = int(4) #int(4) #int(8)
cfg['num_workers'] = int(8)
cfg['log_interval'] = int(5)
cfg['log_tensorboard'] = True
cfg['log_n_figures'] = int(0)
cfg['save_weights_every'] = int(1)
cfg['dataset'] = "generic"
cfg['data_dir'] = outputdir
cfg['dynamic_inputs'] = dynamic_inputs
cfg['target_variables'] = ["Q95"]
cfg['clip_targets_to_zero'] = ["Q95"]
cfg['static_attributes'] = static_attributes
conf_filename = os.path.join(outputdir, 'basins.yml')
with open(conf_filename, 'wb') as f:
    yaml.dump(cfg, f)

# ##################################### #
# 4 - Basins list                       #
# ##################################### #

basin_filename = os.path.join(outputdir, 'basins.txt')
with open(basin_filename, 'w') as f:
    for stn in station_ids:
        f.write(str(stn))
        f.write('\n')
