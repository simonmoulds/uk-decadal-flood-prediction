#!/usr/bin/env python3

import os
import sys
import shutil
import pickle
import pandas as pd
import xarray as xr
import pyarrow as pa
import pyarrow.parquet as pq
import datetime
from pathlib import Path
from ruamel.yaml import YAML

# import matplotlib.pyplot as plt
import torch
from neuralhydrology.evaluation import metrics
from neuralhydrology.nh_run import start_run, eval_run

# # TESTING
# config = 'config/config.yml'
# nh_config = 'results/exp2/analysis/yr2/nh-input/basins.yml'
# aggregation_period = 'yr2'
# outputdir = 'results/exp2/analysis/hindcast/lstm'

# Get command line arguments
config = sys.argv[1]
nh_config = sys.argv[2]
aggregation_period = sys.argv[3]
outputdir = sys.argv[4]

# Load neuralhydrology configuration
yaml = YAML() #typ = 'safe')
nh_cfg = yaml.load(Path(nh_config))
cfg = yaml.load(Path(config))
aggr_period_info = [d for d in cfg['aggregation_period'] if d['name'] == aggregation_period]
if len(aggr_period_info) != 1:
    ValueError
else:
    aggr_period_info = dict(aggr_period_info[0])
lead_time = str(aggr_period_info['lead_time'])
lead_time = [int(i) for i in lead_time.split(':')]
if len(lead_time) > 1:
    lead_time = [i for i in range(lead_time[0], lead_time[1] + 1)]

# Where to put predictions
prediction_outputdir = os.path.join(outputdir, aggregation_period)

# Set up while-loop
leave_out = len(lead_time)
train_year_start = 1961
train_year_end = 1979
test_year_start = train_year_end + leave_out
max_test_year_start = 2006
output_list = []
while test_year_start <= max_test_year_start:
    nh_cfg['train_start_date'] = "01/12/" + str(train_year_start)
    nh_cfg['train_end_date'] = "01/12/" + str(train_year_end)
    nh_cfg['validation_start_date'] = "01/12/" + str(test_year_start)
    nh_cfg['validation_end_date'] = nh_cfg['validation_start_date']
    nh_cfg['test_start_date'] = "01/12/" + str(test_year_start)
    nh_cfg['test_end_date'] = nh_cfg['test_start_date']

    # Write unique experiment name
    nh_cfg['experiment_name'] = 'prediction_' + str(test_year_start)

    # Write config to temporary yaml file
    tmp_conf_filename = os.path.join('/tmp', 'config_' + str(test_year_start) + '.yml')
    with open(tmp_conf_filename, 'wb') as f:
        yaml.dump(nh_cfg, f)

    if torch.cuda.is_available():
        # Run on GPU
        start_run(config_file = Path(tmp_conf_filename))
    else:
        # Run on CPU instead
        start_run(config_file = Path(tmp_conf_filename), gpu = -1)

    # Neural Hydrology assigns a unique name to the output run
    # directory. We find this name by identifying the most recent
    # directory, then evaluate the model output.
    base_run_dir = nh_cfg['run_dir']
    experiment_name = nh_cfg['experiment_name']
    run_dirs = [os.path.join(base_run_dir, d) for d in os.listdir(base_run_dir) if  os.path.isdir(os.path.join(base_run_dir, d)) & d.startswith(experiment_name)]
    run_dirs.sort(key = lambda x: os.path.getmtime(x))
    run_dir = run_dirs[-1]
    eval_run(Path(run_dir), period = 'test')

    # Unpack time series output
    try:
        os.makedirs(outputdir)
    except FileExistsError:
        pass

    n_epochs = nh_cfg['epochs']
    epoch_dirname = 'model_epoch' + str(n_epochs).zfill(3)
    with open(os.path.join(run_dir, 'test', epoch_dirname, 'test_results.p'), 'rb') as fp:
        results = pickle.load(fp)

    station_ids = list(results.keys())
    n_stations = len(station_ids)
    freq = '1AS-DEC'
    for j in range(n_stations):
        stn = station_ids[j]
        df = results[stn][freq]['xr'].to_dataframe()
        df.insert(0, "ID", stn)
        df['year'] = test_year_start
        df = df.reset_index()
        df = df.drop('time_step', axis = 1)
        rowdata_df = pa.Table.from_pandas(df, preserve_index=False)
        pq.write_to_dataset(
            rowdata_df,
            root_path = os.path.join(prediction_outputdir, 'prediction'),
            partition_cols = ['ID', 'date']
        )
        # csv_filename = 'prediction_' + str(stn) + '_' + str(test_year_start) + '.csv'
        # df.to_csv(os.path.join(prediction_outputdir, 'tmp', csv_filename))

    # Tidy up tmp directory
    os.remove(tmp_conf_filename)

    # Update train/test years for next prediction year
    train_year_end += 1
    test_year_start += 1

