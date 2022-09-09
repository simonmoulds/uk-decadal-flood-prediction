#!/usr/bin/env python3

import os
import sys
import shutil
import pickle
import pandas as pd
import xarray as xr
import datetime
from pathlib import Path
from ruamel.yaml import YAML

# import matplotlib.pyplot as plt
import torch
from neuralhydrology.evaluation import metrics
from neuralhydrology.nh_run import start_run, eval_run

# Get command line arguments
nh_config = sys.argv[1]
# aggr_period = sys.argv[2]
outputdir = sys.argv[2]

# Load configuration
yaml = YAML() #typ = 'safe')
cfg = yaml.load(Path(nh_config))

# if torch.cuda.is_available():
#     start_run(config_file = Path(nh_config))
# else:
#     # raise OSError("No GPU available!")
#     start_run(config_file = Path(nh_config), gpu = -1)

# years = [y for y in range(1961, 2007)]
years = [y for y in range(1961, 2007)]
padding = 20
seq_length = int(cfg['seq_length'])
n_test = 1
n_validate = 1
n_chains = len(years) - seq_length - n_test - n_validate - padding
training_start = "01/12/" + str(years[0])
rundir = cfg['run_dir']
forward_chain_rundir = os.path.join(rundir, 'forward_chain_run')
if os.path.isdir(forward_chain_rundir):
    shutil.rmtree(forward_chain_rundir)

for i in range(n_chains):
    cfg['train_start_date'] = training_start
    cfg['train_end_date'] = "01/12/" + str(years[(i + seq_length + padding)])
    cfg['validation_start_date'] = "01/12/" + str(years[(i + seq_length + padding + n_validate)])
    cfg['validation_end_date'] = cfg['validation_start_date']
    cfg['test_start_date'] = "01/12/" + str(years[(i + seq_length + padding + n_validate + n_test)])
    cfg['test_end_date'] = cfg['test_start_date']

    # Write unique experiment name
    cfg['run_dir'] = forward_chain_rundir
    cfg['experiment_name'] = 'chain_' + str(i)

    # TODO write config to temporary yaml file
    tmp_conf_filename = os.path.join('/tmp', 'chain_' + str(i) + '.yml')
    with open(tmp_conf_filename, 'wb') as f:
        yaml.dump(cfg, f)

    if torch.cuda.is_available():
        start_run(config_file = Path(tmp_conf_filename))
    else:
        # raise OSError("No GPU available!")
        start_run(config_file = Path(tmp_conf_filename), gpu = -1)

    # Neural Hydrology assigns a unique name to the output run
    # directory. We find this name by identifying the most recent
    # directory, then evaluate the model output.
    base_run_dir = cfg['run_dir']
    experiment_name = cfg['experiment_name']
    run_dirs = [os.path.join(base_run_dir, d) for d in os.listdir(base_run_dir) if  os.path.isdir(os.path.join(base_run_dir, d)) & d.startswith(experiment_name)]
    # if len(run_dirs) > 0:
    run_dirs.sort(key = lambda x: os.path.getmtime(x))
    run_dir = run_dirs[-1]

    # # Create a symbolic link
    # src = os.path.join(os.getcwd(), run_dir)
    # dst = os.path.join(base_run_dir, 'latest')
    # try:
    #     os.symlink(src, dst)
    # except FileExistsError:
    #     if os.path.islink(dst):
    #         os.unlink(dst)
    #     elif os.path.isdir(dst):
    #         shutil.rmtree(dst)
    #     os.symlink(src, dst)

    eval_run(Path(run_dir), period = 'test')

    # Unpack time series output
    try:
        os.makedirs(outputdir)
    except FileExistsError:
        pass

    n_epochs = cfg['epochs']
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
        csv_filename = 'chain_' + str(i) + '_' + str(stn) + '.csv'
        df.to_csv(os.path.join(outputdir, csv_filename))
