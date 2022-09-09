#!/usr/bin/env python3

import os
import sys
import shutil
import pickle
import numpy as np
import pandas as pd
import xarray as xr
import datetime
from pathlib import Path
from ruamel.yaml import YAML

from sklearn.metrics import mean_squared_error
# import matplotlib.pyplot as plt
# import torch
# from neuralhydrology.evaluation import metrics
# from neuralhydrology.nh_run import start_run, eval_run
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torch.utils.data import random_split
from torch import Tensor
from torch.nn import Linear
from torch.nn import Sigmoid
from torch.nn import Module
from torch.optim import SGD
from torch.nn import MSELoss
from torch.nn.init import xavier_uniform_

# Some helpers from neuralhydrology
from neuralhydrology.utils.config import Config
from neuralhydrology.datautils.utils import load_basin_file
from neuralhydrology.datasetzoo.genericdataset import load_timeseries, load_attributes

# Dataset definition
class DCPDataset(Dataset):
    def __init__(self,
                 cfg: Config,
                 is_train: bool,
                 period: str):
        self.cfg = cfg
        self.period = period
        data_dir = self.cfg.data_dir
        basins = Path(getattr(self.cfg, f'{period}_basin_file'))
        with basins.open('r') as f:
            basins = f.read().splitlines()

        start_dates = getattr(self.cfg, f'{period}_start_date')
        end_dates = getattr(self.cfg, f'{period}_end_date')
        if not isinstance(start_dates, list):
            start_dates = [start_dates]
            end_dates = [end_dates]

        # Read dynamic input data
        tss = []
        scaler = {}
        for basin in basins:
            ts = load_timeseries(data_dir, basin) # TESTING
            ts = ts.reset_index()
            ts['gauge_id'] = basin
            ts['year'] = [tm.year for tm in ts['date']]
            ts_subset = []
            for i in range(len(start_dates)):
                start = start_dates[i]
                end = end_dates[i]
                ts_subset.append(ts[(ts['date'] >= start) & (ts['date'] <= end)])
            ts = pd.concat(ts_subset)
            ts = ts.drop('date', axis = 1)
            scaler[str(basin)] = {'center' : ts.mean(), 'scale' : ts.std()}
            tss.append(ts)
        ts = pd.concat(tss)
        self.index = ts[['gauge_id', 'year']]
        ts = ts.set_index(['gauge_id'])

        # Read attributes
        attr = load_attributes(data_dir, basins)
        attr = attr.reset_index()
        attr = attr.set_index('gauge_id')

        # Merge dynamic and static inputs
        dyn_inputs = self.cfg.dynamic_inputs
        static_inputs = self.cfg.static_attributes
        inputs = pd.merge(
            ts[self.cfg.dynamic_inputs],
            attr[self.cfg.static_attributes],
            how = 'left',
            on = 'gauge_id'
        ).values
        self.y = ts[self.cfg.target_variables].values
        self.X = inputs
        # Ensure floats
        self.y = self.y.astype('float32')
        self.X = self.X.astype('float32')
        # Validate (check for NaN)
        target_idx = np.array([all(~np.isnan(self.y[i,:])) for i in range(len(self.y))])
        input_idx = np.array([all(~np.isnan(self.X[i,:])) for i in range(len(self.X))])
        nan_idx = target_idx & input_idx
        self.y = self.y[nan_idx,:]
        self.X = self.X[nan_idx,:]

    def __len__(self):
        return(len(self.X))

    def __getitem__(self, idx):
        return [self.X[idx], self.y[idx]]


class MLP(Module):
    def __init__(self, n_inputs):
        super(MLP, self).__init__()
        # Input to first hidden layer
        self.hidden1 = Linear(n_inputs, 10)
        xavier_uniform_(self.hidden1.weight)
        self.act1 = Sigmoid()
        # Input to second hidden layer
        self.hidden2 = Linear(10, 8)
        xavier_uniform_(self.hidden1.weight)
        self.act2 = Sigmoid()
        # Third hidden layer and output
        self.hidden3 = Linear(8, 1)
        xavier_uniform_(self.hidden3.weight)

    def forward(self, X):
        # Input to first hidden layer
        X = self.hidden1(X)
        X = self.act1(X)
        # Second hidden layer
        X = self.hidden2(X)
        X = self.act2(X)
        # Third hidden layer and output
        X = self.hidden3(X)
        return X

def prepare_mlp_data(cfg: Config):
    train = DCPDataset(cfg, is_train = True, period = 'train')
    test = DCPDataset(cfg, is_train = False, period = 'test')
    # Prepare data loaders
    train_dl = DataLoader(train, batch_size = cfg.batch_size, shuffle = True)
    test_dl = DataLoader(test, batch_size = 1024, shuffle = False)
    return train_dl, test_dl

def prepare_xgb_data(cfg: Config):
    train = DCPDataset(cfg, is_train = True, period = 'train')
    test = DCPDataset(cfg, is_train = False, period = 'test')
    return train, test

def train_model(train_dl, model):
    criterion = MSELoss()
    # lr is learning rate, which could change by epoch (see neuralhydrology.training)
    optimizer = SGD(model.parameters(), lr = 0.01, momentum = 0.9)
    for epoch in range(100):
        for i, (inputs, targets) in enumerate(train_dl):
            # Clear the gradients
            optimizer.zero_grad()
            # Compute the model output
            yhat = model(inputs)
            # Calculate loss
            loss = criterion(yhat, targets)
            # Credit assignment
            loss.backward()
            # Update model weights
            optimizer.step()

def evaluate_model(test_dl, model):
    predictions, actuals = list(), list()
    for i, (inputs, targets) in enumerate(test_dl):
        # Evaluate the model on the test set
        yhat = model(inputs)
        # Retrieve numpy array
        yhat = yhat.detach().numpy()
        actual = targets.numpy()
        actual = actual.reshape((len(actual), 1))
        # Store
        predictions.append(yhat)
        actuals.append(actual)
    predictions, actuals = vstack(predictions), vstack(actuals)
    # Calculate MSE
    mse = mean_squared_error(actuals, predictions)
    return mse

def predict(row, model):
    # Convert row to data
    row = Tensor([row])
    # Make prediction
    yhat = model(row)
    # Retrieve numpy array
    yhat = yhat.detach().numpy()
    return yhat

# # Get command line arguments
# nh_config = sys.argv[1]
# outputdir = sys.argv[2]

# Testing
nh_config = 'results/analysis/yr2to9_lag/nh-input/basins.yml'
outputdir = '.'

# # Load configuration
# cfg = Config(Path(nh_config))
# train_dl, test_dl = prepare_data(cfg)
# model = MLP(38)                  # not sure where 4 comes from?
# train_model(train_dl, model)
# acc = evaluate_model(test_dl, model) # NaN problem

# For handling NaN, see validate_samples in basedataset.py
# For now, just remove all data points that have an NaN

from xgboost import XGBRegressor
from sklearn.metrics import mean_absolute_error, r2_score

cfg = Config(Path(nh_config))
cfg.update_config(
    {'train_start_date' : '01/12/1961',
     # 'train_end_date' : '01/12/2006',
     'train_end_date' : '01/12/1979',
     'test_start_date' : '01/12/1986',
     'test_end_date' : '01/12/2006',
     # 'static_attributes' : [],
     'dynamic_inputs' : ['P', 'T', 'EA', 'AMV']
     })

train, test = prepare_xgb_data(cfg)

model = XGBRegressor(n_estimators = 500)
model.set_params(early_stopping_rounds = 5) #, eval_metric = [r2_score])
model.fit(train.X, train.y, verbose = False, eval_set = [(test.X, test.y)])
predictions = model.predict(test.X)
# print("Mean Absolute Error : " + str(r2_score(predictions, test.y)))
output = test.index
output['obs'] = test.y
output['exp'] = predictions
output.to_csv('results/test.csv', index = False)

# if torch.cuda.is_available():
#     start_run(config_file = Path(nh_config))
# else:
#     # raise OSError("No GPU available!")
#     start_run(config_file = Path(nh_config), gpu = -1)

# # years = [y for y in range(1961, 2007)]
# years = [y for y in range(1961, 2007)]
# padding = 20
# seq_length = int(cfg['seq_length'])
# n_test = 1
# n_validate = 1
# n_chains = len(years) - seq_length - n_test - n_validate - padding
# training_start = "01/12/" + str(years[0])
# rundir = cfg['run_dir']
# forward_chain_rundir = os.path.join(rundir, 'forward_chain_run')
# if os.path.isdir(forward_chain_rundir):
#     shutil.rmtree(forward_chain_rundir)

# for i in range(n_chains):
#     cfg['train_start_date'] = training_start
#     cfg['train_end_date'] = "01/12/" + str(years[(i + seq_length + padding)])
#     cfg['validation_start_date'] = "01/12/" + str(years[(i + seq_length + padding + n_validate)])
#     cfg['validation_end_date'] = cfg['validation_start_date']
#     cfg['test_start_date'] = "01/12/" + str(years[(i + seq_length + padding + n_validate + n_test)])
#     cfg['test_end_date'] = cfg['test_start_date']

#     # Write unique experiment name
#     cfg['run_dir'] = forward_chain_rundir
#     cfg['experiment_name'] = 'chain_' + str(i)

#     # TODO write config to temporary yaml file
#     tmp_conf_filename = os.path.join('/tmp', 'chain_' + str(i) + '.yml')
#     with open(tmp_conf_filename, 'wb') as f:
#         yaml.dump(cfg, f)

#     if torch.cuda.is_available():
#         start_run(config_file = Path(tmp_conf_filename))
#     else:
#         # raise OSError("No GPU available!")
#         start_run(config_file = Path(tmp_conf_filename), gpu = -1)

#     # Neural Hydrology assigns a unique name to the output run
#     # directory. We find this name by identifying the most recent
#     # directory, then evaluate the model output.
#     base_run_dir = cfg['run_dir']
#     experiment_name = cfg['experiment_name']
#     run_dirs = [os.path.join(base_run_dir, d) for d in os.listdir(base_run_dir) if  os.path.isdir(os.path.join(base_run_dir, d)) & d.startswith(experiment_name)]
#     # if len(run_dirs) > 0:
#     run_dirs.sort(key = lambda x: os.path.getmtime(x))
#     run_dir = run_dirs[-1]

#     # # Create a symbolic link
#     # src = os.path.join(os.getcwd(), run_dir)
#     # dst = os.path.join(base_run_dir, 'latest')
#     # try:
#     #     os.symlink(src, dst)
#     # except FileExistsError:
#     #     if os.path.islink(dst):
#     #         os.unlink(dst)
#     #     elif os.path.isdir(dst):
#     #         shutil.rmtree(dst)
#     #     os.symlink(src, dst)

#     eval_run(Path(run_dir), period = 'test')

#     # Unpack time series output
#     try:
#         os.makedirs(outputdir)
#     except FileExistsError:
#         pass

#     n_epochs = cfg['epochs']
#     epoch_dirname = 'model_epoch' + str(n_epochs).zfill(3)
#     with open(os.path.join(run_dir, 'test', epoch_dirname, 'test_results.p'), 'rb') as fp:
#         results = pickle.load(fp)

#     station_ids = list(results.keys())
#     n_stations = len(station_ids)
#     freq = '1AS-DEC'

#     for j in range(n_stations):
#         stn = station_ids[j]
#         df = results[stn][freq]['xr'].to_dataframe()
#         df.insert(0, "ID", stn)
#         csv_filename = 'chain_' + str(i) + '_' + str(stn) + '.csv'
#         df.to_csv(os.path.join(outputdir, csv_filename))
