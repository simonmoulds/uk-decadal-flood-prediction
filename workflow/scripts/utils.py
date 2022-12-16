#!/usr/bin/env python3

import os
import subprocess
import glob
import re
import datetime
import pandas as pd
import numpy as np
import netCDF4
import iris
import iris.pandas
import xarray
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from functools import reduce
from calendar import monthrange
from tqdm import tqdm

from esmvalcore import *

# TESTING
# from workflow.scripts.esmvalcore import *

VALID_Y_NAMES = ['latitude', 'lat']
VALID_X_NAMES = ['longitude', 'lon']
VALID_TIME_NAMES = ['t', 'time']


def _parse_filepath(fpath):
    # Parse CMIP-formatted filename to get key id variables
    fname = os.path.basename(fpath)
    fname = fname.split('_')
    project = fname[0]
    model = fname[1]
    mip = fname[2]
    exp = fname[3]
    if project == 'CMIP5':
        ens = fname[4]
        init_year = exp.replace('decadal', '')
    elif project == 'CMIP6':
        sub_exp = fname[4]
        sub_exp = sub_exp.split('-')
        ens = sub_exp[1]
        init_year = sub_exp[0].replace('s','')
        # exp = sub_exp
    out = {
        'project': project,
        'model': model,
        'mip': mip,
        'experiment': exp,
        'ensemble': ens,
        'init_year': int(init_year)
    }
    return out


def _compute_mean_season(ds,
                         varname,
                         init_year,
                         season=[12, 1, 2, 3],
                         start_year=2,
                         end_year=9):

    # Compute the mean value for DJFM season.
    #
    # Args:
    #   ds        : xarray dataset
    #   varname   : string. Variable name.
    #   init_year : int. Initialization year of `ds`
    #   season    : list.
    #   start_year: int.
    #   end_year  : int.
    #
    # Returns:
    #   xarray dataset.

    n_months = len(season)
    if not all([month in range(1, 13) for month in season]):
        raise ValueError("`season` contains invalid months")

    # Sometimes Iris "fixes" a non-monotonic time
    # dimension by adding another dimension
    # called 'dim_0' - account for this here
    dimnames = list(ds.dims.keys())
    time_dimname = 'time'
    if 'time' not in dimnames:
        if len(dimnames) == 1:
            time_dimname = dimnames[0]
        else:
            raise ValueError("Ambiguous time dimension")

    # This is now necessary because we have included additional months
    ds = ds.where(ds.time.dt.month.isin(season), drop=True)

    ds['counter'] = ((time_dimname,), [1, ] * len(ds.time))
    ds = ds.groupby('season_year').sum(time_dimname)
    # Limit dataset to complete seasons
    complete_index = ds['counter'].values == n_months
    ds = ds.sel(season_year=complete_index)
    # Divide by number of months to get mean value
    ds[varname] = ds[varname] / n_months
    # Add lead time
    # [N.B. season_year (added automatically by Iris)
    # defines the year of final month in the season]
    ds['lead_time'] = (('season_year',), ds.season_year.values - int(init_year))

    def _is_yr2to9(year):
        return (year >= start_year) & (year <= end_year)

    ds = ds.sel(season_year=_is_yr2to9(ds['lead_time']))
    return ds[varname]


def _get_variable_name(ds):
    # Get variable name from dataset.
    varnames = [k for k, _ in ds.data_vars.items()]
    if 'nao' in varnames:
        return 'nao'
    elif 'ea' in varnames:
        return 'ea'
    elif 'amv' in varnames:
        return 'amv'
    elif 'uk_temp' in varnames:
        return 'uk_temp'
    elif 'european_precip' in varnames:
        return 'european_precip'
    elif 'uk_precip' in varnames:
        return 'uk_precip'
    elif 'uk_precip_field' in varnames:
        return 'uk_precip_field'
    elif 'precip_field' in varnames:
        return 'precip_field'
    elif 'temp_field' in varnames:
        return 'temp_field'


def _index_xarray_to_dataframe(x, varname, metadata):
    df = x.to_dataframe()
    df = df.rename(
        {varname: 'value'}, axis="columns"
    ).reset_index()
    df['project'] = metadata['project']
    df['source_id'] = metadata['model']
    df['mip'] = metadata['mip']
    df['experiment'] = metadata['experiment']
    df['member'] = metadata['ensemble']
    df['init_year'] = metadata['init_year']
    df['variable'] = varname
    df = df.sort_values(['init_year', 'season_year'])
    cols = [
        'project', 'source_id', 'mip', 'experiment',
        'member', 'init_year', 'season_year', 'variable', 'value'
    ]
    df = df[cols]
    return df


def _index_xarray_to_dataframe_new(x, varname, metadata):
    df = x.to_dataframe()
    df = df.rename(
        {varname: 'value'}, axis="columns"
    ).reset_index()
    df['project'] = metadata['project']
    df['source_id'] = metadata['model']
    df['mip'] = metadata['mip']
    df['experiment'] = metadata['experiment']
    df['member'] = metadata['ensemble']
    df['init_year'] = metadata['init_year']
    df['variable'] = varname
    # Fix time values
    df['year'] = [tm.year for tm in df['time']]
    df['month'] = [tm.month for tm in df['time']]
    df = df.sort_values(['init_year', 'year', 'month']) #, 'season_year'])
    cols = [
        'year', 'month', # TODO change to ID
        'project', 'source_id', 'mip', 'experiment',
        'member', 'init_year', 'variable', 'value'
    ]
    # try:
    df = df[cols]
    # except:
    #     # print(df)
    #     # print(varname)
    #     # print(x)
    #     raise KeyError
    return df


def _spatial_xarray_to_dataframe_new(x, varname, metadata):
    df = x.to_dataframe()
    df = df.reset_index()
    keep_cols = ['ID', 'time', varname]
    # keep_cols = ['ID', 'season_year', varname]
    drop_cols = [nm for nm in df if nm not in keep_cols]
    # keep_cols = ['time', 'ID', 'season_year', varname]
    # drop_cols = [nm for nm in df if nm not in keep_cols]
    df = df.drop(drop_cols, axis=1)
    df = df.set_index(keep_cols)
    df = df.reset_index()
    df = df.rename(
        {varname: 'value'}, axis="columns"
    ).reset_index()
    df['project'] = metadata['project']
    df['source_id'] = metadata['model']
    df['mip'] = metadata['mip']
    df['experiment'] = metadata['experiment']
    df['member'] = metadata['ensemble']
    df['init_year'] = metadata['init_year']
    df['variable'] = varname
    # Convert time into year month
    df['year'] = [tm.year for tm in df['time']]
    df['month'] = [tm.month for tm in df['time']]
    df = df.sort_values(['init_year', 'year', 'month']) #, 'season_year'])
    cols = [
        'ID', 'year', 'month', # TODO change to ID
        'project', 'source_id', 'mip', 'experiment',
        'member', 'init_year', 'variable', 'value'
    ]
    df = df[cols]
    return df


# def _spatial_xarray_to_dataframe_new(x, varname, metadata):
#     df = x.to_dataframe()
#     df = df.reset_index()
#     keep_cols = ['ID', 'season_year', varname]
#     drop_cols = [nm for nm in df if nm not in keep_cols]
#     # keep_cols = ['time', 'ID', 'season_year', varname]
#     # drop_cols = [nm for nm in df if nm not in keep_cols]
#     df = df.drop(drop_cols, axis=1)
#     df = df.set_index(keep_cols)
#     df = df.reset_index()
#     df = df.rename(
#         {varname: 'value'}, axis="columns"
#     ).reset_index()
#     df['project'] = metadata['project']
#     df['source_id'] = metadata['model']
#     df['mip'] = metadata['mip']
#     df['experiment'] = metadata['experiment']
#     df['member'] = metadata['ensemble']
#     df['init_year'] = metadata['init_year']
#     df['variable'] = varname
#     df = df.sort_values(['init_year', 'season_year'])
#     cols = [
#         'ID',                    # TODO change to ID
#         'project', 'source_id', 'mip', 'experiment',
#         'member', 'init_year', 'season_year', 'variable', 'value'
#     ]
#     df = df[cols]
#     return df


def _spatial_xarray_to_dataframe(x, varname, metadata):
    df = x.to_dataframe()
    df = df[[varname]]
    index_names = df.index.names
    drop_index = [nm for nm in index_names if nm not in ['season_year'] + VALID_X_NAMES + VALID_Y_NAMES]
    keep_index = [nm for nm in index_names if nm not in drop_index]
    df = df.reset_index()
    df = df.drop(drop_index, axis=1)
    df = df.set_index(keep_index)
    # df = xr.to_dataframe()[[xr.name]]
    df = df.reset_index()
    # df = df.rename({xr.name : 'var'}, axis=1)
    # Check time index
    # time_idx = [col for col in df if col in VALID_TIME_NAMES]
    # df = df.rename({time_idx[0] : 'time'}, axis=1)
    # Do we have spatial coordinates?
    lon_idx = [col for col in df if col in VALID_X_NAMES]
    lat_idx = [col for col in df if col in VALID_Y_NAMES]
    if (len(lon_idx) == 1) & (len(lat_idx) == 1):
        df = df.rename({lon_idx[0] : 'lon'}, axis=1)
        df = df.rename({lat_idx[0] : 'lat'}, axis=1)
        lats = ['N' + str(abs(x)) if x > 0 else 'S' + str(abs(x)) for x in df['lat']]
        lons = ['E' + str(abs(x)) if x > 0 else 'W' + str(abs(x)) for x in df['lon']]
        coords = [varname + '_' + lats[i] + '_' + lons[i] for i in range(len(lats))]
        df['variable'] = coords
        df = df.drop(['lat', 'lon'], axis=1)

    df = df.rename(
        {varname: 'value'}, axis="columns"
    ).reset_index()
    df['project'] = metadata['project']
    df['source_id'] = metadata['model']
    df['mip'] = metadata['mip']
    df['experiment'] = metadata['experiment']
    df['member'] = metadata['ensemble']
    df['init_year'] = metadata['init_year']
    # df['variable'] = varname

    df = df.sort_values(['init_year', 'season_year'])
    cols = [
        'project', 'source_id', 'mip', 'experiment',
        'member', 'init_year', 'season_year', 'variable', 'value'
    ]
    df = df[cols]
    return df


def _ensemble_field_preprocessor(inputdir, outputdir, station_ids, lat_coords, lon_coords, config):
    with open(config, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # Collect ESMValTool output files
    rootdir = inputdir
    recipe_output_dirs = \
        config['ensemble_data']['cmip5']['subdirectory'] \
        + config['ensemble_data']['cmip6']['subdirectory']

    fs = []
    for recipe in recipe_output_dirs:
        datadir = os.path.join(rootdir, recipe)
        prec_field_datadir = os.path.join(
            datadir,
            'work/precip_field/precip_field'
        )
        temp_field_datadir = os.path.join(
            datadir,
            'work/temp_field/temp_field'
        )
        recipe_fs = \
            glob.glob(prec_field_datadir + "/*.nc") \
            + glob.glob(temp_field_datadir + "/*.nc")
        fs += recipe_fs

    # TODO do not select seasons from data in EsmValTool
    # This will allow us to carry out the analysis across all seasons
    lats = xarray.DataArray(lat_coords, dims="ID")
    lons = xarray.DataArray(lon_coords, dims="ID")
    for i in tqdm(range(len(fs))):
        # f = '/Users/simonmoulds/projects/decadal-flood-prediction/data-raw/esmvaltool_output/recipe_s20_cmip6_autogen_20221214_120531/work/precip_field/precip_field/CMIP6_NorCPM1_Amon_dcppA-hindcast_s2010-r9i2p1f1_pr_gn_2011-2019.nc'
        f = fs[i]
        metadata = _parse_filepath(f)
        ds = xarray.open_dataset(f)
        lon_coord = _get_xarray_longitude_name(ds)
        ds.coords[lon_coord] = (ds.coords[lon_coord] + 180) % 360 - 180
        ds = ds.sortby(ds[lon_coord])
        varname = _get_variable_name(ds)
        # OLD:
        # x = _compute_mean_season(
        #     ds, varname, metadata['init_year'],
        #     season=[12, 1, 2, 3]
        # )
        # x = x.sel(lat=lats, lon=lons, method='nearest')
        x = ds.sel(lat=lats, lon=lons, method='nearest')
        x = x.assign_coords(ID=station_ids)
        df = _spatial_xarray_to_dataframe_new(x, varname, metadata)
        tbl = pa.Table.from_pandas(df, preserve_index=False)
        pq.write_to_dataset(
            tbl,
            root_path = os.path.join(outputdir, 'ensemble-forecast-field'),
            partition_cols = ['source_id', 'member', 'init_year', 'variable']
        )
        # OLD:
        # # Antecedent conditions [only for prec and temp]
        # if varname in ['precip_field', 'temp_field']:
        #     x = _compute_mean_season(
        #         ds, varname, metadata['init_year'], season=[9, 10, 11]
        #     )
        #     x = x.sel(lat=lats, lon=lons, method='nearest')
        #     x = x.assign_coords(ID=station_ids)
        #     df = _spatial_xarray_to_dataframe_new(x, varname, metadata)
        #     tbl = pa.Table.from_pandas(df, preserve_index=False)
        #     pq.write_to_dataset(
        #         tbl,
        #         root_path = os.path.join(outputdir, 'ensemble-forecast-field'),
        #         partition_cols = ['source_id', 'member', 'init_year', 'variable']
        #     )


def _ensemble_preprocessor(inputdir, outputdir, config):
    # Preprocessor for ensemble dataset which is
    # created by the two ESMValTOol jobs specified
    # in `02_run-esmvaltool-job.sh`
    with open(config, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # Collect ESMValTool output files
    # output_dir = config['output_data']['root']
    # output_dir = "data"
    # rootdir = config['input_data_root']
    rootdir = inputdir
    recipe_output_dirs = \
        config['ensemble_data']['cmip5']['subdirectory'] \
        + config['ensemble_data']['cmip6']['subdirectory']

    fs = []
    for recipe in recipe_output_dirs:
        datadir = os.path.join(rootdir, recipe)
        nao_datadir = os.path.join(datadir, 'work/nao/nao')
        ea_datadir = os.path.join(datadir, 'work/ea/ea')
        amv_datadir = os.path.join(datadir, 'work/amv/amv')
        prec_datadir = os.path.join(
            datadir,
            'work/european_precip/european_precip'
        )
        uk_prec_datadir = os.path.join(
            datadir,
            'work/uk_precip/uk_precip'
        )
        uk_temp_datadir = os.path.join(
            datadir,
            'work/uk_temp/uk_temp'
        )
        recipe_fs = glob.glob(nao_datadir + "/*.nc") \
            + glob.glob(ea_datadir + "/*.nc") \
            + glob.glob(amv_datadir + "/*.nc") \
            + glob.glob(prec_datadir + "/*.nc") \
            + glob.glob(uk_prec_datadir + "/*.nc") \
            + glob.glob(uk_temp_datadir + "/*.nc")
        fs += recipe_fs

    # Process files, sorting them on the basis of the filename
    for i in tqdm(range(len(fs))):
        # f = '/Users/simonmoulds/projects/decadal-flood-prediction/data-raw/esmvaltool_output/recipe_s20_cmip6_autogen_20221214_120531/work/uk_precip/uk_precip/CMIP6_NorCPM1_Amon_dcppA-hindcast_s2010-r9i2p1f1_pr_gn_2011-2019.nc'
        f = fs[i]
        metadata = _parse_filepath(f)
        ds = xarray.open_dataset(f)
        varname = _get_variable_name(ds)
        # OLD:
        # x = _compute_mean_season(ds, varname, metadata['init_year'], season=[12, 1, 2, 3])
        # df = _index_xarray_to_dataframe(x, varname, metadata)
        df = _index_xarray_to_dataframe_new(ds, varname, metadata)
        tbl = pa.Table.from_pandas(df, preserve_index=False)
        pq.write_to_dataset(
            tbl,
            root_path=os.path.join(outputdir, 'ensemble-forecast'),
            partition_cols=['source_id', 'member', 'init_year', 'variable']
        )
        # OLD:
        # # Antecedent conditions [only for prec and temp]
        # if varname in ['european_precip', 'uk_precip', 'uk_temp']:
        #     try:
        #         x = _compute_mean_season(
        #             ds, varname, metadata['init_year'], season=[9, 10, 11]
        #         )
        #     except ValueError:
        #         continue

        #     df = _index_xarray_to_dataframe(x, varname, metadata)
        #     df['variable'] = [var + '_antecedent' for var in df['variable']]
        #     tbl = pa.Table.from_pandas(df, preserve_index=False)
        #     pq.write_to_dataset(
        #         tbl,
        #         root_path=os.path.join(outputdir, 'ensemble-forecast'),
        #         partition_cols=['source_id', 'member', 'init_year', 'variable']
        #     )

    # # Now handle variables with a spatial dimension
    # fs = []
    # for recipe in recipe_output_dirs:
    #     datadir = os.path.join(rootdir, recipe)
    #     uk_prec_field_datadir = os.path.join(
    #         datadir,
    #         'work/uk_precip_field/uk_precip_field'
    #     )
    #     recipe_fs = glob.glob(uk_prec_field_datadir + "/*.nc")
    #     fs += recipe_fs

    # for i in tqdm(range(len(fs))):
    #     # f = '/Users/simonmoulds/projects/decadal-flood-prediction/data-raw/esmvaltool_output/recipe_s20_cmip6_autogen_20221116_114454/work/uk_precip_field/uk_precip_field/CMIP6_NorCPM1_Amon_dcppA-hindcast_s2010-r9i2p1f1_pr_gn_2011-2019.nc'
    #     f = fs[i]
    #     metadata = _parse_filepath(f)
    #     ds = xarray.open_dataset(f)
    #     varname = _get_variable_name(ds)
    #     x = _compute_mean_season(ds, varname, metadata['init_year'])
    #     df = _spatial_xarray_to_dataframe(x, varname, metadata)
    #     tbl = pa.Table.from_pandas(df, preserve_index=False)
    #     pq.write_to_dataset(
    #         tbl,
    #         root_path = os.path.join(outputdir, 'ensemble-forecast'),
    #         partition_cols = ['source_id', 'member', 'init_year', 'variable']
    #     )
    #     # Antecedent conditions [only for prec and temp]
    #     if varname in ['uk_precip_field']:
    #         x = _compute_mean_season(
    #             ds, varname, metadata['init_year'], season=[9, 10, 11]
    #         )
    #         df = _spatial_xarray_to_dataframe(x, varname, metadata)
    #         df['variable'] = [var + '_antecedent' for var in df['variable']]
    #         tbl = pa.Table.from_pandas(df, preserve_index=False)
    #         pq.write_to_dataset(
    #             tbl,
    #             root_path = os.path.join(outputdir, 'ensemble-forecast'),
    #             partition_cols = ['source_id', 'member', 'init_year', 'variable']
    #         )


def _coord_names(x):
    # Coordinate names from an Iris cube
    return tuple([coord.name() for coord in x.coords()])


def _get_time_name(x):
    coord_nms = _coord_names(x)
    return [nm for nm in coord_nms if nm in VALID_TIME_NAMES][0]


def _get_latitude_name(x):
    coord_nms = _coord_names(x)
    return [nm for nm in coord_nms if nm in VALID_Y_NAMES][0]


def _get_longitude_name(x):
    coord_nms = _coord_names(x)
    return [nm for nm in coord_nms if nm in VALID_X_NAMES][0]


def _get_xarray_longitude_name(x):
    coord_nms = list(x.coords)
    return [nm for nm in coord_nms if nm in VALID_X_NAMES][0]


def _matlab_mod(a, b):
    # Analog of Matlab's mod() function
    return a - b * np.floor(a / b)


def _has_positive_longitude(x):
    lon_name = _get_longitude_name(x)
    lon = x.coord(lon_name).points
    positive_lon = np.all(lon >= 0)
    return positive_lon


def _transform_longitude(xmin, xmax, positive_lon):
    # Transform standard (-180 to 180) longitude
    # value to (0 to 360) value, as is commonly
    # (always?) used by CMIP models
    if positive_lon:
        xmin = _matlab_mod(xmin, 360.)
        xmax = _matlab_mod(xmax, 360.)
    if xmin > xmax:
        xmin -= 360
    return xmin, xmax

def _regrid_cube(source, target):
    # Ensure grid coords have the same name
    target_lon_name = _get_longitude_name(target)
    target_lat_name = _get_latitude_name(target)
    source_lon_name = _get_longitude_name(source)
    source_lat_name = _get_latitude_name(source)
    source.coord(source_lon_name).rename(target_lon_name)
    source.coord(source_lat_name).rename(source_lat_name)

    # Make sure coord_systems are the same (this feels a bit hacky...)
    for coord_nm in [target_lat_name, target_lon_name]:
        source.coord(coord_nm).coord_system = target.coord(coord_nm).coord_system

    # source.coord(source_lat_name).attributes.pop('actual_range')
    # target.coord(target_lat_name).attributes.pop('actual_range')
    # source.coord(source_lon_name).attributes.pop('actual_range')
    # target.coord(target_lon_name).attributes.pop('actual_range')

    # Perform the regridding
    regrid_source = source.regrid(target, iris.analysis.Linear())
    return regrid_source


def _extract_subgrid(x, xmin, xmax, ymin, ymax):
    x_copy = x.copy()
    lon_name = _get_longitude_name(x_copy)
    lat_name = _get_latitude_name(x_copy)
    if (x_copy.coord(lon_name).bounds is None):
        x_copy.coord(lon_name).guess_bounds()
    if (x_copy.coord(lat_name).bounds is None):
        x_copy.coord(lat_name).guess_bounds()
    # positive_lon = _has_positive_longitude(x_copy)
    # xmin, xmax = _transform_longitude(xmin, xmax, positive_lon)
    box = x_copy.intersection(
        longitude=(xmin, xmax),
        latitude=(ymin, ymax)
    )
    return box

# def _extract_field(x, xmin, xmax, ymin, ymax):
#     x_copy = x.copy()
#     subset = x_copy.intersection(longitude=(xmin, xmax), latitude=(ymin, ymax))
#     return subset

def _extract_ts(x, xmin, xmax, ymin, ymax):
    # # work on a copy
    # x_copy = x.copy()
    # lon_name = _get_longitude_name(x_copy)
    # lat_name = _get_latitude_name(x_copy)
    # if (x_copy.coord(lon_name).bounds is None):
    #     x_copy.coord(lon_name).guess_bounds()
    # if (x_copy.coord(lat_name).bounds is None):
    #     x_copy.coord(lat_name).guess_bounds()
    # positive_lon = _has_positive_longitude(x_copy)
    # xmin, xmax = _transform_longitude(xmin, xmax, positive_lon)
    # box = x_copy.intersection(
    #     longitude=(xmin, xmax),
    #     latitude=(ymin, ymax)
    # )
    box = _extract_subgrid(x, xmin, xmax, ymin, ymax)
    grid_areas = iris.analysis.cartography.area_weights(box)
    ts = box.collapsed(
        ['latitude','longitude'],
        iris.analysis.MEAN,
        weights=grid_areas
    )
    return ts

# # NINO indices require SST
# def _extract_nino1_ts(x):
#     nino1 = _extract_ts(x, -90, -80, -10, -5)
#     return nino1

# def _extract_nino2_ts(x):
#     nino2 = _extract_ts(x, -90, -80, -5, 0)
#     return nino2

# def _extract_nino12_ts(x):
#     nino12 = _extract_ts(x, -90, -80, -10, 0)
#     return nino12

# def _extract_nino3_ts(x):
#     nino3 = _extract_ts(x, -150, -90, -5, 5)
#     return nino3

# def _extract_nino34_ts(x):
#     nino34 = _extract_ts(x, -170, -120, -5, 5)
#     return nino34

# def _extract_nino4_ts(x):
#     nino4 = _extract_ts(x, 160, -150, -5, 5)
#     return nino4

# def _extract_iod_ts(x):
#     iod_west = _extract_ts(x, 50, 70, -10, 10)
#     iod_east = _extract_ts(x, 90, 110, -10, 0)
#     return iod_west - iod_east

# def _extract_pdv_ts(x):
#     tropical = _extract_ts(x, -160, -110, -10, 6)
#     northern = _extract_ts(x, -180, -145, 30, 45)
#     return tropical - northern

def _extract_nao_ts(x):
    iceland_ts = _extract_ts(x, -25, -16, 63, 70)
    azores_ts = _extract_ts(x, -28, -20, 36, 40)
    nao = azores_ts - iceland_ts
    return nao

def _extract_ea_ts(x):
    ycoord = 52.5
    xcoord = -27.5
    positive_lon = _has_positive_longitude(x)
    if positive_lon:
        xcoord = _matlab_mod(xcoord, 360.)
    xr = xarray.DataArray.from_iris(x)
    ea = xr.sel(lat = ycoord, lon = xcoord, method = 'nearest')
    ea = ea.to_iris()
    return ea

def _extract_amv_ts(x):
    atlantic_ts = _extract_ts(x, -80, 0, 0, 60)
    globe_ts = _extract_ts(x, -180, 180, -60, 60)
    amv = atlantic_ts - globe_ts
    return amv

def _extract_atlantic_ts(x):
    atlantic_ts = _extract_ts(x, -80, 0, 0, 60)
    return atlantic_ts

def _extract_globe_ts(x):
    globe_ts = _extract_ts(x, -180, 180, -60, 60)
    return globe_ts

def _extract_european_precip_ts(x):
    europe_prec = _extract_ts(x, -10, 25, 55, 70)
    return europe_prec

def _extract_uk_precip_ts(x):
    uk_prec = _extract_ts(x, -8, 2, 50, 59)
    return uk_prec

def _extract_uk_precip_field(x):
    uk_precip_field = _extract_subgrid(x, -8, 2, 50, 59)
    return uk_precip_field

def _extract_uk_temp_ts(x):
    uk_temp = _extract_ts(x, -8, 2, 50, 59)
    return uk_temp

def _convert_ncdc(ncdc_filename):
    # First of all we find the rows with the date
    with open(ncdc_filename, 'r') as f:
        daterows = []
        daterows_index = []
        for i, v in enumerate(f):
            if re.match('^[0-9]+\\s+[0-9]{4}', v.strip()):
                daterows_index.append(i)
                daterows.append(v.strip().replace('\n', ''))

    # Compute dimension variables
    lat_values = np.arange(-87.5, 90, 5)
    lon_values = np.arange(-177.5, 180, 5)
    time_values = []
    time_bnds_values = []
    time_units = 'days since 1800-01-01'
    time_calendar = 'standard'
    for mon_year in daterows:
        month, year = mon_year.split()
        month, year = int(month), int(year)
        # Mid-point
        tm = netCDF4.date2num(
            datetime.datetime(year, month, 15, 0, 0),
            time_units, time_calendar
        )
        # Lower bound
        tm0 = netCDF4.date2num(
            datetime.datetime(year, month, 1, 0, 0),
            time_units, time_calendar
        )
        # Upper bound (account for end of year)
        next_month = month + 1
        next_year = year
        if month == 12:
            next_month = 1
            next_year = year + 1
        tm1 = netCDF4.date2num(
            datetime.datetime(next_year, next_month, 1, 0, 0),
            time_units, time_calendar
        )
        time_values.append(tm)
        time_bnds_values.append((tm0, tm1))

    # Set up netCDF file
    ds = netCDF4.Dataset('/tmp/ncdc-merged-sfc-mntp.nc', 'w')
    ds.createDimension('time', None)
    ds.createDimension('lat', len(lat_values))
    ds.createDimension('lon', len(lon_values))
    ds.createDimension('nv', 2)

    var = ds.createVariable('time', 'i4', ('time',))
    var.units = time_units
    var.calendar = time_calendar
    var[:] = np.array(time_values)

    var = ds.createVariable('time_bnds', 'i4', ('time', 'nv'))
    var.units = time_units
    var.calendar = time_calendar
    var[:] = np.array(time_bnds_values)

    var = ds.createVariable('lat', np.float32, ('lat',))
    var.units = 'degrees_north'
    var[:] = lat_values

    var = ds.createVariable('lon', np.float32, ('lon',))
    var.units = 'degrees_east'
    var[:] = lon_values

    # F4_FILLVAL = netCDF4.default_fillvals['f4']
    var = ds.createVariable('tempanomaly', np.float32, ('time', 'lat', 'lon'), fill_value=np.nan)
    var.units = 'K'

    with open(ncdc_filename, 'r') as f:

        time_index = -1
        lat_index = 0
        for i, v in enumerate(f):
            if i in daterows_index:
                time_index += 1
                lat_index = 0
            else:
                vals = np.array([float(x) for x in v.strip().split()])
                vals[vals <= -9998.] = np.nan #F4_FILLVAL
                var[time_index, lat_index, :] = vals / 100
                lat_index += 1
    ds.close()

def _extract_index(ds, index_name, column_name=None):
    # if index_name == "enso1":
    #     ds_index = _extract_enso1_ts(ds)
    # elif index_name == "enso2":
    #     ds_index = _extract_enso2_ts(ds)
    # elif index_name == "enso12":
    #     ds_index = _extract_enso12_ts(ds)
    # elif index_name == "enso3":
    #     ds_index = _extract_enso3_ts(ds)
    # elif index_name == "enso34":
    #     ds_index = _extract_enso34_ts(ds)
    # elif index_name == "enso4":
    #     ds_index = _extract_enso4_ts(ds)
    # elif index_name == "iod":
    #     ds_index = _extract_iod_ts(ds)
    # elif index_name == "pdv":
    #     ds_index = _extract_pdv_ts(ds)
    # elif index_name == "nao":
    if index_name == "nao":
        ds_index = _extract_nao_ts(ds)
    elif index_name == "ea":
        ds_index = _extract_ea_ts(ds)
    elif index_name == "amv":
        ds_index = _extract_amv_ts(ds)
    elif index_name == "european_precip":
        ds_index = _extract_european_precip_ts(ds)
    elif index_name == "uk_precip":
        ds_index = _extract_uk_precip_ts(ds)
    elif index_name == "uk_precip_field":
        ds_index = _extract_uk_precip_field(ds)
    elif index_name == "uk_temp":
        ds_index = _extract_uk_temp_ts(ds)
    elif index_name == "atlantic":
        ds_index = _extract_atlantic_ts(ds)
    elif index_name == "globe":
        ds_index = _extract_globe_ts(ds)
    else:
        raise "Index not yet implemented"

    if column_name is None:
        column_name = index_name

    # Alternative to below
    xr = xarray.DataArray.from_iris(ds_index)
    xr.name = column_name
    df = xr.to_dataframe()
    time_idx = [col for col in df if col in VALID_TIME_NAMES]
    if len(time_idx) > 0:
        df = df.set_index(time_idx[0], append=True)
    df = df[[column_name]]
    index_names = df.index.names
    drop_index = [nm for nm in index_names if nm not in VALID_TIME_NAMES + VALID_X_NAMES + VALID_Y_NAMES]
    keep_index = [nm for nm in index_names if nm not in drop_index]
    df = df.reset_index()
    df = df.drop(drop_index, axis=1)
    df = df.set_index(keep_index)

    # df = xr.to_dataframe()[[xr.name]]
    df = df.reset_index()
    # df = df.rename({xr.name : 'var'}, axis=1)
    # Check time index
    time_idx = [col for col in df if col in VALID_TIME_NAMES]
    df = df.rename({time_idx[0]: 'time'}, axis=1)
    # Do we have spatial coordinates?
    lon_idx = [col for col in df if col in VALID_X_NAMES]
    lat_idx = [col for col in df if col in VALID_Y_NAMES]
    if (len(lon_idx) == 1) & (len(lat_idx) == 1):
        df = df.rename({lon_idx[0] : 'lon'}, axis=1)
        df = df.rename({lat_idx[0] : 'lat'}, axis=1)
        lats = ['N' + str(abs(x)) if x > 0 else 'S' + str(abs(x)) for x in df['lat']]
        lons = ['E' + str(abs(x)) if x > 0 else 'W' + str(abs(x)) for x in df['lon']]
        coords = [xr.name + '_' + lats[i] + '_' + lons[i] for i in range(len(lats))]
        df['coord'] = coords
        df = df.pivot(index='time', columns=['coord'], values=xr.name)
        df = df.reset_index()

    tm = df['time'].dt.to_pydatetime()
    df['year'] = [tm.year for tm in tm]
    df['month'] = [tm.month for tm in tm]
    df = df.drop('time', axis=1)
    # Check for any duplicated values which we have to remove
    _, idx = np.unique(tm, return_index = True)
    df = df.iloc[idx]
    df.reset_index()

    # # Convert to data frame
    # df = iris.pandas.as_data_frame(ds_index)
    # df = df.rename(columns={0 : column_name})
    # df.index.name = 'time'
    # df = df.reset_index(level='time')
    # tm = df['time'].values
    # df['year'] = [tm.year for tm in tm]
    # df['month'] = [tm.month for tm in tm]
    # df = df.drop('time', axis=1)
    # # Check for any duplicated values which we have to remove
    # _, idx = np.unique(tm, return_index = True)
    # df = df.iloc[idx]
    # df.reset_index()

    return df


def _parse_observed_inputs(config, inputdir):
    # Must be a way to avoid repeating this
    with open(config, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    observed_root = inputdir
    # TESTING
    # config = {
    #     'observed_data': {
    #         'hadslp2r': {'subdirectory': 'HadSLP2r'},
    #         'gpcc': {'subdirectory': 'GPCC'},
    #         'hadcrut4': {'subdirectory': 'HadCRUT4'},
    #         'hadisst': {'subdirectory': 'HadISST'},
    #         'giss': {'subdirectory': 'GISS'},
    #         'ncdc': {'subdirectory': 'NCDC'}}}
    # observed_root = '/Users/simonmoulds/projects/decadal-flood-prediction/data-raw/observed_data'
    # END TESTING
    hadslp2r_filename = os.path.join(
        observed_root,
        config['observed_data']['hadslp2r']['subdirectory'],
        'slp.mnmean.real.nc'
    )
    gpcc_filename = os.path.join(
        observed_root,
        config['observed_data']['gpcc']['subdirectory'],
        'precip.mon.total.2.5x2.5.v2018.nc'
    )
    hadcrut4_filename = os.path.join(
        observed_root,
        config['observed_data']['hadcrut4']['subdirectory'],
        'HadCRUT.4.6.0.0.median.nc'
    )
    hadsst_filename = os.path.join(
        observed_root,
        config['observed_data']['hadisst']['subdirectory'],
        'HadISST_sst.nc'
    )
    giss_filename = os.path.join(
        observed_root,
        config['observed_data']['giss']['subdirectory'],
        'gistemp1200_GHCNv4_ERSSTv5.nc'
    )
    ncdc_filename = os.path.join(
        observed_root,
        config['observed_data']['ncdc']['subdirectory'],
        'ncdc-merged-sfc-mntp.dat'
    )
    out = {
        'HadSLP2r': hadslp2r_filename,
        'GPCC': gpcc_filename,
        'HadCRUT4': hadcrut4_filename,
        'HadSST': hadsst_filename,
        'GISS': giss_filename,
        'NCDC': ncdc_filename
    }
    return out


def _observed_field_preprocessor(inputdir, outputdir, station_ids, lat_coords, lon_coords, config):

    input_filenames = _parse_observed_inputs(config, inputdir)
    hadslp2r_filename = input_filenames['HadSLP2r']
    gpcc_filename = input_filenames['GPCC']
    hadcrut4_filename = input_filenames['HadCRUT4']
    hadsst_filename = input_filenames['HadSST']
    giss_filename = input_filenames['GISS']
    ncdc_filename = input_filenames['NCDC']

    def myfun(x, varname):
        df = x.to_dataframe()
        df = df.reset_index()
        keep_cols = ['time', 'ID', varname]
        drop_cols = [nm for nm in df if nm not in keep_cols]
        df = df.drop(drop_cols, axis=1)
        df = df.set_index(['time', 'ID'])
        df = df.reset_index()
        df = df.rename(
            {varname: 'value'}, axis="columns"
        )
        df['year'] = [tm.year for tm in df['time']]
        df['month'] = [tm.month for tm in df['time']]
        df = df[['ID', 'year', 'month', 'value']]
        return df

    # ################# #
    # UK rainfall field
    # ################# #
    # import iris
    # gpcc_filename = '/Users/simonmoulds/projects/decadal-flood-prediction/data-raw/observed_data/GPCC/precip.mon.total.2.5x2.5.v2018.nc'
    gpcc_filename_interp = os.path.join(
        outputdir,
        os.path.basename(os.path.splitext(gpcc_filename)[0]) + '_interp.nc'
    )
    # FIXME - ideally use another dataset without missing data
    subprocess.run(['cdo', 'fillmiss', gpcc_filename, gpcc_filename_interp])

    source = iris.load_cube(gpcc_filename_interp, 'precip')
    # This is OK because the original data has 0-360 longitude as well
    target = regional_stock_cube({
        'start_longitude': 0,
        'end_longitude': 355,
        'step_longitude': 5,
        'start_latitude': 90,
        'end_latitude': -90,
        'step_latitude': -5
    })
    ds_regrid = _regrid_cube(source, target)
    ds = xarray.DataArray.from_iris(ds_regrid)

    # Ensure the dataset has longitudes from -180 to +180
    lon_coord = _get_xarray_longitude_name(ds)
    ds.coords[lon_coord] = (ds.coords[lon_coord] + 180) % 360 - 180
    ds = ds.sortby(ds[lon_coord])

    lats = xarray.DataArray(lat_coords, dims="ID")
    lons = xarray.DataArray(lon_coords, dims="ID")
    x = ds.sel(lat=lats, lon=lons, method='nearest')

    x = x.assign_coords(ID=station_ids)
    gpcc = myfun(x, 'precip')
    gpcc['days_in_month'] = gpcc.apply(
        lambda x: monthrange(int(x['year']), int(x['month']))[1],
        axis=1
    )
    gpcc['value'] /= gpcc['days_in_month']
    gpcc = gpcc[['ID', 'year', 'month', 'value']]
    prec_df = gpcc.rename(columns={'value': 'precip_field'})

    # ################# #
    # Temperature fields
    # ################# #

    # Average over three datasets [check: do we need to regrid?]

    # 1 - HadCRUT4
    ds = iris.load_cube(hadcrut4_filename, 'temperature_anomaly')
    ds = xarray.DataArray.from_iris(ds)
    lon_coord = _get_xarray_longitude_name(ds)
    ds.coords[lon_coord] = (ds.coords[lon_coord] + 180) % 360 - 180
    ds = ds.sortby(ds[lon_coord])
    x = ds.sel(latitude=lats, longitude=lons, method='nearest')
    x = x.assign_coords(ID=station_ids)
    hct4 = myfun(x, 'temperature_anomaly')
    hct4 = hct4.rename(columns={'value': 'HadCRUT4'})

    # 2 - GISS
    ds = iris.load_cube(giss_filename, 'tempanomaly')
    ds = xarray.DataArray.from_iris(ds)
    lon_coord = _get_xarray_longitude_name(ds)
    ds.coords[lon_coord] = (ds.coords[lon_coord] + 180) % 360 - 180
    ds = ds.sortby(ds[lon_coord])
    x = ds.sel(lat=lats, lon=lons, method='nearest')
    x = x.assign_coords(ID=station_ids)
    giss = myfun(x, 'tempanomaly')
    giss = giss.rename(columns={'value': 'GISS'})

    # NCDC [needs to be converted to netCDF first of all]
    _convert_ncdc(ncdc_filename)
    ds = iris.load_cube('/tmp/ncdc-merged-sfc-mntp.nc', 'tempanomaly')
    ds = xarray.DataArray.from_iris(ds)
    lon_coord = _get_xarray_longitude_name(ds)
    ds.coords[lon_coord] = (ds.coords[lon_coord] + 180) % 360 - 180
    ds = ds.sortby(ds[lon_coord])
    x = ds.sel(lat=lats, lon=lons, method='nearest')
    x = x.assign_coords(ID=station_ids)
    ncdc = myfun(x, 'tempanomaly')
    ncdc = ncdc.rename(columns={'value': 'NCDC'})

    temp_df = reduce(lambda left, right: pd.merge(left, right), [hct4, giss, ncdc])
    temp_df['temp_field'] = temp_df[['HadCRUT4', 'GISS', 'NCDC']].mean(axis=1)
    temp_df = temp_df[['ID', 'year', 'month', 'temp_field']]

    # Merge precipitation and temperature
    combined_df = pd.merge(prec_df, temp_df)

    combined_df = combined_df.reset_index()
    combined_df = pd.melt(
        combined_df,
        id_vars=['ID', 'year', 'month'],
        value_vars=['temp_field', 'precip_field']
    )
    # Save dataset
    table = pa.Table.from_pandas(combined_df, preserve_index=False)
    # table = pa.Table.from_pandas(combined_df)
    # pq.write_table(table, os.path.join(outputdir, 'obs_field.parquet'))
    pq.write_to_dataset(
        table,
        root_path = os.path.join(outputdir, 'observed-field'),
        partition_cols = ['ID']
    )


def _observed_preprocessor(inputdir, outputdir, config):

    input_filenames = _parse_observed_inputs(config, inputdir)
    hadslp2r_filename = input_filenames['HadSLP2r']
    gpcc_filename = input_filenames['GPCC']
    hadcrut4_filename = input_filenames['HadCRUT4']
    hadsst_filename = input_filenames['HadSST']
    giss_filename = input_filenames['GISS']
    ncdc_filename = input_filenames['NCDC']
    with open(config, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # ####################### #
    # Sea-surface temperature
    # ####################### #

    # ds = iris.load_cube(hadsst_filename, 'sst')
    # enso1_df = _extract_index(ds, "enso1")
    # enso2_df = _extract_index(ds, "enso2")
    # enso3_df = _extract_index(ds, "enso3")
    # enso34_df = _extract_index(ds, "enso34")
    # enso4_df = _extract_index(ds, "enso4")
    # iod_df = _extract_index(ds, "iod")
    # pdv_df = _extract_index(ds, "pdv")

    # ####################### #
    # Mean sea-level pressure #
    # ####################### #
    ds = iris.load_cube(hadslp2r_filename, 'slp')
    nao_df = _extract_index(ds, "nao")
    ea_df = _extract_index(ds, "ea")

    # ################# #
    # European rainfall #
    # ################# #
    ds = iris.load_cube(gpcc_filename, 'precip')
    europe_precip_df = _extract_index(ds, "european_precip")
    europe_precip_df['days_in_month'] = europe_precip_df.apply(
        lambda x: monthrange(int(x['year']), int(x['month']))[1],
        axis=1
    )
    europe_precip_df['european_precip'] /= europe_precip_df['days_in_month']
    europe_precip_df = europe_precip_df.drop('days_in_month', axis=1)

    # ################# #
    # UK rainfall       #
    # ################# #
    ds = iris.load_cube(gpcc_filename, 'precip')
    uk_precip_df = _extract_index(ds, "uk_precip")
    uk_precip_df['days_in_month'] = uk_precip_df.apply(
        lambda x: monthrange(int(x['year']), int(x['month']))[1],
        axis=1
    )
    uk_precip_df['uk_precip'] /= uk_precip_df['days_in_month']
    uk_precip_df = uk_precip_df.drop('days_in_month', axis=1)

    # # ################# #
    # # UK rainfall field
    # # ################# #
    # # import iris
    # # gpcc_filename = '/Users/simonmoulds/projects/decadal-flood-prediction/data-raw/observed_data/GPCC/precip.mon.total.2.5x2.5.v2018.nc'
    # source = iris.load_cube(gpcc_filename, 'precip')
    # # target = iris.load_cube(hadslp2r_filename, 'slp')
    # # from esmvalcore import *
    # target = regional_stock_cube({
    #     'start_longitude': 0,
    #     'end_longitude': 355,
    #     'step_longitude': 5,
    #     'start_latitude': 90,
    #     'end_latitude': -90,
    #     'step_latitude': -5
    # })
    # ds_regrid = _regrid_cube(source, target)
    # # xr = xarray.DataArray.from_iris(ds_regrid)
    # # xr.to_netcdf("../../test.nc")
    # uk_precip_field_df = _extract_index(ds_regrid, 'uk_precip_field')
    # uk_precip_field_df['days_in_month'] = uk_precip_field_df.apply(
    #     lambda x: monthrange(int(x['year']), int(x['month']))[1],
    #     axis=1
    # )
    # cols = [nm for nm in uk_precip_field_df if nm.startswith('uk_precip_field')]
    # for col in cols:
    #     uk_precip_field_df[col] /= uk_precip_field_df['days_in_month']

    # # ################### #
    # # Gridded UK rainfall #
    # # ################### #
    # gpcc_filename = '/Users/simonmoulds/projects/decadal-flood-prediction/data-raw/GPCC/precip.mon.total.2.5x2.5.v2018.nc'
    # ds = iris.load_cube(gpcc_filename, 'precip')
    # uk_precip_grid = _extract_subgrid(ds, -8, 2, 50, 60)
    # uk_precip_arr = xarray.DataArray.from_iris(uk_precip_grid)
    # df = xarray.DataArray.to_dataframe(uk_precip_arr).reset_index()
    # df['lat'] = [str(lat) + 'N' if lat >= 0 else str(abs(lat)) + 'S' for lat in df['lat']]
    # df['lon'] = [str(lon) + 'E' if lon >= 0 else str(abs(lon)) + 'W' for lon in df['lon']]
    # df['latlon'] = df['lat'] + '_' + df['lon']
    # df = df.drop(['lat', 'lon'], axis = 1)
    # df = df.set_index(['time', 'latlon'])
    # df = df.pivot_table(index = 'time', columns = 'latlon', values = 'precip')
    # # EOF analysis
    # coslat = np.cos(np.deg2rad(uk_precip_arr.coords['lat'].values)).clip(0, 1)
    # wgts = np.sqrt(coslat)[..., np.newaxis]
    # solver = Eof(uk_precip_arr, weights = wgts)
    # variance_explained = list(np.cumsum(np.array(solver.varianceFraction())))

    # ################################# #
    # Atlantic Multidecadal Oscillation #
    # ################################# #

    # This is slightly more complicated because it is
    # the average value across three datasets

    # HADCRUT4_FILENAME = 'HadCRUT4/4.6.0.0.anomalies.{member:d}.nc'

    # 1 - HadCRUT4
    ds = iris.load_cube(hadcrut4_filename, 'temperature_anomaly')
    hct4_atlantic_df = _extract_index(ds, "atlantic", "hct4_atlantic_temp")
    hct4_globe_df = _extract_index(ds, "globe", "hct4_globe_temp")
    hct4_df = pd.merge(hct4_atlantic_df, hct4_globe_df)
    hct4_uk_df = _extract_index(ds, "uk_temp", "hct4_uk_temp")

    # 2 - GISS
    ds = iris.load_cube(giss_filename, 'tempanomaly')
    giss_atlantic_df = _extract_index(ds, "atlantic", "giss_atlantic_temp")
    giss_globe_df = _extract_index(ds, "globe", "giss_globe_temp")
    giss_df = pd.merge(giss_atlantic_df, giss_globe_df)
    giss_uk_df = _extract_index(ds, "uk_temp", "giss_uk_temp")

    # NCDC [needs to be converted to netCDF first of all]
    _convert_ncdc(ncdc_filename)
    ds = iris.load_cube('/tmp/ncdc-merged-sfc-mntp.nc', 'tempanomaly')
    ncdc_atlantic_df = _extract_index(ds, "atlantic", "ncdc_atlantic_temp")
    ncdc_globe_df = _extract_index(ds, "globe", "ncdc_globe_temp")
    ncdc_df = pd.merge(ncdc_atlantic_df, ncdc_globe_df)
    ncdc_uk_df = _extract_index(ds, "uk_temp", "ncdc_uk_temp")

    # Now merge datasets
    temp_df = reduce(
        lambda left, right: pd.merge(left, right),
        [hct4_df, giss_df, ncdc_df, hct4_uk_df, giss_uk_df, ncdc_uk_df]
    )
    # Take mean temperature for Atlantic and global regions
    temp_datasets = ['hct4', 'giss']#, 'ncdc']
    atlantic_cols = [ds + '_atlantic_temp' for ds in temp_datasets]
    globe_cols = [ds + '_globe_temp' for ds in temp_datasets]
    uk_cols = [ds + '_uk_temp' for ds in temp_datasets]

    temp_df['mean_atlantic_temp'] = temp_df[atlantic_cols].mean(axis=1)
    temp_df['mean_globe_temp'] = temp_df[globe_cols].mean(axis=1)
    temp_df['amv'] = temp_df['mean_atlantic_temp'] - temp_df['mean_globe_temp']
    temp_df['uk_temp'] = temp_df[uk_cols].mean(axis=1)

    # Now merge all time series
    obs_df = reduce(
        lambda left, right: pd.merge(left, right, how='outer'),
        [temp_df, nao_df, ea_df, europe_precip_df, uk_precip_df]
    )
    obs_df = obs_df[['year', 'month', 'nao', 'ea', 'amv', 'european_precip', 'uk_precip', 'uk_temp']]
    # Pivot to long format
    obs_df = obs_df.reset_index()
    obs_df = pd.melt(
        obs_df,
        id_vars=['year', 'month'],
        # value_vars=['nao', 'ea', 'amv','european_precip', 'uk_precip', 'uk_temp']
        value_vars=[nm for nm in obs_df if nm not in ['year', 'month', 'index']]
    )
    table = pa.Table.from_pandas(obs_df)
    # Save dataset
    pq.write_table(table, os.path.join(outputdir, 'obs.parquet'))
