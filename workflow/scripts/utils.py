#!/usr/bin/env python3

import os
import shutil
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


def _compute_mean_djfm(ds, varname, init_year, start_year=2, end_year=9):
    # Compute the mean value for DJFM season.
    #
    # Args:
    #   ds        : xarray dataset
    #   varname   : string. Variable name.
    #   init_year : int. Initialization year of `ds`
    #
    # Returns:
    #   xarray dataset.

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

    ds['counter'] = ((time_dimname,), [1, ] * len(ds.time))
    ds = ds.groupby('season_year').sum(time_dimname)
    # Limit dataset to complete [boreal winter] seasons
    complete_index = ds['counter'].values == 4
    ds = ds.sel(season_year=complete_index)
    # Divide by four to get mean value
    ds[varname] = ds[varname] / 4
    # Add lead time
    # [N.B. season_year (added automatically by Iris)
    # defines the year of final month in the season]
    ds['lead_time'] = (('season_year',), ds.season_year.values - int(init_year))
    def _is_yr2to9(year): #, init_year):
        return (year >= start_year) & (year <= end_year)
    ds = ds.sel(season_year=_is_yr2to9(ds['lead_time']))
    # ds = ds.drop(['counter','lead_time'])
    # ds = ds.mean('season_year')
    return ds[varname]#.values


def _get_variable_name(ds):
    # Get variable name from dataset.
    varnames = [k for k, _ in ds.data_vars.items()]
    if 'nao' in varnames:
        return 'nao'
    elif 'ea' in varnames:
        return 'ea'
    elif 'amv' in varnames:
        return 'amv'
    elif 'european_precip' in varnames:
        return 'european_precip'
    elif 'uk_precip' in varnames:
        return 'uk_precip'
    elif 'uk_temp' in varnames:
        return 'uk_temp'


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

    # Output directories
    # decadal_output_path = os.path.join(
    #     output_dir, 'ensemble-forecast-decadal'
    # )
    # try:
    #     shutil.rmtree(decadal_output_path)
    # except OSError:
    #     pass
    # annual_output_path = os.path.join(
    #     output_dir, 'ensemble-forecast'
    # )
    # try:
    #     shutil.rmtree(annual_output_path)
    # except OSError:
    #     pass

    # Process files, sorting them on the basis of the filename
    data = []
    for i in tqdm(range(len(fs))):
        f = fs[i]
        meta = _parse_filepath(f)
        ds = xarray.open_dataset(f)
        varname = _get_variable_name(ds)
        val = _compute_mean_djfm(ds, varname, meta['init_year'])
        # print(val.to_dataframe())
        # rowdata = {
        #     'project': meta['project'],
        #     'source_id': meta['model'],
        #     'mip': meta['mip'],
        #     'experiment': meta['experiment'],
        #     'member': meta['ensemble'],
        #     'init_year': meta['init_year'],
        #     'variable' : varname,
        #     'value': val
        # }
        # rowdata = pd.DataFrame(rowdata)
        # # rowdata = pd.DataFrame(rowdata, index=[i])
        rowdata = val.to_dataframe()
        rowdata = rowdata.rename(
            {varname: 'value'}, axis="columns"
        ).reset_index()
        rowdata['project'] = meta['project']
        rowdata['source_id'] = meta['model']
        rowdata['mip'] = meta['mip']
        rowdata['experiment'] = meta['experiment']
        rowdata['member'] = meta['ensemble']
        rowdata['init_year'] = meta['init_year']
        rowdata['variable'] = varname

        rowdata = rowdata.sort_values(['init_year', 'season_year'])
        cols = [
            'project', 'source_id', 'mip', 'experiment',
            'member', 'init_year', 'season_year', 'variable', 'value'
        ]
        rowdata = rowdata[cols]
        rowdata_tbl = pa.Table.from_pandas(
            rowdata,
            preserve_index=False
        )
        pq.write_to_dataset(
            rowdata_tbl,
            root_path = os.path.join(outputdir, 'ensemble-forecast'),
            partition_cols = ['source_id', 'member', 'init_year', 'variable']
        )
        # NOT USED
        # # Take mean
        # cols.remove('season_year')
        # group_cols = cols[:]
        # group_cols.remove('value')
        # rowdata_mean = rowdata.groupby(group_cols).mean().reset_index()
        # rowdata_mean = rowdata_mean[cols]
        # rowdata_mean_tbl = pa.Table.from_pandas(
        #     rowdata_mean,
        #     preserve_index=False
        # )
        # pq.write_to_dataset(
        #     rowdata_mean_tbl,
        #     root_path=decadal_output_path,
        #     partition_cols=['source_id', 'member', 'init_year', 'variable']
        # )

    #     data.append(rowdata)

    # df = pd.concat(data)
    # df = df.sort_values(
    #     ['project', 'source_id', 'init_year', 'member'],
    #     axis=0, ignore_index=True
    # )
    # table = pa.Table.from_pandas(df)
    # pq.write_table(table, os.path.join(output_dir, 'ensemble.parquet'))


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
    positive_lon = _has_positive_longitude(x_copy)
    xmin, xmax = _transform_longitude(xmin, xmax, positive_lon)
    box = x_copy.intersection(
        longitude=(xmin, xmax),
        latitude=(ymin, ymax)
    )
    return box

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

def _extract_index(ds, index_name, column_name = None):
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

    # Convert to data frame
    df = iris.pandas.as_data_frame(ds_index)
    df = df.rename(columns={0 : column_name})
    df.index.name = 'time'
    df = df.reset_index(level='time')
    tm = df['time'].values
    df['year'] = [tm.year for tm in tm]
    df['month'] = [tm.month for tm in tm]
    df = df.drop('time', axis=1)
    # Check for any duplicated values which we have to remove
    _, idx = np.unique(tm, return_index = True)
    df = df.iloc[idx]
    df.reset_index()
    return df

# def _observed_preprocessor(config, outputdir):
def _observed_preprocessor(inputdir, outputdir, config):
    with open(config, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # observed_root = config['input_data_root']
    observed_root = inputdir
    print(observed_root)
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
    # hadsst_filename = "data-raw/HadISST/HadISST_sst.nc"
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
    # output_dir = config['output_data']['root']

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
        lambda left, right: pd.merge(left, right),
        [temp_df, nao_df, ea_df, europe_precip_df, uk_precip_df]
    )
    obs_df = obs_df[['year','month','nao','ea','amv','european_precip','uk_precip','uk_temp']]
    # Pivot to long format
    obs_df = obs_df.reset_index()
    obs_df = pd.melt(
        obs_df,
        id_vars=['year', 'month'],
        value_vars=['nao', 'ea', 'amv','european_precip', 'uk_precip', 'uk_temp']
    )
    table = pa.Table.from_pandas(obs_df)
    # Save dataset
    pq.write_table(table, os.path.join(outputdir, 'obs.parquet'))
