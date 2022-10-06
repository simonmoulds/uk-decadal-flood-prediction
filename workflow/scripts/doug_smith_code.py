#!/usr/bin/env python3

import numpy as np

def corr_cross_validate(fcst,obs,leave_out_add=0):
   """
   Compute cross validated anomaly correlation
   leave_out_add = no of additional points to leave out
   """
   ntimes = len(obs)
   corr = np.zeros(ntimes)
   use_data = np.ones(ntimes)
   for i in range(ntimes):
       leave_out_start = np.max([0,i-leave_out_add])
       leave_out_end = np.min([ntimes-1,i+leave_out_add])
       use_data[leave_out_start:leave_out_end] = 0.0
       corr_this = pearsonr(fcst[(use_data > 0.5)],obs[(use_data > 0.5)])
       corr[i] = corr_this[0]
       use_data[:] = 1.0
   return corr

def get_ts(fields,index,do_detrend = False,do_norm = False):
   """
   Get timeseries
   """
   years = fields.coord('year').points
   try:
       realization_dim_no, = fields.coord_dims(fields.coord('realization'))
   except:
       realization_dim_no = -1
   if index == 'NAO':
       iceland_box = fields.intersection(longitude=(-25,-16),latitude=(63,70))
       grid_areas_iceland = iris.analysis.cartography.area_weights(iceland_box)
       iceland = iceland_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_iceland)
       azores_box = fields.intersection(longitude=(-28,-20),latitude=(36,40))
       grid_areas_azores = iris.analysis.cartography.area_weights(azores_box)
       azores = azores_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_azores)
       ts = azores.data - iceland.data
   elif index == 'AO':
       arctic_box = fields.intersection(longitude=(0,360),latitude=(60,90))
       grid_areas_arctic = iris.analysis.cartography.area_weights(arctic_box)
       arctic = arctic_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_arctic)
       midlat_box = fields.intersection(longitude=(0,360),latitude=(30,50))
       grid_areas_midlat = iris.analysis.cartography.area_weights(midlat_box)
       midlat = midlat_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_midlat)
       ts = midlat.data - arctic.data
   elif index == 'SAM':
       southern_box = fields.intersection(longitude=(0,360),latitude=(63,67))
       grid_areas_southern = iris.analysis.cartography.area_weights(southern_box)
       southern = southern_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_southern)
       midlat_box = fields.intersection(longitude=(0,360),latitude=(38,42))
       grid_areas_midlat = iris.analysis.cartography.area_weights(midlat_box)
       midlat = midlat_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_midlat)
       ts = midlat.data - southern.data
   elif index == 'AMV':
       atl_box = fields.intersection(longitude=(280,360),latitude=(0,60))
       grid_areas_atl = iris.analysis.cartography.area_weights(atl_box)
       atl = atl_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_atl)
       gbl_box = fields.intersection(longitude=(0,360),latitude=(-60,60))
       grid_areas_gbl = iris.analysis.cartography.area_weights(gbl_box)
       gbl = gbl_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_gbl)
       ts = atl.data - gbl.data
   elif index == 'PDV':
       tropical_box = fields.intersection(longitude=(200,250),latitude=(-10,6))
       grid_areas_tropical = iris.analysis.cartography.area_weights(tropical_box)
       tropical = tropical_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_tropical)
       northern_box = fields.intersection(longitude=(180,215),latitude=(30,45))
       grid_areas_northern = iris.analysis.cartography.area_weights(northern_box)
       northern = northern_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_northern)
       ts = tropical.data - northern.data
       #ts = northern.data
   elif index == 'ipo':
       northern_box = fields.intersection(longitude=(140,215),latitude=(25,45))
       grid_areas_northern = iris.analysis.cartography.area_weights(northern_box)
       northern = northern_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_northern)
       southern_box = fields.intersection(longitude=(150,200),latitude=(-50,-15))
       grid_areas_southern = iris.analysis.cartography.area_weights(southern_box)
       southern = southern_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_southern)
       middle_box = fields.intersection(longitude=(170,270),latitude=(-10,10))
       grid_areas_middle = iris.analysis.cartography.area_weights(middle_box)
       middle = middle_box.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas_middle)
       ts = middle.data - 0.5*(northern.data + southern.data)
   elif index == 'SahelP':
       sahel = fields.intersection(longitude=(-16,36),latitude=(10,20))
       grid_areas = iris.analysis.cartography.area_weights(sahel)
       av = sahel.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
       ts = av.data
   elif index == 'NEuropeP':
       neurope = fields.intersection(longitude=(-10,25),latitude=(55,70))
       grid_areas = iris.analysis.cartography.area_weights(neurope)
       av = neurope.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
       ts = av.data
   elif index == 'MDR_slp':
       neurope = fields.intersection(longitude=(-80,-20),latitude=(10,20))
       grid_areas = iris.analysis.cartography.area_weights(neurope)
       av = neurope.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
       ts = av.data
   elif index == 'SPG':
       neurope = fields.intersection(longitude=(-60,-10),latitude=(50,65))
       grid_areas = iris.analysis.cartography.area_weights(neurope)
       av = neurope.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
       ts = av.data
   elif index == 'nino3.4':
       neurope = fields.intersection(longitude=(-170,-120),latitude=(-5,5))
       #neurope = fields.intersection(longitude=(-200,-150),latitude=(-5,5)) # nino4
       grid_areas = iris.analysis.cartography.area_weights(neurope)
       av = neurope.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
       ts = av.data
   elif index == 'GlobalT':
       grid_areas = iris.analysis.cartography.area_weights(fields)
       av = fields.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
       ts = av.data
   if realization_dim_no >= 0:
       ens_num = fields.coord('realization').points
       if realization_dim_no == 0:
           for i,ens in enumerate(ens_num):
               ts_this = ts[i,:]
               # detrend
               if do_detrend:
                   ts[i,:] = signal.detrend(ts_this)
               # normalise
               if do_norm:
                   ts[i,:] = (ts_this - np.mean(ts_this)) / np.std(ts_this)
       elif realization_dim_no == 1:
           for i,ens in enumerate(ens_num):
               ts_this = ts[:,i]
               # detrend
               if do_detrend:
                   ts[:,i] = signal.detrend(ts_this)
               # normalise
               if do_norm:
                   ts[:,i] = (ts_this - np.mean(ts_this)) / np.std(ts_this)
       else:
           print('realization coord error')
           import pdb; pdb.set_trace()
   else:
       # detrend
       if do_detrend:
           ts = signal.detrend(ts)
       # normalise
       if do_norm:
           ts = (ts - np.mean(ts)) / np.std(ts)
   return ts, years

def get_all_fields_ens(datadir_obs,cmip6dir,cmip5dir,var,lead,season,do_detrend = False,do_norm = False):
   obs_fields = get_fields(datadir_obs,'obs_'+var+lead+season+'_fields.nc')
   init_fields_cmip6 = get_fields(cmip6dir,'init_*'+var+lead+season+'_fields_??.nc')
   init_fields_cmip5 = get_fields(cmip5dir,'init_*'+var+lead+season+'_fields_??.nc')
   init_fields_extra = get_fields(cmip6dir,'extra_members/init_*'+var+lead+season+'_fields_??.nc')
   # regrid
   init_fields_cmip5.coord('latitude').coord_system = obs_fields.coord('latitude').coord_system
   init_fields_cmip5.coord('longitude').coord_system = obs_fields.coord('longitude').coord_system
   init_fields_cmip5 = init_fields_cmip5.regrid(obs_fields[0], iris.analysis.Linear())
   init_fields_cmip6 = init_fields_cmip6.regrid(obs_fields[0], iris.analysis.Linear())
   init_fields_extra = init_fields_extra.regrid(obs_fields[0], iris.analysis.Linear())
   # combine cmip5 and cmip6
   # ...common years
   init_fields_cmip6 = init_fields_cmip6.extract(iris.Constraint(year=init_fields_cmip5.coord('year').points))
   init_fields_extra = init_fields_extra.extract(iris.Constraint(year=init_fields_cmip5.coord('year').points))
   init_fields_cmip5 = init_fields_cmip5.extract(iris.Constraint(year=init_fields_cmip6.coord('year').points))
   # ...concatenate
   init_fields_cmip5 = convert_units(init_fields_cmip5, var)
   init_fields_cmip6.rename(init_fields_cmip5.standard_name)
   init_fields_extra.rename(init_fields_cmip5.standard_name)
   init_fields_cmip5.data = init_fields_cmip5.data.astype(init_fields_cmip6.dtype)
   init_fields_cmip5.cell_methods = None
   cubelist = iris.cube.CubeList([init_fields_cmip6,init_fields_cmip5,init_fields_extra])
   cubelist[0].metadata = cubelist[1].metadata
   cubelist[1].coord('year').var_name=cubelist[0].coord('year').var_name
   cubelist[2].metadata = cubelist[1].metadata
   init_fields = cubelist.concatenate_cube()
   # common years
   obs_fields = obs_fields.extract(iris.Constraint(year=init_fields.coord('year').points))
   init_fields = init_fields.extract(iris.Constraint(year=obs_fields.coord('year').points))
   # detrend
   if do_detrend:
       init_fields, fits = detrend(init_fields)
       obs_fields, fits = detrend(obs_fields)
   # normalise
   if do_norm:
       obs_fields = normalise(obs_fields)
       init_fields = normalise(init_fields)
   return obs_fields,init_fields

def convert_units(cube, var):
   """
   Make units the same for all models
   """
   #import pdb; pdb.set_trace()
   if var == 'tas':
       if cube.units != 'K':
           print('*** check units:', cube.units)
           #import pdb; pdb.set_trace()
   elif var == 'precip':
       if cube.units == 'kg m-2 s-1':
           cube.units = 'mm day^-1'
           cube.data = cube.data * 86400.0
       elif cube.units != 'mm day^-1':
           print('*** check units', cube.units)
           #import pdb; pdb.set_trace()
   elif var == 'mslp':
       if cube.units == 'Pa':
           cube.units = 'hPa'
           #cube.data = cube.data * 0.01
       elif cube.units != 'hPa':
           print('*** check units', cube.units)
           import pdb; pdb.set_trace()
   return cube

def get_fields(dirname,filename):
   """
   Get fields
   """
   #import pdb; pdb.set_trace()
   if("obs" in filename):
       fields = iris.load_cube(dirname+filename)
       if not fields[0].coords('year'):
           iris.coord_categorisation.add_year(fields, 'time')
   else:
       #fields = iris.load_cube(dirname+filename,callback=realization_metadata)
       fields_cubelist = iris.load(dirname+filename,callback=realization_metadata)
       # Unify and merge into a cube
       # ...common years
       yr1 = np.max([c.coord('year').points[0] for c in fields_cubelist])
       yr2 = np.min([c.coord('year').points[-1] for c in fields_cubelist])
       for i,cube in enumerate(fields_cubelist):
           fields_cubelist[i] = cube.extract(iris.Constraint(year = lambda y: yr1 <= y <= yr2))
       #import pdb; pdb.set_trace()
       # ...other things
       years_coord = fields_cubelist[0].coord('year')
       time_coord = fields_cubelist[0].coord('time').copy()
       for i,cube in enumerate(fields_cubelist):
           cube.attributes = None
           try:
               cube.remove_coord('height')
           except:
               print('No height coordinate')
           try:
               cube.remove_coord('month') # because DPS4 has numbers
           except:
               print('No month coordinate')
           if cube.coord('year') != years_coord:
               print('Year error')
               import pdb; pdb.set_trace()
           cube.remove_coord('time')
           #cube.add_aux_coord(time_coord,1)
           #cube.coord('time').points = time_coord.points
           #iris.util.promote_aux_coord_to_dim_coord(cube,'time')
           #cube.remove_coord('year')
       #import pdb; pdb.set_trace()
       #check_cubes_concatenate(fields_cubelist[0],fields_cubelist[1])
       fields = fields_cubelist.concatenate_cube()
   try:
       fields.coord('latitude').guess_bounds()
       fields.coord('longitude').guess_bounds()
   except:
       dummy = 0
   return fields

def get_fields_nao_match(fields, obsdir, cmip6dir, cmip5dir, lead, do_smooth, do_sample):
    """Extract subset of fields closest to ensemble mean forecast NAO
    """
    # SM: These first few lines obviously retrieve data for the analysis
    obs_fields_nao, init_fields_nao = get_all_fields_ens(
        obsdir, cmip6dir, cmip5dir, 'mslp', lead, '_djfm'
    )
    # SM: Based on the subsequent code, we can conclude that obs_nao is 1D
    # list/array, such that len(obs_nao) == nyear, and init_nao is a 2D
    # array with columns corresponding to years and rows to ensemble members
    obs_nao, year_obs_nao = get_ts(obs_fields_nao, 'NAO')
    init_nao, year_init_nao = get_ts(init_fields_nao, 'NAO')
    init_nao_em = []
    year_init_nao_em = []

    # SM: Calculate the ensemble mean
    for year in year_init_nao:
        year_init_nao_em.append(year)
        init_nao_em.append(np.mean(init_nao[:,(year_init_nao == year)]))

    # SM: list -> numpy array
    year_init_nao_em = np.asarray(year_init_nao_em)
    init_nao_em = np.asarray(init_nao_em)
    # SM: Saving a copy of the original ensemble mean
    init_nao_em_orig = init_nao_em.copy()

    # Smoothing [i.e. creating the lagged ensemble - see Methods]
    if do_smooth:
        # smooth
        istart = 3
        for iyr in range(istart, len(year_init_nao)):
            # Take the mean of the current year, year-1, year-2, year-3
            init_nao_em[iyr] = np.mean(
                [init_nao_em_orig[iyr], init_nao_em_orig[iyr-1],
                 init_nao_em_orig[iyr-2], init_nao_em_orig[iyr-3]]
            )

        # Remove data before `istart`
        init_nao_em = init_nao_em[istart:]
        year_init_nao_em = year_init_nao_em[istart:]
        #init_nao = init_nao[:,istart:] ; year_init_nao = year_init_nao[istart:]
        obs_nao = obs_nao[istart:] ; year_obs_nao = year_obs_nao[istart:]

    # Standardization - this is not described in Smith et al. 2019,
    # but it means the variance of the two is now the same (i.e. 1)
    init_nao_em = (
        (init_nao_em - np.mean(init_nao_em))
        / np.std(init_nao_em)
    )
    obs_nao = (
        (obs_nao - np.mean(obs_nao))
        / np.std(obs_nao)
    )
    # The nao_corr describes the variance in the observed accounted
    # for by the forecast, which is assumed to equate to the signal
    # in the observed [?]
    nao_corr = corr_cross_validate(init_nao_em, obs_nao, leave_out_add=1)
    for iyr,year in enumerate(year_obs_nao):
        init_nao_em[iyr] = init_nao_em[iyr] * nao_corr[iyr]

    # Select members based on NAO forecast
    nens = init_nao.shape[0]
    for ens in range(0,nens):
        init_nao[ens,:] = (init_nao[ens,:] - np.mean(init_nao[ens,:])) / np.std(init_nao[ens,:])

    years = year_init_nao_em
    init_sample_cubelist = iris.cube.CubeList()
    init_all_cubelist = iris.cube.CubeList()
    if do_smooth:
        nens_tot = (istart+1)*nens
    else:
        nens_tot = nens

    if do_sample:
        nsample = 20
    else:
        nsample = nens_tot

    # SM: This section seems to loop through each
    for iyr, year_this in enumerate(years):
        if do_smooth:
            init_this_cubelist = iris.cube.CubeList() ; error2 = []
            for iyr2,yr2 in enumerate([year_this, year_this-1, year_this-2, year_this-3]):

                fields_this = fields.extract(iris.Constraint(year=yr2))
                fields_this.coord('realization').points = range(iyr2*nens, (iyr2+1)*nens)
                fields_this.coord('year').points = year_this
                error2.append(abs(init_nao[:,(year_init_nao==yr2)] - init_nao_em[iyr]).data)
                init_this_cubelist.append(fields_this)

            fields_this = init_this_cubelist.concatenate_cube()
            error_this = np.asarray(error2).reshape(nens_tot)
        else:
            fields_this = fields.extract(iris.Constraint(year=year_this))
            error_this = abs(init_nao[:,iyr] - init_nao_em[iyr]).data

        error_this_sorted = error_this.copy() ; error_this_sorted.sort()
        if nsample < nens_tot:
            thresh = error_this_sorted[nsample]
        else:
            thresh = max(error_this_sorted)
        realization_this = fields_this.coord('realization').points
        realization_sample = realization_this[(error_this <= thresh)]
        realization_sample = realization_sample[:nsample]
        init_sample_this = fields_this.extract(iris.Constraint(realization=realization_sample))

        init_sample_this.coord('realization').points=range(nsample)
        init_sample_cubelist.append(init_sample_this)
        init_all_cubelist.append(fields_this)

    init_sample = init_sample_cubelist.merge_cube()
    init_all = init_all_cubelist.merge_cube()
    return init_sample, init_all
