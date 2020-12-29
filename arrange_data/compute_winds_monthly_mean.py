#compute monthly means of u and v
import numpy as np
import xarray as xr
import pandas as pd
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '/storage/shared/glusterfs/acrcc/users/ys916780/Data_for_Sol/output'
winds = []
for Y in np.arange(1981, 2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + '_' + str(Y)+ '_uv200', engine='cfgrib', chunks={'step': 10},
				     backend_kwargs={'indexpath': ''})
		ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		winds.append(ds)

winds = xr.concat(winds, dim='year')

winds.to_netcdf('../data/winds200.nc')

RUTA = '~/datos/UV_S4Hindcasts_0108_200hPalevels_24Hourly_'
winds = []
for Y in np.arange(1981, 1993):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib', chunks={'step': 10},
				     backend_kwargs={'indexpath': ''})
		ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
		ds = ds.sel(**{'IsobaricInhPa': [200]})
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		winds.append(ds)


RUTA = '~/datos/UV_S4Hindcasts_0108_levels_24Hourly_'
for Y in np.arange(1993, 2014):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib', chunks={'step': 10},
				     backend_kwargs={'indexpath': ''})
		ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
		ds = ds.sel(**{'isobaricInhPa': [200]}).squeeze()
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		winds.append(ds)

RUTA = '~/UV_S4Hindcasts_0108_levels_24Hourly_'

for Y in np.arange(2014, 2018):
	ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib', chunks={'step': 10},
			     backend_kwargs={'indexpath': ''})
	ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
	ds = ds.sel(**{'isobaricInhPa': [200]}).squeeze()
	ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
	ds = ds.rename({'step': 'time'}).set_coords(['time'])
	ds = ds.groupby('time.month').mean(dim='time')
	winds.append(ds)

winds =xr.concat(winds, dim='year')

winds.to_netcdf('~/datos/data/winds200_djf.nc')

RUTA = '~/datos/UV_S4Hindcasts_50hPa_24Hourly_'
winds = []
for Y in np.arange(1981,1993):
	ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib', chunks={'step': 10})
	ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
	ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
	ds = ds.rename({'step': 'time'}).set_coords(['time'])
	ds = ds.groupby('time.month').mean(dim='time', skipna='True')
	winds.append(ds)
winds = xr.concat(winds, dim='year')
winds.to_netcdf('~/datos/data/winds50_asondjf_81-92.nc4')

RUTA = '~/datos/UV_S4Hindcasts_0108_levels_24Hourly_'
winds = []
for Y in np.arange(1993, 2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib', chunks={'step': 10})
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		ds = ds.sel(**{'isobaricInhPa': [50]}).squeeze()
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time', skipna='True')
		winds.append(ds)
winds = xr.concat(winds, dim='year')

winds.to_netcdf('~/datos/data/winds50_djf_93-17.nc4')

RUTA = '~/datos/UV_S4Hindcasts_50hPa_2904hs_'
winds = []
for Y in np.arange(1993, 2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib', chunks={'step':10})
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		ds = ds.sel(**{'isobaricInhPa': [50]}).squeeze()
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time', skipna='True')
		winds.append(ds)
winds = xr.concat(winds, dim='year')
winds.to_netcdf('~/datos/data/winds50_ason_93-17.nc4')


