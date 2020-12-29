#compute monthly means of u and v from S4 models
import numpy as np
import xarray as xr
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '/storage/silver/acrcc/ys916780/Data_for_Sol/output'
winds = []
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + '_' + str(Y)+ '_uv200', engine='cfgrib', chunks={'step': 10},
				     backend_kwargs={'indexpath': ''})
		ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		winds.append(ds)
winds =xr.concat(winds, dim='year')
winds.to_netcdf('~/datos/data/fogt/winds200_aug_nov.nc4')

RUTA = '~/datos/UV_S4Hindcasts_0108_200hPalevels_24Hourly_'
winds = []
for Y in np.arange(1981,1993):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib', chunks={'step': 10},
				     backend_kwargs={'indexpath': ''})
		ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		winds.append(ds)


RUTA = '~/datos/UV_S4Hindcasts_0108_levels_24Hourly_'
for Y in np.arange(1993,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib', chunks={'step': 10},
				     backend_kwargs={'indexpath': ''})
		ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
		ds = ds.sel(**{'isobaricInhPa': [200]}).squeeze()
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		winds.append(ds)

winds =xr.concat(winds, dim='year')
winds.to_netcdf('~/datos/data/fogt/winds200_dec_feb.nc4')


