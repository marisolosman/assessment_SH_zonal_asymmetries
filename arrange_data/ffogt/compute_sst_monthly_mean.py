#compute monthly sst for S4 forecasts
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '/storage/silver/acrcc/ys916780/Data_for_Sol/output'

sst = []
#sst3 = [] #ninio 3 5N-5S 150W-90W 
#sst4 = [] #ninio 4 5N-5S 160E-150W
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + '_' + str(Y) + '_sst', engine='cfgrib', chunks={'step':10},
				     backend_kwargs={'indexpath':''})
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step':'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		sst.append(ds)
#
sst =xr.concat(sst, dim='year')
sst.to_netcdf('~/datos/data/fogt/sst_s4_aug_nov.nc4')
RUTA = '/storage/silver/acrcc/vg140344/'

sst = []
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + 'SST_S4Hindcast_0108_24Hourly_' + str(Y) + '0801.grib',
					 engine='cfgrib', chunks={'step':10},
				     backend_kwargs={'indexpath':''})
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step':'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		sst.append(ds)

sst =xr.concat(sst, dim='year')
sst.to_netcdf('/home/users/vg140344/datos/data/fogt/sst_s4_dec_feb.nc4')

