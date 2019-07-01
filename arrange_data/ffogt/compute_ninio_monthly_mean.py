#compute monthly el ninio indexes
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '/storage/shared/glusterfs/acrcc/users/ys916780/Data_for_Sol/output'

#sst = []
#for Y in np.arange(1981,2018):
#	if Y != 2002:
#		ds = xr.open_dataset(RUTA + '_' + str(Y) + '_sst', engine='cfgrib', chunks={'step':10},
#				     backend_kwargs={'indexpath':''})
#		ds = ds.sel(**{'latitude':slice(5, -5), 'longitude':slice(190, 240)}).mean(
#			dim=['longitude', 'latitude']).compute()
#		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
#		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
#		ds = ds.rename({'step':'time'}).set_coords(['time'])
#		ds = ds.groupby('time.month').mean(dim='time')
#		sst.append(ds)
#
#sst =xr.concat(sst, dim='year')
#
#sst.to_netcdf('./data/sst.nc')

#ninio 3 5N-5S 150W-90W 
sst = []
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + '_' + str(Y) + '_sst', engine='cfgrib', chunks={'step':10},
				     backend_kwargs={'indexpath':''})
		ds = ds.sel(**{'latitude':slice(5, -5), 'longitude':slice(210, 270)}).mean(
			dim=['longitude', 'latitude']).compute()
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step':'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		sst.append(ds)

sst =xr.concat(sst, dim='year')

sst.to_netcdf('~/datos/data/sst_ninio3.nc')

#ninio 4 5N-5S 160E-150W
sst = []
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + '_' + str(Y) + '_sst', engine='cfgrib', chunks={'step':10},
				     backend_kwargs={'indexpath':''})
		ds = ds.sel(**{'latitude':slice(5, -5), 'longitude':slice(160, 210)}).mean(
			dim=['longitude', 'latitude']).compute()
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step':'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		sst.append(ds)

sst =xr.concat(sst, dim='year')

sst.to_netcdf('~/datos/data/sst_ninio4.nc')


