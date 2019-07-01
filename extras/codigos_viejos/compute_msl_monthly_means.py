import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

RUTA = '/storage/shared/glusterfs/acrcc/users/ys916780/Data_for_Sol/output'
msl = []
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + '_' + str(Y) + '_slp', engine='cfgrib', chunks={'step':10},
				     backend_kwargs={'indexpath':''})
		ds = ds.sel(**{'latitude':slice(0, -30), 'longitude':slice(100, 300)}).compute()
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step':'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		msl.append(ds)

msl =xr.concat(msl, dim='year')

msl.to_netcdf('./data/msl.nc')


