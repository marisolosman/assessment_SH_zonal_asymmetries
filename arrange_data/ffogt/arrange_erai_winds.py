#concatenate HGT at 200hPa and 50 hPa era interim data from 1981 to 2018
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/ERAI/'
#200 hPa
winds200 = []
for Y in np.arange(1981, 2019):
	ds = xr.open_dataset(RUTA + 'WINDS_erai_' + str(Y) + '.grib',
				engine='cfgrib', backend_kwargs={'indexpath':''})
	ds = ds.sel(isobaricInhPa=200)
	ds = ds.drop(['number', 'isobaricInhPa', 'step'])
	winds200.append(ds)

winds200 =xr.concat(winds200, dim='time')
winds200.time.values = np.arange(0, len(winds200.time.values))
winds200.to_netcdf('~/datos/data/winds_erai_200.nc4')
#50hPa
winds50 = []
for Y in np.arange(1981, 2019):
	ds = xr.open_dataset(RUTA + 'WINDS_erai_' + str(Y) + '.grib',
				engine='cfgrib', backend_kwargs={'indexpath':''})
	ds = ds.sel(isobaricInhPa=50)
	ds = ds.drop(['number', 'isobaricInhPa', 'step'])
	winds50.append(ds)

winds50 =xr.concat(winds50, dim='time')
winds50.time.values = np.arange(0, len(winds50.time.values))
winds50.to_netcdf('~/datos/data/winds_erai_50.nc4')

