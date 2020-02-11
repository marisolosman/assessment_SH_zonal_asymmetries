#concatenate HGT at 200hPa and 50 hPa era interim data from 1981 to 2018
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/ERAI/'
#200 hPa
hgt = []
for Y in np.arange(1981, 2019):
	ds = xr.open_dataset(RUTA + 'H200_erai_' + str(Y) + '.grib',
				engine='cfgrib', backend_kwargs={'indexpath':''})
	ds = ds.drop(['number', 'isobaricInhPa', 'step'])
	hgt.append(ds)

hgt =xr.concat(hgt, dim='time')
hgt.time.values = np.arange(0, len(hgt.time.values))
hgt.z.values = hgt.z.values /10
hgt.to_netcdf('~/datos/data/hgt_erai_200.nc4')
#50hPa
hgt = []
for Y in np.arange(1981, 2019):
	ds = xr.open_dataset(RUTA + 'HGT50_erai_' + str(Y) + '.grib',
				engine='cfgrib', backend_kwargs={'indexpath':''})
	ds = ds.drop(['number', 'isobaricInhPa', 'step'])
	hgt.append(ds)

hgt =xr.concat(hgt, dim='time')
hgt.z.values = hgt.z.values /10
hgt.time.values = np.arange(0, len(hgt.time.values))
hgt.to_netcdf('~/datos/data/hgt_erai_50.nc4')

