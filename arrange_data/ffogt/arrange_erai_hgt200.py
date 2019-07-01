#compute monthly el ninio 3.4 index
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '/storage/shared/glusterfs/acrcc/users/vg140344/ERAI/'
hgt = []
for Y in np.arange(1981, 2019):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + 'H200_erai_' + str(Y) + '.grib',
				     engine='cfgrib', backend_kwargs={'indexpath':''})
		ds = ds.drop(['number', 'isobaricInhPa', 'step'])
		hgt.append(ds)

hgt =xr.concat(hgt, dim='time')
hgt.time.values = np.arange(0, len(hgt.time.values))
hgt.to_netcdf('~/datos/data/hgt_erai_200.nc4')

