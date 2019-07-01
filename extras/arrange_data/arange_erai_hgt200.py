#compute monthly el ninio 3.4 index
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

RUTA = '/storage/shared/glusterfs/acrcc/users/vg140344/ERAI/'
hgt = []
for Y in np.arange(1981, 2019):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + 'H200_erai_' + str(Y) + '.grib',
				     engine='cfgrib', backend_kwargs={'indexpath':''})
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		#ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step':'time'}).set_coords(['time'])
		hgt.append(ds)

hgt =xr.concat(hgt, dim='time')
print("guardo")
hgt.to_netcdf('./data/hgt_erai_200.nc')

