#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plots
import os
from windspharm.xarray import VectorWind
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '/storage/shared/glusterfs/acrcc/users/vg140344/data/'
ds = xr.open_dataset(RUTA + 'winds200.nc', chunks={'longitude':10})
ds = ds.stack(realiz=['year', 'number'])
ds = ds.transpose('month', 'realiz', 'latitude', 'longitude')
ds = ds.drop(['IC', 'isobaricInhPa'])
ds = ds.reset_index('realiz')
month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']
S_monthly = []
uchi_monthly = []
vchi_monthly = []
#for i in np.arange(0,4):
#	#create vectorwind instance
#	w = VectorWind(ds['u'][i, :, :, :], ds['v'][i, :, :, :])
#	eta = w.absolutevorticity()
#	div = w.divergence()
#	uchi, vchi = w.irrotationalcomponent()
#	etax, etay = w.gradient(eta)
#	etax.attrs['units'] = 'm**-1 s**-1'
#	etay.attrs['units'] = 'm**-1 s**-1'
#	# Combine the components to form the Rossby wave source term.
#	S = eta * -1. * div - (uchi * etax + vchi * etay)
#	S *= 1e11 
#	S_monthly.append(S)
#	uchi_monthly.append(uchi)
#	vchi_monthly.append(vchi)
#
#S_monthly = xr.concat(S_monthly, dim='month')
#S_monthly['month'] = np.array([8, 9, 10, 11])
#S_monthly.to_netcdf(RUTA + 'monthly_RWS.nc') 
#uchi_monthly = xr.concat(uchi_monthly, dim='month')
#uchi_monthly['month'] = np.array([8, 9, 10, 11])
#uchi_monthly.to_netcdf(RUTA + 'monthly_uchi.nc') 
#vchi_monthly = xr.concat(vchi_monthly, dim='month')
#vchi_monthly['month'] = np.array([8, 9, 10, 11])
#vchi_monthly.to_netcdf(RUTA + 'monthly_vchi.nc') 

S_seasonal = []
uchi_seasonal = []
vchi_seasonal = []

for i in np.arange(0,2):
	ds1 = ds.isel(month=range(i, i +3)).mean(dim='month')
	#create vectorwind instance
	w = VectorWind(ds1['u'][:, :, :], ds1['v'][:, :, :])
	eta = w.absolutevorticity()
	div = w.divergence()
	uchi, vchi = w.irrotationalcomponent()
	etax, etay = w.gradient(eta)
	etax.attrs['units'] = 'm**-1 s**-1'
	etay.attrs['units'] = 'm**-1 s**-1'
	# Combine the components to form the Rossby wave source term.
	#S = eta * -1. * div - (uchi * etax + vchi * etay)
	#S *= 1e11 
	#S_seasonal.append(S)
	uchi_seasonal.append(uchi)
	vchi_seasonal.append(vchi)

#S_seasonal = xr.concat(S_seasonal, dim='season')
#S_seasonal['season'] = np.array(['ASO', 'SON'])
#S_seasonal.to_netcdf(RUTA + 'seasonal_RWS.nc') 
uchi_seasonal = xr.concat(uchi_seasonal, dim='season')
uchi_seasonal['season'] = np.array(['ASO', 'SON'])
uchi_seasonal.to_netcdf(RUTA + 'seasonal_uchi.nc') 
vchi_seasonal = xr.concat(vchi_seasonal, dim='season')
vchi_seasonal['season'] = np.array(['ASO', 'SON'])
vchi_seasonal.to_netcdf(RUTA + 'seasonal_vchi.nc') 

