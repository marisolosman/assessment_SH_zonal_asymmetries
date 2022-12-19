#compute monthly el ninio indexes for S4 forecasts
import numpy as np
import xarray as xr
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/fogt/'


ds = xr.open_dataset(RUTA + 'sst_s4_aug_feb_inc_2002.nc4')

ds_34 = ds.sel(**{'latitude': slice(5, -5), 'longitude': slice(190, 240)}).mean(
		dim=['longitude', 'latitude']).compute()
ds_3 = ds.sel(**{'latitude': slice(5, -5), 'longitude': slice(210, 270)}).mean(
	dim=['longitude', 'latitude']).compute()
ds_4 = ds.sel(**{'latitude': slice(5, -5), 'longitude': slice(160, 210)}).mean(
	dim=['longitude', 'latitude']).compute()

ds_34.to_netcdf('/home/users/vg140344/datos/data/fogt/sst_ninio34_s4_aug_feb_inc_2002.nc4')
ds_3.to_netcdf('/home/users/vg140344/datos/data/fogt/sst_ninio3_s4_aug_feb_inc_2002.nc4')
ds_4.to_netcdf('/home/users/vg140344/datos/data/fogt/sst_ninio4_s4_aug_feb_inc_2002.nc4')

