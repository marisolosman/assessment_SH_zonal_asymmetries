#compute monthly PoV index for S4
import numpy as np
import xarray as xr
import eofdata
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/'
ds = xr.open_dataset(RUTA + 'monthly_hgt50_aug_feb_inc_2002.nc')
ds = ds.sel(**{'latitude': slice(-60, -90)}).mean(dim=['longitude', 'latitude']).compute()

[lamb, v, PC] = eofdata.eofdata(ds.z.values, 3)
print(v[:, 0])
SPV_monthly_index = -1* PC[0, :]

ds_new = xr.Dataset({'SPV_index': xr.DataArray(SPV_monthly_index)})

ds_new.to_netcdf(RUTA + 'fogt/SPV_index_inc_2002.nc4')


