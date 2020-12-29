#compute monthly PoV index S4
import numpy as np
import xarray as xr
import eofdata
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/'
ds = xr.open_dataset(RUTA + 'monthly_hgt50_aug_feb.nc4')
ds = ds.sel(**{'latitude': slice(-60, -90)}).mean(dim=['longitude', 'latitude']).compute()
ds.to_netcdf(RUTA + 'fogt/SPV_s4.nc4')


