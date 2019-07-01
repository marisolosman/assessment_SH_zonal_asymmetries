#compute monthly el ninio 3.4 index
import numpy as np
import xarray as xr
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/'
#ds2 = xr.concat([ds, ds1], dim='month')
#ds2 = ds2.loc[{'month': [8, 9, 10, 11, 12, 1, 2]}]
hgt50 = xr.open_dataset(RUTA + 'monthly_hgt50_aug_feb.nc')
#

hgt50 = hgt50.sel(**{'latitude':slice(-60, -90)}).mean(dim=['longitude','latitude']).compute()

hgt.to_netcdf(RUTA + 'PV_monthly_s4.nc4')


