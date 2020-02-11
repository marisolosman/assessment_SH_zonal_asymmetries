import datetime
import numpy as np
import xarray as xr
from eofs.standard import Eof
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '~/datos/data/' 
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
hgt_s4 = xr.open_dataset(RUTA + FILE_HGT_S4)
hgt_s4 = hgt_s4.sel(**{'latitude':slice(-20, -90)})
FILE_HGT_ERAI = 'hgt_erai_200.nc4'
hgt_erai = xr.open_dataset(RUTA + FILE_HGT_ERAI)
hgt_erai.time.values = hgt_erai.valid_time.values
# Discard data from 2002-2003
hgt_erai = hgt_erai.sel(time= np.logical_or(hgt_erai.time.values <=np.datetime64('2002-07-31'), hgt_erai.time.values>=np.datetime64('2003-08-01')))
hgt_erai = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')})
hgt_erai = hgt_erai.sel(**{'latitude':slice(-20, -90)})

season = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lmonth = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec']
#loop over seasons, selec data and performs boxplot
SAM_erai = np.zeros([5, 36])
SAM_s4 = np.zeros([5, len(hgt_s4.realiz.values)])
eof_erai = np.zeros([5, len(hgt_erai.latitude.values)])
eof_s4 = np.zeros([5, len(hgt_s4.latitude.values)])
sign_s4 = np.array([1, 1, 1, 1, -1])
sign_erai = np.array([1, -1, -1, -1, 1])
for i in np.arange(0, 5):
	aux = hgt_erai['z'].resample(time='QS-' + lmonth[i]).mean(dim='time',skipna=True)
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	aux = aux.sel(time= np.logical_and(aux['time.month'] == mes, aux['time.year']!=2002))
	X_zm = aux.mean(dim='longitude')
	X_an = X_zm - X_zm.mean(dim='time')
	solver = Eof(X_an.values)
	pcs = solver.pcs(npcs=1, pcscaling=1)
	eof_erai[i, :] = solver.eofs(neofs=1)[0,:]
	SAM_erai[i, :] = sign_erai[i] * pcs[:, 0]
	hgt_s4_smean = np.nanmean(np.nanmean(hgt_s4.z.values[i:i + 3, :, :, :], axis=0), axis=2)
	hgt_s4_smean = hgt_s4_smean - np.nanmean(hgt_s4_smean, axis=0)
	solver = Eof(hgt_s4_smean)
	pcs = solver.pcs(npcs=1, pcscaling=1)
	eof_s4[i, :] = solver.eofs(neofs=1)[0,:]
	SAM_s4[i, :] = sign_s4[i] * pcs[:, 0]

time = np.concatenate([np.arange(1981,2002), np.arange(2003,2018)])
ds = xr.Dataset({'SAM_index': xr.DataArray(SAM_erai, coords=[('season', season),('year', time)])})

ds.to_netcdf(RUTA + 'fogt/SAM_index_erai.nc4')

ds1 = xr.Dataset({'SAM_index': xr.DataArray(SAM_s4, coords=[('season', season),('realiz', np.arange(SAM_s4.shape[1]))])})

ds1.to_netcdf(RUTA + 'fogt/SAM_index_s4.nc4')
