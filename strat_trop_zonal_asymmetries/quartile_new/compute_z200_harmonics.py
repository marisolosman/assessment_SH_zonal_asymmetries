#decomposed z200* into the first 5 harmonics
import numpy as np
import xarray as xr
import os
import doharmonic

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4, chunks={'latitude':10})
hgt = hgt.sel(latitude=slice(-20, -90)).compute()
hgt = hgt - hgt.mean(dim='longitude')
[nmonth, nrealiz, nlat, nlon] = np.shape(hgt.z.values)

#harms = np.zeros([5, nmonth, nrealiz, nlat, nlon])
#varianza = np.zeros([5, nmonth, nrealiz, nlat])
#for i in np.arange(7):
#	for j in np.arange(nrealiz):
#		for k in np.arange(nlat):
#			exp_var, y, c = doharmonic.doharmonic(hgt.z.values[i, j, k, :], 1, 5)
#			harms[:, i, j, k, :] = np.transpose(y[:, :])
#			varianza [:, i, j, k] = exp_var[:]
#
#ds = xr.Dataset({'Harmonics': (['harmonic', 'month', 'realiz', 'latitude', 'longitude'], harms),
#		 'Varianza': (['harmonic', 'month', 'realiz', 'latitude'], varianza)},
#		 coords={'month': month,
#			 'realiz': np.arange(nrealiz),
#			 'longitude': ('longitude', hgt.longitude.values),
#			 'latitude': ('latitude', hgt.latitude.values),
#			  'harmonic': np.arange(5)})
#ds.to_netcdf(PATH_DATA_2 + 'hgt200_monthly_harmonics.nc4')
harms = np.zeros([2, nmonth, nlat, nlon])
hgt = hgt.mean(dim='realiz')
varianza = np.zeros([2, nmonth, nlat])
for i in np.arange(7):
	for k in np.arange(nlat):
		exp_var, y, c = doharmonic.doharmonic(hgt.z.values[i, k, :], 1, 2)
		harms[:, i, k, :] = np.transpose(y[:, :])
		varianza [:, i, k] = exp_var[:]
ds = xr.Dataset({'Harmonics': (['harmonic', 'month', 'latitude', 'longitude'], harms),
		 'Varianza': (['harmonic', 'month', 'latitude'], varianza)},
		 coords={'month': month,
			 'longitude': ('longitude', hgt.longitude.values),
			 'latitude': ('latitude', hgt.latitude.values),
			  'harmonic': np.arange(2)})
ds.to_netcdf(PATH_DATA_2 + 'hgt200_monthly_clim_harmonics.nc4')
#
###Seasonal Z*
#harms = np.zeros([5, 5, nrealiz, nlat, nlon])
#varianza = np.zeros([5, 5, nrealiz, nlat])
#
#for i in np.arange(0, 5):
#	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
#	for j in np.arange(nrealiz):
#		for k in np.arange(nlat):
#			exp_var, y, c = doharmonic.doharmonic(hgt_s.z.values[j, k, :], 1, 5)
#			harms[:, i, j, k, :] = np.transpose(y)
#			varianza[:, i, j, k] = exp_var
#ds = xr.Dataset({'Harmonic1': (['harmonic', 'season', 'realiz', 'latitude', 'longitude'], harms),
#		 'Varianza': (['harmonic', 'season', 'realiz', 'latitude'], varianza)},
#		 coords={'season': seas,
#			 'realiz': np.arange(nrealiz),
#			 'longitude': ('longitude', hgt.longitude.values),
#			 'latitude': ('latitude', hgt.latitude.values),
#			 'harmonic': ('harmonic', np.arange(5))})
#ds.to_netcdf(PATH_DATA_2 + 'hgt200_seasonal_harmonics.nc4')
#
