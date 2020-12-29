#composites of z200 and plumb fluxes for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import doharmonic

#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT_S4 = 'monthly_hgt50_aug_feb.nc4'

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

#hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4, chunks={'latitude':10})
#hgt = hgt.sel(latitude=slice(-20, -90)).compute()
#hgt = hgt - hgt.mean(dim='longitude')
#harm1 = np.zeros_like(hgt.z.values)
#harm2 = np.zeros_like(hgt.z.values)
#[nmonth, nrealiz, nlat, nlon] = np.shape(hgt.z.values)
#for i in np.arange(7):
#	for j in np.arange(nrealiz):
#		for k in np.arange(nlat):
#			varianza, y, c = doharmonic.doharmonic(hgt.z.values[i, j, k, :], 1, 3)
#			harm1[i, j, k, :] = y[:, 0]
#			harm2[i, j, k, :] = y[:, 1]
#ds = xr.Dataset({'Harmonic1': (['month', 'realiz', 'latitude', 'longitude'], harm1),
#		 'Harmonic2': (['month', 'realiz', 'latitude', 'longitude'], harm2)},
#		 coords={'month': month,
#			 'realiz': np.arange(nrealiz),
#			 'longitude': ('longitude', hgt.longitude.values),
#			 'latitude': ('latitude', hgt.latitude.values)})
#ds.to_netcdf(PATH_DATA_2 + 'hgt50_monthly_harmonics.nc4')
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4, chunks={'latitude':10})
hgt = hgt.sel(latitude=slice(-20, -90)).compute()
hgt = hgt - hgt.mean(dim='longitude')
hgt = hgt.mean(dim='realiz')
harm1 = np.zeros_like(hgt.z.values)
harm2 = np.zeros_like(hgt.z.values)
[nmonth, nlat, nlon] = np.shape(hgt.z.values)
for i in np.arange(7):
	for k in np.arange(nlat):
		varianza, y, c = doharmonic.doharmonic(hgt.z.values[i, k, :], 1, 2)
		harm1[i, k, :] = y[:, 0]
		harm2[i, k, :] = y[:, 1]
ds = xr.Dataset({'Harmonic1': (['month', 'latitude', 'longitude'], harm1),
		 'Harmonic2': (['month', 'latitude', 'longitude'], harm2)},
		 coords={'month': month,
			 'longitude': ('longitude', hgt.longitude.values),
			 'latitude': ('latitude', hgt.latitude.values)})
ds.to_netcdf(PATH_DATA_2 + 'hgt50_monthly_climatology_harmonics.nc4')


### Seasonal harmonics
#harm1 = np.zeros([5, nrealiz, nlat, nlon])
#harm2 = np.zeros([5, nrealiz, nlat, nlon])
#varianza = np.zeros([5, 2, nrealiz, nlat])
#
#for i in np.arange(0, 5):
#	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
#	for j in np.arange(nrealiz):
#		for k in np.arange(nlat):
#			exp_var, y, c = doharmonic.doharmonic(hgt_s.z.values[j, k, :], 1, 3)
#			harm1[i, j, k, :] = y[:, 0]
#			harm2[i, j, k, :] = y[:, 1]
#			varianza [i, :, j, k] = exp_var[0:2]
#ds = xr.Dataset({'Harmonic1': (['season', 'realiz', 'latitude', 'longitude'], harm1),
#		 'Harmonic2': (['season', 'realiz', 'latitude', 'longitude'], harm2),
#		 'Varianza': (['season', 'harmonic', 'realiz', 'latitude'], varianza)},
#		 coords={'season': seas,
#			 'realiz': np.arange(nrealiz),
#			 'longitude': ('longitude', hgt.longitude.values),
#			 'latitude': ('latitude', hgt.latitude.values),
#			 'harmonic': ('harmonic', [1, 2])})
#ds.to_netcdf(PATH_DATA_2 + 'hgt50_seasonal_harmonics.nc4')
#
#		
