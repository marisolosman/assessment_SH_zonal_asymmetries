#composites of z200 and plumb fluxes for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots_stereo as plots
import random
import doharmonic

def TestCorrelation(var1, var2):
	var1 = np.ravel(var1)
	var2 = np.ravel(var2)
	correlacion = np.corrcoef(var1, var2)[0, 1]
	distrib = np.empty([10000])
	for i in range(10000):
		j = [random.randint(0, var1.shape[0] - 1) for k in range(var1.shape[0])]
		distrib[i] = np.corrcoef(var1[j], var2[j])[0, 1]
	distrib = np.sort(distrib)
	return correlacion, distrib[499], distrib[9499]
def ComputeAsymmetry(field, pattern):
	pattern = np.ravel(pattern)
	pattern = np.tile(pattern[np.newaxis, :], (field.shape[0], 1))
	field = np.reshape(field,[field.shape[0], field.shape[1]*field.shape[2]])
	index = np.squeeze(np.sum(field * pattern, axis=1))
	return index
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')

index_SPV_normal = np.logical_and(PV_index.SPV_index > PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'), PV_index.SPV_index < PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'))

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))
#enso during weak PoV

index_ninio_WPV = np.logical_and(index_ninio_all.values, index_SPV_upper.values)
index_ninia_WPV = np.logical_and(index_ninia_all.values, index_SPV_upper.values)

#enso during strong PoV
index_ninio_SPV = np.logical_and(index_ninio_all.values, index_SPV_lower.values)
index_ninia_SPV = np.logical_and(index_ninia_all.values, index_SPV_lower.values)

index_normal = np.logical_and(index_SPV_normal, index_normal_all)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

#hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4, chunks={'latitude':10})
#hgt = hgt.sel(latitude=slice(-20, -90)).compute()
#hgt = hgt - hgt.mean(dim='longitude')
#[nmonth, nrealiz, nlat, nlon] = np.shape(hgt.z.values)

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

harms = xr.open_dataset(PATH_DATA_2 + 'hgt200_monthly_harmonics.nc4')# , chunks={'latitude':10})

for i in np.arange(0, 7):
	for j in np.arange(5):
		print(np.mean(harms.Varianza.values[j, i, :, :], axis=0))
		var_normal = np.mean(harms.Harmonics.values[j, i, :, :, :], axis=0)
		var_ninio_WPV = np.mean(harms.Harmonics.values[j, i, index_ninio_WPV, :, :], axis=0)
		var_ninia_WPV = np.mean(harms.Harmonics.values[j, i, index_ninia_WPV, :, :], axis=0)	
		var_ninio_SPV = np.mean(harms.Harmonics.values[j, i, index_ninio_SPV, :, :], axis=0)
		var_ninia_SPV = np.mean(harms.Harmonics.values[j, i, index_ninia_SPV, :, :], axis=0)	
		var_ninio_all = np.mean(harms.Harmonics.values[j, i, index_ninio_all, :, :], axis=0)
		var_ninia_all = np.mean(harms.Harmonics.values[j, i, index_ninia_all, :, :], axis=0)	
		tit = 'Composites S4 Z* 200hPa Wave-' + str(j + 1) + 'Conditioned - SPoV - ' + month[i]
		filename = FIG_PATH + 'z200W' + str(j + 1) + '_composites_ENSO_' + month[i] +'_SPoV.png'
		plots.PlotEnsoCompositesPoV(var_ninio_all-var_normal, var_ninia_all-var_normal,
					    var_ninio_WPV-var_normal, var_ninia_WPV-var_normal,
					    var_ninio_SPV-var_normal, var_ninia_SPV-var_normal,
					    var_normal, harms.latitude, harms.longitude, tit, filename)
harms = xr.open_dataset(PATH_DATA_2 + 'hgt200_seasonal_harmonics.nc4')
#had a typo and harmonics are under variable Harmonic1
for i in np.arange(0, 5):
	for j in np.arange(0, 5):
		var_normal = np.mean(harms.Harmonic1.values[j, i, :, :, :], axis=0)
		var_ninio_WPV = np.mean(harms.Harmonic1.values[j, i, index_ninio_WPV, :, :], axis=0)
		var_ninia_WPV = np.mean(harms.Harmonic1.values[j, i, index_ninia_WPV, :, :], axis=0)	
		var_ninio_SPV = np.mean(harms.Harmonic1.values[j, i, index_ninio_SPV, :, :], axis=0)
		var_ninia_SPV = np.mean(harms.Harmonic1.values[j, i, index_ninia_SPV, :, :], axis=0)	
		var_ninio_all = np.mean(harms.Harmonic1.values[j, i, index_ninio_all, :, :], axis=0)
		var_ninia_all = np.mean(harms.Harmonic1.values[j, i, index_ninia_all, :, :], axis=0)
		tit = 'Composites S4 Z* 200hPa Wave-' + str(j + 1) +  'Conditioned - SPoV - ' + seas[i]
		filename = FIG_PATH + 'z200W' + str(j + 1) + '_composites_ENSO_' + seas[i] +'_SPoV.png'
		plots.PlotEnsoCompositesPoV(var_ninio_all-var_normal, var_ninia_all-var_normal,
					   var_ninio_WPV-var_normal, var_ninia_WPV-var_normal,
					   var_ninio_SPV-var_normal, var_ninia_SPV-var_normal,
					   var_normal, harms.latitude, harms.longitude, tit, filename)
	
