#composites of z200 and plumb fluxes for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots
import random
import plumb_flux
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
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_WINDS = 'fogt/monthly_winds50_aug_feb.nc4'
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
#search for years with strong PV
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

winds = xr.open_dataset(PATH_DATA + FILE_WINDS, chunks={'latitude':10})
winds = winds.drop_vars('v')
winds = winds.sel(latitude=slice(20, -90)).compute()
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')

for i in np.arange(0, 7):
	var_normal = np.nanmean(winds.u.values[i, :, :, :], axis=0)
	var_ninio_WPV = np.nanmean(winds.u.values[i, index_ninio_WPV, :, :], axis=0)
	var_ninia_WPV = np.nanmean(winds.u.values[i, index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.nanmean(winds.u.values[i, index_ninio_SPV, :, :], axis=0)
	var_ninia_SPV = np.nanmean(winds.u.values[i, index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.nanmean(winds.u.values[i, index_ninio_all, :, :], axis=0)
	var_ninia_all = np.nanmean(winds.u.values[i, index_ninia_all, :, :], axis=0)	
	var ={'z1': var_ninio_all-var_normal,# 'px1': px_ninio_all, 'py1': py_ninio_all,
	      'z2': var_ninia_all-var_normal,# 'px2': px_ninia_all, 'py2': py_ninia_all,
	      'z3': var_ninio_WPV-var_normal,# 'px3': px_ninio_WPV, 'py3': py_ninio_WPV,
	      'z4': var_ninia_WPV-var_normal,# 'px4': px_ninia_WPV, 'py4': py_ninia_WPV,
	      'z5': var_ninio_SPV-var_normal,# 'px5': px_ninio_SPV, 'py5': py_ninio_SPV,
	      'z6': var_ninia_SPV-var_normal}# 'px6': px_ninia_SPV, 'py6': py_ninia_SPV}
	tit = 'Composites S4 U 50hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'U50_composites_ENSO_' + month[i] +'_SPoV_q.png'
	plots.PlotEnsoCompositesPoVU(var, winds.latitude, winds.longitude, tit, filename)

for i in np.arange(0, 5):
	hgt_s = winds.isel(month=range(i, i+3)).mean(dim='month')
	var_ninio_WPV = np.nanmean(hgt_s.u.values[index_ninio_WPV, :, :], axis=0)
	var_normal = np.nanmean(hgt_s.u.values[:, :, :], axis=0)
	var_ninia_WPV = np.nanmean(hgt_s.u.values[index_ninia_WPV, :, :], axis=0)
	var_ninio_SPV = np.nanmean(hgt_s.u.values[index_ninio_SPV, :, :], axis=0)
	var_ninia_SPV = np.nanmean(hgt_s.u.values[index_ninia_SPV, :, :], axis=0)
	var_ninio_all = np.nanmean(hgt_s.u.values[index_ninio_all, :, :], axis=0)
	var_ninia_all = np.nanmean(hgt_s.u.values[index_ninia_all, :, :], axis=0)
	tit = 'Composites S4 U 50hPa Conditioned - SPoV - ' + seas[i]
	filename = FIG_PATH + 'U50_composites_ENSO_' + seas[i] +'_SPoV_q.png'
	var ={'z1': var_ninio_all-var_normal,# 'px1': px_ninio_all, 'py1': py_ninio_all,
	      'z2': var_ninia_all-var_normal,# 'px2': px_ninia_all, 'py2': py_ninia_all,
	      'z3': var_ninio_WPV-var_normal,# 'px3': px_ninio_WPV, 'py3': py_ninio_WPV,
	      'z4': var_ninia_WPV-var_normal,# 'px4': px_ninia_WPV, 'py4': py_ninia_WPV,
	      'z5': var_ninio_SPV-var_normal,# 'px5': px_ninio_SPV, 'py5': py_ninio_SPV,
	      'z6': var_ninia_SPV-var_normal}# 'px6': px_ninia_SPV, 'py6': py_ninia_SPV}
	plots.PlotEnsoCompositesPoVU(var, winds.latitude, winds.longitude, tit, filename)

