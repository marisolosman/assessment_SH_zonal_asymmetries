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
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_WINDS = 'fogt/winds200_aug_feb.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
winds = xr.open_dataset(PATH_DATA + FILE_WINDS)
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
winds_clm = winds.mean(dim='realiz')
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)


#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.90, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.10, dim='dim_0', interpolation='linear')

# compute EN-LA composites conditioned on PV anomalies
ninio34_WPV = ninio34.sel(dim_0 = index_SPV_upper.values)
hgt_WPV = hgt.sel(realiz = index_SPV_upper.values)
winds_WPV = winds.sel(realiz = index_SPV_upper.values)
ninio34_SPV = ninio34.sel(dim_0 = index_SPV_lower.values)
hgt_SPV = hgt.sel(realiz = index_SPV_lower.values)
winds_SPV = winds.sel(realiz = index_SPV_lower.values)

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during weak PoV
index_ninio_WPV = ninio34_WPV.ninio34_index >= ninio34_WPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_WPV = ninio34_WPV.ninio34_index <= ninio34_WPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_WPV = np.logical_and(ninio34_WPV.ninio34_index < ninio34_WPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_index > ninio34_WPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during strong PoV
index_ninio_SPV = ninio34_SPV.ninio34_index >= ninio34_SPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_SPV = ninio34_SPV.ninio34_index <= ninio34_SPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_SPV = np.logical_and(ninio34_SPV.ninio34_index < ninio34_SPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_index > ninio34_SPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))


month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_ninio_WPV = np.mean(hgt_WPV.z.values[i, index_ninio_WPV.values, :, :], axis=0)
	var_normal_WPV = np.mean(hgt_WPV.z.values[i, index_normal_WPV.values, :, :], axis=0)
	px_ninio_WPV, py_ninio_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_WPV-var_normal_WPV, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var_ninia_WPV = np.mean(hgt_WPV.z.values[i, index_ninia_WPV.values, :, :], axis=0)	
	px_ninia_WPV, py_ninia_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_WPV-var_normal_WPV, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var_ninio_SPV = np.mean(hgt_SPV.z.values[i, index_ninio_SPV.values, :, :], axis=0)
	var_normal_SPV = np.mean(hgt_SPV.z.values[i, index_normal_SPV.values, :, :], axis=0)
	px_ninio_SPV, py_ninio_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_SPV-var_normal_SPV, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)

	var_ninia_SPV = np.mean(hgt_SPV.z.values[i, index_ninia_SPV.values, :, :], axis=0)	
	px_ninia_SPV, py_ninia_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_SPV-var_normal_SPV, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var_ninio_all = np.mean(hgt.z.values[i, index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.mean(hgt.z.values[i, index_normal_all.values, :, :], axis=0)
	px_ninio_all, py_ninio_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_all-var_normal_all, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var_ninia_all = np.mean(hgt.z.values[i, index_ninia_all.values, :, :], axis=0)	
	px_ninia_all, py_ninia_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_all-var_normal_all, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var ={'z1': var_ninio_all-var_normal_all, 'px1': px_ninio_all, 'py1': py_ninio_all,
	      'z2': var_ninia_all-var_normal_all, 'px2': px_ninia_all, 'py2': py_ninia_all,
	      'z3': var_ninio_WPV-var_normal_WPV, 'px3': px_ninio_WPV, 'py3': py_ninio_WPV,
	      'z4': var_ninia_WPV-var_normal_WPV, 'px4': px_ninia_WPV, 'py4': py_ninia_WPV,
	      'z5': var_ninio_SPV-var_normal_SPV, 'px5': px_ninio_SPV, 'py5': py_ninio_SPV,
	      'z6': var_ninia_SPV-var_normal_SPV, 'px6': px_ninia_SPV, 'py6': py_ninia_SPV}
	tit = 'Composites S4 Z* and Plumb Fluxes 200hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'z200_pf_composites_ENSO_' + month[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoVZPF(var, hgt.latitude, hgt.longitude, tit, filename)

for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	winds_clm_s = winds_clm.isel(month=range(i, i+3)).mean(dim='month')
	hgt_s_WPV = hgt_s.sel(realiz=index_SPV_upper.values)
	hgt_s_SPV = hgt_s.sel(realiz=index_SPV_lower.values)
	var_ninio_WPV = np.mean(hgt_s_WPV.z.values[index_ninio_WPV.values, :, :], axis=0)
	px_ninio_WPV, py_ninio_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_WPV-var_normal_WPV, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var_normal_WPV = np.mean(hgt_s_WPV.z.values[index_normal_WPV.values, :, :], axis=0)
	var_ninia_WPV = np.mean(hgt_s_WPV.z.values[index_ninia_WPV.values, :, :], axis=0)
	px_ninia_WPV, py_ninia_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninia_WPV-var_normal_WPV, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	
	var_ninio_SPV = np.mean(hgt_s_SPV.z.values[index_ninio_SPV.values, :, :], axis=0)
	var_normal_SPV = np.mean(hgt_s_SPV.z.values[index_normal_SPV.values, :, :], axis=0)
	px_ninio_SPV, py_ninio_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_SPV-var_normal_SPV, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var_ninia_SPV = np.mean(hgt_s_SPV.z.values[index_ninia_SPV.values, :, :], axis=0)
	px_ninia_SPV, py_ninia_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_SPV-var_normal_SPV, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var_ninio_all = np.mean(hgt_s.z.values[index_ninio_all.values, :, :], axis=0)
	px_ninio_all, py_ninio_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_all-var_normal_all, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	var_normal_all = np.mean(hgt_s.z.values[index_normal_all.values, :, :], axis=0)
	var_ninia_all = np.mean(hgt_s.z.values[index_ninia_all.values, :, :], axis=0)
	px_ninia_all, py_ninia_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_all-var_normal_all, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
	tit = 'Composites S4 Z* and Plumb Fluxes 200hPa Conditioned - SPoV - ' + seas[i]
	filename = FIG_PATH + 'z200_pf_composites_ENSO_' + seas[i] +'_SPoV.png'
	var ={'z1': var_ninio_all-var_normal_all, 'px1': px_ninio_all, 'py1': py_ninio_all,
	      'z2': var_ninia_all-var_normal_all, 'px2': px_ninia_all, 'py2': py_ninia_all,
	      'z3': var_ninio_WPV-var_normal_WPV, 'px3': px_ninio_WPV, 'py3': py_ninio_WPV,
	      'z4': var_ninia_WPV-var_normal_WPV, 'px4': px_ninia_WPV, 'py4': py_ninia_WPV,
	      'z5': var_ninio_SPV-var_normal_SPV, 'px5': px_ninio_SPV, 'py5': py_ninio_SPV,
	      'z6': var_ninia_SPV-var_normal_SPV, 'px6': px_ninia_SPV, 'py6': py_ninia_SPV}
	plots.PlotEnsoCompositesPoVZPF(var, hgt.latitude, hgt.longitude, tit, filename)

