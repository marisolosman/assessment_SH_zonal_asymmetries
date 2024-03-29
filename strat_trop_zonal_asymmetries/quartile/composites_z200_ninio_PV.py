#composites of EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')

#search for years with normal PV
index_SPV_normal = np.logical_and(PV_index.SPV_index > PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'), PV_index.SPV_index < PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'))



# compute EN-LA composites conditioned on PV anomalies
ninio34_WPV = ninio34.sel(dim_0 = index_SPV_upper.values)
hgt_WPV = hgt.sel(realiz = index_SPV_upper.values)
ninio34_SPV = ninio34.sel(dim_0 = index_SPV_lower.values)
hgt_SPV = hgt.sel(realiz = index_SPV_lower.values)

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))
#enso during weak PoV

index_ninio_WPV = np.logical_and(index_ninio_all.values, index_SPV_upper.values)
index_ninia_WPV = np.logical_and(index_ninia_all.values, index_SPV_upper.values)
index_normal_WPV = np.logical_and(index_normal_all.values, index_SPV_upper.values)

#enso during strong PoV
index_ninio_SPV = np.logical_and(index_ninio_all.values, index_SPV_lower.values)
index_ninia_SPV = np.logical_and(index_ninia_all.values, index_SPV_lower.values)
index_normal_SPV = np.logical_and(index_normal_all.values, index_SPV_lower.values)

index_normal = np.logical_and(index_normal_all.values, index_SPV_normal.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_ninio_WPV = np.mean(hgt.z.values[i, index_ninio_WPV, :, :], axis=0)
	var_normal_WPV = np.mean(hgt.z.values[i, index_normal_WPV, :, :], axis=0)
	var_ninia_WPV = np.mean(hgt.z.values[i, index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.mean(hgt.z.values[i, index_ninio_SPV, :, :], axis=0)
	var_normal_SPV = np.mean(hgt.z.values[i, index_normal_SPV, :, :], axis=0)
	var_ninia_SPV = np.mean(hgt.z.values[i, index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.mean(hgt.z.values[i, index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.mean(hgt.z.values[i, index_normal, :, :], axis=0)
	var_ninia_all = np.mean(hgt.z.values[i, index_ninia_all.values, :, :], axis=0)	
	tit = 'Composites S4 Z* 200hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'z200_composites_ENSO_' + month[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoV(var_ninio_all-var_normal_all, var_ninia_all- var_normal_all,
				    var_ninio_WPV - var_normal_all, var_ninia_WPV - var_normal_all,
				    var_ninio_SPV - var_normal_all, var_ninia_SPV - var_normal_all,
				    hgt.latitude, hgt.longitude, tit, filename)
for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	#hgt_s_WPV = hgt_s.sel(realiz=index_SPV_upper.values)
	#hgt_s_SPV = hgt_s.sel(realiz=index_SPV_lower.values)
	var_ninio_WPV = np.mean(hgt_s.z.values[index_ninio_WPV, :, :], axis=0)
	var_normal_WPV = np.mean(hgt_s.z.values[index_normal_WPV, :, :], axis=0)
	var_ninia_WPV = np.mean(hgt_s.z.values[index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.mean(hgt_s.z.values[index_ninio_SPV, :, :], axis=0)
	var_normal_SPV = np.mean(hgt_s.z.values[index_normal_SPV, :, :], axis=0)
	var_ninia_SPV = np.mean(hgt_s.z.values[index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.mean(hgt_s.z.values[index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.mean(hgt_s.z.values[index_normal, :, :], axis=0)
	var_ninia_all = np.mean(hgt_s.z.values[index_ninia_all.values, :, :], axis=0)	
	tit = 'Composites S4 Z* 200hPa Conditioned - SPoV - ' + seas[i]
	filename = FIG_PATH + 'z200_composites_ENSO_' + seas[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoV(var_ninio_all-var_normal_all, var_ninia_all- var_normal_all,
				    var_ninio_WPV - var_normal_all, var_ninia_WPV - var_normal_all,
				    var_ninio_SPV - var_normal_all, var_ninia_SPV - var_normal_all,
				    hgt.latitude, hgt.longitude, tit, filename)


