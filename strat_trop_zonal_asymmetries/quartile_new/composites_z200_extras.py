#composites of remaining combinations
import numpy as np
import xarray as xr
import os
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
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
#search for years with normal PV
index_SPV_normal = np.logical_and(PV_index.SPV_index > PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'), PV_index.SPV_index < PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'))

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))
#enso during weak PoV

index_ninio_normal = np.logical_and(index_ninio_all.values, index_SPV_normal.values)
index_ninio_SPVPWPV = np.logical_and(index_ninio_all.values,
				     np.logical_or(index_SPV_upper.values, index_SPV_lower.values))

index_ninia_normal = np.logical_and(index_ninia_all.values, index_SPV_normal.values)
index_ninia_SPVPWPV = np.logical_and(index_ninia_all.values,
				     np.logical_or(index_SPV_upper.values, index_SPV_lower.values))

hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_ninio_normal = np.mean(hgt.z.values[i, index_ninio_normal, :, :], axis=0)
	var_ninio_SPVPWPV = np.mean(hgt.z.values[i, index_ninio_SPVPWPV, :, :], axis=0)	
	var_normal = np.mean(hgt.z.values[i, :, :, :], axis=0)
	var_ninia_normal = np.mean(hgt.z.values[i, index_ninia_normal, :, :], axis=0)
	var_ninia_SPVPWPV = np.mean(hgt.z.values[i, index_ninia_SPVPWPV, :, :], axis=0)	
	tit = 'Composites S4 Z* 200hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'z200_composites_ENSO_' + month[i] +'_SPoV_extras.png'
	plots.PlotEnsoCompositesPoVex(var_ninio_normal - var_normal, var_ninia_normal - var_normal,
				    var_ninio_SPVPWPV - var_normal, var_ninia_SPVPWPV - var_normal,
				    hgt.latitude, hgt.longitude, tit, filename)
for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	var_ninio_normal = np.mean(hgt_s.z.values[index_ninio_normal, :, :], axis=0)
	var_ninio_SPVPWPV = np.mean(hgt_s.z.values[index_ninio_SPVPWPV, :, :], axis=0)	
	var_normal = np.mean(hgt_s.z.values[ :, :, :], axis=0)
	var_ninia_normal = np.mean(hgt_s.z.values[index_ninia_normal, :, :], axis=0)
	var_ninia_SPVPWPV = np.mean(hgt_s.z.values[index_ninia_SPVPWPV, :, :], axis=0)	
	tit = 'Composites S4 Z* 200hPa Conditioned - SPoV - ' + seas[i]
	filename = FIG_PATH + 'z200_composites_ENSO_' + seas[i] +'_SPoV_extras.png'
	plots.PlotEnsoCompositesPoVex(var_ninio_normal - var_normal, var_ninia_normal - var_normal,
				    var_ninio_SPVPWPV - var_normal, var_ninia_SPVPWPV - var_normal,
				    hgt.latitude, hgt.longitude, tit, filename)

index_SPV_normal = np.logical_and(index_SPV_lower.values, index_normal_all.values)
index_SPV_ENPLN = np.logical_and(index_SPV_lower.values,
				     np.logical_or(index_ninio_all.values, index_ninia_all.values))
index_WPV_normal = np.logical_and(index_SPV_upper.values, index_normal_all.values)
index_WPV_ENPLN = np.logical_and(index_SPV_upper.values,
				     np.logical_or(index_ninio_all.values, index_ninia_all.values))
for i in np.arange(0, 7):
	var_SPV_normal = np.mean(hgt.z.values[i, index_SPV_normal, :, :], axis=0)
	var_SPV_ENPLN = np.mean(hgt.z.values[i, index_SPV_ENPLN, :, :], axis=0)	
	var_WPV_normal = np.mean(hgt.z.values[i, index_WPV_normal, :, :], axis=0)
	var_WPV_ENPLN = np.mean(hgt.z.values[i, index_WPV_ENPLN, :, :], axis=0)	
	var_normal = np.mean(hgt.z.values[i, :, :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa Conditioned - ENSO - ' + month[i]
	filename = FIG_PATH + 'z200_composites_SPoV2_' + month[i] +'_ENSO_extras.png'
	plots.PlotPoVCompositesENSOex(var_SPV_normal - var_normal, var_WPV_normal - var_normal,
				    var_SPV_ENPLN - var_normal, var_WPV_ENPLN - var_normal,
				    hgt.latitude, hgt.longitude, tit, filename)
for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	var_SPV_normal = np.mean(hgt_s.z.values[index_SPV_normal, :, :], axis=0)
	var_SPV_ENPLN = np.mean(hgt_s.z.values[index_SPV_ENPLN, :, :], axis=0)	
	var_WPV_normal = np.mean(hgt_s.z.values[index_WPV_normal, :, :], axis=0)
	var_WPV_ENPLN = np.mean(hgt_s.z.values[index_WPV_ENPLN, :, :], axis=0)	
	var_normal = np.mean(hgt_s.z.values[:, :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa Conditioned - ENSO - ' + seas[i]
	filename = FIG_PATH + 'z200_composites_SPoV2_' + seas[i] +'_"NSO_extras.png'
	plots.PlotPoVCompositesENSOex(var_SPV_normal - var_normal, var_WPV_normal - var_normal,
				    var_SPV_ENPLN - var_normal, var_WPV_ENPLN - var_normal,
				    hgt.latitude, hgt.longitude, tit, filename)

