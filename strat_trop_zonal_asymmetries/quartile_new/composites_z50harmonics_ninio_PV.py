#composites of z200 and plumb fluxes for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots_stereo as plots
import random
import doharmonic

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'
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

harms = xr.open_dataset(PATH_DATA_2 + 'hgt50_monthly_harmonics.nc4')
harm1 = harms.Harmonic1.values
harm2 = harms.Harmonic2.values

for i in np.arange(0, 7):
	var_normal = np.mean(harm1[i, :, :, :], axis=0)
	var_ninio_WPV = np.mean(harm1[i, index_ninio_WPV, :, :], axis=0)
	var_ninia_WPV = np.mean(harm1[i, index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.mean(harm1[i, index_ninio_SPV, :, :], axis=0)
	var_ninia_SPV = np.mean(harm1[i, index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.mean(harm1[i, index_ninio_all, :, :], axis=0)
	var_ninia_all = np.mean(harm1[i, index_ninia_all, :, :], axis=0)	
	tit = 'Composites S4 Z* 50hPa Wave-1 Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'z50W1_composites_ENSO_' + month[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoV(var_ninio_all-var_normal, var_ninia_all-var_normal,
				    var_ninio_WPV-var_normal, var_ninia_WPV-var_normal,
				    var_ninio_SPV-var_normal, var_ninia_SPV-var_normal,
				    var_normal, harms.latitude, harms.longitude, tit, filename)
	var_normal = np.mean(harm2[i, :, :, :], axis=0)
	var_ninio_WPV = np.mean(harm2[i, index_ninio_WPV, :, :], axis=0)
	var_ninia_WPV = np.mean(harm2[i, index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.mean(harm2[i, index_ninio_SPV, :, :], axis=0)
	var_ninia_SPV = np.mean(harm2[i, index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.mean(harm2[i, index_ninio_all, :, :], axis=0)
	var_ninia_all = np.mean(harm2[i, index_ninia_all, :, :], axis=0)	
	tit = 'Composites S4 Z* 50hPa Wave-2 Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'z50W2_composites_ENSO_' + month[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoV(var_ninio_all-var_normal, var_ninia_all-var_normal,
				    var_ninio_WPV-var_normal, var_ninia_WPV-var_normal,
				    var_ninio_SPV-var_normal, var_ninia_SPV-var_normal,
				    var_normal, harms.latitude, harms.longitude, tit, filename)

harms = xr.open_dataset(PATH_DATA_2 + 'hgt50_seasonal_harmonics.nc4')
harm1 = harms.Harmonic1.values
harm2 = harms.Harmonic2.values
for i in np.arange(0, 5):
	var_normal = np.mean(harm1[i, :, :, :], axis=0)
	var_ninio_WPV = np.mean(harm1[i, index_ninio_WPV, :, :], axis=0)
	var_ninia_WPV = np.mean(harm1[i, index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.mean(harm1[i, index_ninio_SPV, :, :], axis=0)
	var_ninia_SPV = np.mean(harm1[i, index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.mean(harm1[i, index_ninio_all, :, :], axis=0)
	var_ninia_all = np.mean(harm1[i, index_ninia_all, :, :], axis=0)	
	tit = 'Composites S4 Z* 50hPa Wave-1 Conditioned - SPoV - ' + seas[i]
	filename = FIG_PATH + 'z50W1_composites_ENSO_' + seas[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoV(var_ninio_all-var_normal, var_ninia_all-var_normal,
				    var_ninio_WPV-var_normal, var_ninia_WPV-var_normal,
				    var_ninio_SPV-var_normal, var_ninia_SPV-var_normal,
				    var_normal, harms.latitude, harms.longitude, tit, filename)
	var_normal = np.mean(harm2[i, :, :, :], axis=0)
	var_ninio_WPV = np.mean(harm2[i, index_ninio_WPV, :, :], axis=0)
	var_ninia_WPV = np.mean(harm2[i, index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.mean(harm2[i, index_ninio_SPV, :, :], axis=0)
	var_ninia_SPV = np.mean(harm2[i, index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.mean(harm2[i, index_ninio_all, :, :], axis=0)
	var_ninia_all = np.mean(harm2[i, index_ninia_all, :, :], axis=0)	
	tit = 'Composites S4 Z* 50hPa Wave-2 Conditioned - SPoV - ' + seas[i]
	filename = FIG_PATH + 'z50W2_composites_ENSO_' + seas[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoV(var_ninio_all-var_normal, var_ninia_all-var_normal,
				    var_ninio_WPV-var_normal, var_ninia_WPV-var_normal,
				    var_ninio_SPV-var_normal, var_ninia_SPV-var_normal,
				    var_normal, harms.latitude, harms.longitude, tit, filename)
		
