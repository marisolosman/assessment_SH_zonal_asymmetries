#composites of z50* harmonics for PV events conditioned on EN events 
import numpy as np
import xarray as xr
import os
import plots_stereo as plots
import doharmonic

#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for EN years 
index_EN = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for LN years
index_LN = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

index_neutral = np.logical_and(ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'), ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'))

#PV during all phases
index_SPV_all = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_all = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(PV_index.SPV_index < PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'), PV_index.SPV_index > PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'))

#PV  during EN
index_SPV_EN = np.logical_and(index_SPV_all.values, index_EN.values)
index_WPV_EN = np.logical_and(index_WPV_all.values, index_EN.values)

#PV  during EN
index_SPV_LN = np.logical_and(index_SPV_all.values, index_LN.values)
index_WPV_LN = np.logical_and(index_WPV_all.values, index_LN.values)

index_normal = np.logical_and(index_normal_all.values, index_neutral.values)


month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

harms = xr.open_dataset(PATH_DATA_2 + 'hgt200_monthly_harmonics.nc4')
for i in np.arange(0, 7):
	for j in np.arange(5):
		var_normal = np.mean(harms.Harmonics.values[j, i, :, :, :], axis=0)
		var_WPV_EN = np.mean(harms.Harmonics.values[j, i, index_WPV_EN, :, :], axis=0)
		var_WPV_LN = np.mean(harms.Harmonics.values[j, i, index_WPV_LN, :, :], axis=0)	
		var_SPV_EN = np.mean(harms.Harmonics.values[j, i, index_SPV_EN, :, :], axis=0)	
		var_SPV_LN = np.mean(harms.Harmonics.values[j, i, index_SPV_LN, :, :], axis=0)
		var_WPV_all = np.mean(harms.Harmonics.values[j, i, index_WPV_all, :, :], axis=0)	
		var_SPV_all = np.mean(harms.Harmonics.values[j, i, index_SPV_all, :, :], axis=0)	
		tit = 'Composites S4 Z* 200hPa Wave-' + str(j + 1) + ' Conditioned - ENSO - ' + month[i]
		filename = FIG_PATH + 'z200W' + str(j + 1) + '_composites_SPoV2_' + month[i] +'_ENSO.png'
		plots.PlotPoVCompositesENSO(var_SPV_all - var_normal, var_WPV_all - var_normal,
				    	    var_SPV_EN - var_normal, var_WPV_EN - var_normal,
				    	    var_SPV_LN - var_normal, var_WPV_LN - var_normal,
				    	    var_normal, harms.latitude, harms.longitude, tit, filename)

harms = xr.open_dataset(PATH_DATA_2 + 'hgt200_seasonal_harmonics.nc4')
for i in np.arange(0, 5):
	for j in np.arange(5):
		var_normal = np.mean(harms.Harmonic1.values[j, i, :, :, :], axis=0)
		var_WPV_EN = np.mean(harms.Harmonic1.values[j, i, index_WPV_EN, :, :], axis=0)
		var_WPV_LN = np.mean(harms.Harmonic1.values[j, i, index_WPV_LN, :, :], axis=0)	
		var_SPV_EN = np.mean(harms.Harmonic1.values[j, i, index_SPV_EN, :, :], axis=0)	
		var_SPV_LN = np.mean(harms.Harmonic1.values[j, i, index_SPV_LN, :, :], axis=0)
		var_WPV_all = np.mean(harms.Harmonic1.values[j, i, index_WPV_all, :, :], axis=0)	
		var_SPV_all = np.mean(harms.Harmonic1.values[j, i, index_SPV_all, :, :], axis=0)	
		tit = 'Composites S4 Z* 200hPa Wave-' + str(j + 1) + ' Conditioned - ENSO - ' + seas[i]
		filename = FIG_PATH + 'z200W' + str(j + 1) + '_composites_SPoV2_' + seas[i] +'_ENSO.png'
		plots.PlotPoVCompositesENSO(var_SPV_all - var_normal, var_WPV_all - var_normal,
					    var_SPV_EN - var_normal, var_WPV_EN - var_normal,
					    var_SPV_LN - var_normal, var_WPV_LN - var_normal,
					    var_normal, harms.latitude, harms.longitude, tit, filename)
		
