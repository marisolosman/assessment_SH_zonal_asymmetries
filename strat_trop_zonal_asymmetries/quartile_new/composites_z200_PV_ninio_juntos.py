#composites of PV events conditioned on ENSO phase
import numpy as np
import xarray as xr
import os
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/datos/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for EN years 
index_EN = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for LN years
index_LN = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

# compute SPoV composites conditioned on ENSO phase
PV_index_EN = PV_index.SPV_index.sel(dim_0 = index_EN.values)
PV_index_LN = PV_index.SPV_index.sel(dim_0 = index_LN.values)

#PV during all phases
index_SPV_all = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_all = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#PV  during EN
index_SPV_EN = np.logical_and(index_SPV_all.values, index_EN.values)
index_WPV_EN = np.logical_and(index_WPV_all.values, index_EN.values)
#PV  during LN
index_SPV_LN =  np.logical_and(index_SPV_all.values, index_LN.values)
index_WPV_LN = np.logical_and(index_WPV_all.values, index_LN.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_WPV_EN = np.mean(hgt.z.values[i, index_WPV_EN, :, :], axis=0)
	var_WPV_LN = np.mean(hgt.z.values[i, index_WPV_LN, :, :], axis=0)	
	var_SPV_EN = np.mean(hgt.z.values[i, index_SPV_EN, :, :], axis=0)
	var_SPV_LN = np.mean(hgt.z.values[i, index_SPV_LN, :, :], axis=0)	
	var_WPV_all = np.mean(hgt.z.values[i, index_WPV_all.values, :, :], axis=0)
	var_SPV_all = np.mean(hgt.z.values[i, index_SPV_all.values, :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa Weak SPoV - Strong SPoV Conditioned - ENSO - ' + month[i]
	filename = FIG_PATH + 'z200_composites_SPoV_' + month[i] +'_ENSO_juntos.png'
	plots.PlotPoVCompositesDiffENSO(var_WPV_all-var_SPV_all, var_WPV_EN - var_SPV_EN,
				    var_WPV_LN - var_SPV_LN, hgt.latitude, hgt.longitude,
				    tit, filename)
#for i in np.arange(0, 5):
#	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
#	hgt_s_EN = hgt_s.sel(realiz=index_EN.values)
#	hgt_s_LN = hgt_s.sel(realiz=index_LN.values)
#	var_WPV_EN = np.mean(hgt_s_EN.z.values[index_WPV_EN.values, :, :], axis=0)
#	var_WPV_LN = np.mean(hgt_s_LN.z.values[index_WPV_LN.values, :, :], axis=0)	
#	var_SPV_EN = np.mean(hgt_s_EN.z.values[index_SPV_EN.values, :, :], axis=0)
#	var_SPV_LN = np.mean(hgt_s_LN.z.values[index_SPV_LN.values, :, :], axis=0)	
#	var_WPV_all = np.mean(hgt_s.z.values[index_WPV_all.values, :, :], axis=0)
#	var_SPV_all = np.mean(hgt_s.z.values[index_SPV_all.values, :, :], axis=0)
#	tit = 'Composites S4 Z* 50hPa Strong SPoV - Weak SPoV Conditioned - ENSO - ' + seas[i]
#	filename = FIG_PATH + 'z50_composites_SPoV_' + seas[i] +'_ENSO_q.png'
#	plots.PlotPoVCompositesDiffENSO(var_SPV_all-var_WPV_all, var_SPV_EN - var_WPV_EN,
#				    var_SPV_LN - var_WPV_LN, hgt.latitude, hgt.longitude,
#				    tit, filename)
#
#
