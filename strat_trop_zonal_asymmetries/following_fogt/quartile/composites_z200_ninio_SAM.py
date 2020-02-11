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
FILE_SAM_S4 = 'fogt/SAM_monthly_index_s4.nc4'

hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
SAM_index =  xr.open_dataset(PATH_DATA + FILE_SAM_S4)

#search for years with negative SAM
index_SAM_negative = SAM_index.SAM_index <= SAM_index.SAM_index.quantile(0.25, dim='realiz', interpolation='linear')
#search for years with positive SAM
index_SAM_positive = SAM_index.SAM_index >= SAM_index.SAM_index.quantile(0.75, dim='realiz', interpolation='linear')
#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']

for i in np.arange(0, 7):

	# compute EN-LA composites conditioned on SAM anomalies
	ninio34_WSAM = ninio34.sel(dim_0 = index_SAM_negative.values[i, :])
	hgt_WSAM = hgt.sel(realiz = index_SAM_negative.values[i, :])
	ninio34_SSAM = ninio34.sel(dim_0 = index_SAM_positive.values[i, :])
	hgt_SSAM = hgt.sel(realiz = index_SAM_positive.values[i, :])

	#enso during weak SAM
	index_ninio_WSAM = ninio34_WSAM.ninio34_index >= ninio34_WSAM.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
	index_ninia_WSAM = ninio34_WSAM.ninio34_index <= ninio34_WSAM.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
	index_normal_WSAM = np.logical_and(ninio34_WSAM.ninio34_index < ninio34_WSAM.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34_WSAM.ninio34_index > ninio34_WSAM.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))
	#enso during strong SAM
	index_ninio_SSAM = ninio34_SSAM.ninio34_index >= ninio34_SSAM.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
	index_ninia_SSAM = ninio34_SSAM.ninio34_index <= ninio34_SSAM.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
	index_normal_SSAM = np.logical_and(ninio34_SSAM.ninio34_index < ninio34_SSAM.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34_SSAM.ninio34_index > ninio34_SSAM.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))
	var_ninio_WSAM = np.mean(hgt_WSAM.z.values[i, index_ninio_WSAM.values, :, :], axis=0)
	var_normal_WSAM = np.mean(hgt_WSAM.z.values[i, index_normal_WSAM.values, :, :], axis=0)
	var_ninia_WSAM = np.mean(hgt_WSAM.z.values[i, index_ninia_WSAM.values, :, :], axis=0)	
	var_ninio_SSAM = np.mean(hgt_SSAM.z.values[i, index_ninio_SSAM.values, :, :], axis=0)
	var_normal_SSAM = np.mean(hgt_SSAM.z.values[i, index_normal_SSAM.values, :, :], axis=0)
	var_ninia_SSAM = np.mean(hgt_SSAM.z.values[i, index_ninia_SSAM.values, :, :], axis=0)	
	var_ninio_all = np.mean(hgt.z.values[i, index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.mean(hgt.z.values[i, index_normal_all.values, :, :], axis=0)
	var_ninia_all = np.mean(hgt.z.values[i, index_ninia_all.values, :, :], axis=0)	
	tit = 'Composites S4 Z* 200hPa Conditioned - SAM - ' + month[i]
	filename = FIG_PATH + 'z200_composites_ENSO_' + month[i] +'_SAM_q.png'
	plots.PlotEnsoCompositesSAM(var_ninio_all-var_normal_all, var_ninia_all- var_normal_all,
				    var_ninio_WSAM - var_normal_WSAM, var_ninia_WSAM - var_normal_WSAM,
				    var_ninio_SSAM - var_normal_SSAM, var_ninia_WSAM - var_normal_SSAM,
				    hgt.latitude, hgt.longitude, tit, filename)

FILE_SAM_S4 = 'fogt/SAM_index_s4.nc4'

SAM_index =  xr.open_dataset(PATH_DATA + FILE_SAM_S4)

#search for years with negative SAM
index_SAM_negative = SAM_index.SAM_index <= SAM_index.SAM_index.quantile(0.25, dim='realiz', interpolation='linear')
#search for years with positive SAM
index_SAM_positive = SAM_index.SAM_index >= SAM_index.SAM_index.quantile(0.75, dim='realiz', interpolation='linear')
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 5):
	# compute EN-LA composites conditioned on SAM anomalies
	ninio34_WSAM = ninio34.sel(dim_0 = index_SAM_negative.values[i, :])
	ninio34_SSAM = ninio34.sel(dim_0 = index_SAM_positive.values[i, :])

	#enso during weak SAM
	index_ninio_WSAM = ninio34_WSAM.ninio34_index >= ninio34_WSAM.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
	index_ninia_WSAM = ninio34_WSAM.ninio34_index <= ninio34_WSAM.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
	index_normal_WSAM = np.logical_and(ninio34_WSAM.ninio34_index < ninio34_WSAM.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34_WSAM.ninio34_index > ninio34_WSAM.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))
	#enso during strong SAM
	index_ninio_SSAM = ninio34_SSAM.ninio34_index >= ninio34_SSAM.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
	index_ninia_SSAM = ninio34_SSAM.ninio34_index <= ninio34_SSAM.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
	index_normal_SSAM = np.logical_and(ninio34_SSAM.ninio34_index < ninio34_SSAM.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34_SSAM.ninio34_index > ninio34_SSAM.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))

	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	hgt_s_WSAM = hgt_s.sel(realiz=index_SAM_negative.values[i, :])
	hgt_s_SSAM = hgt_s.sel(realiz=index_SAM_positive.values[i, :])
	var_ninio_WSAM = np.mean(hgt_s_WSAM.z.values[index_ninio_WSAM.values, :, :], axis=0)
	var_normal_WSAM = np.mean(hgt_s_WSAM.z.values[index_normal_WSAM.values, :, :], axis=0)
	var_ninia_WSAM = np.mean(hgt_s_WSAM.z.values[index_ninia_WSAM.values, :, :], axis=0)	
	var_ninio_SSAM = np.mean(hgt_s_SSAM.z.values[index_ninio_SSAM.values, :, :], axis=0)
	var_normal_SSAM = np.mean(hgt_s_SSAM.z.values[index_normal_SSAM.values, :, :], axis=0)
	var_ninia_SSAM = np.mean(hgt_s_SSAM.z.values[index_ninia_SSAM.values, :, :], axis=0)	
	var_ninio_all = np.mean(hgt_s.z.values[index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.mean(hgt_s.z.values[index_normal_all.values, :, :], axis=0)
	var_ninia_all = np.mean(hgt_s.z.values[index_ninia_all.values, :, :], axis=0)	
	tit = 'Composites S4 Z* 200hPa Conditioned - SAM - ' + seas[i]
	filename = FIG_PATH + 'z200_composites_ENSO_' + seas[i] +'_SAM_q.png'
	plots.PlotEnsoCompositesSAM(var_ninio_all-var_normal_all, var_ninia_all- var_normal_all,
				    var_ninio_WSAM - var_normal_WSAM, var_ninia_WSAM - var_normal_WSAM,
				    var_ninio_SSAM - var_normal_SSAM, var_ninia_WSAM - var_normal_SSAM,
				    hgt.latitude, hgt.longitude, tit, filename)


