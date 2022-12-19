#composites of Ks at 200hPa for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import metpy
from metpy import calc
import plots
import sys

NAME = sys.argv[1]
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_WINDS_S4 = 'fogt/monthly_winds200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
winds = xr.open_dataset(PATH_DATA + FILE_WINDS_S4, chunks={'realiz': -1, 'month': -1,
							   'latitude': 25, 'longitude': 25})
winds = winds['u']
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)
[nmonth, nrealiz, nlats, nlons] = winds.shape
#defino beta
BETA = np.transpose(np.tile(2 * 7.27e-5 * np.cos(winds.latitude.values * np.pi / 180),
			    [nmonth, nrealiz, nlons, 1]), [0, 1, 3, 2])
BETA /= 6.378e6
N2H2_term = 4 * 4e-4 * 6.4e7

F = np.transpose(np.tile(2 * 7.27e-5 * np.sin(winds.latitude.values * np.pi / 180),
			    [nmonth, nrealiz, nlons, 1]), [0, 1, 3, 2])

Rt = np.transpose(np.tile(6.378e6 * np.cos(winds.latitude.values * np.pi / 180),
			  [nmonth, nrealiz, nlons, 1]), [0, 1, 3, 2])
Umean = winds.mean(dim='latitude')
Umean = np.transpose(np.tile(Umean.values, [nlats, 1, 1, 1]), [1, 2, 0, 3])
Uyy = metpy.calc.second_derivative(winds, axis=2,
				   x=6.378e6 * np.cos(winds.latitude.values * np.pi /180))
if NAME == "Ks":
	Ks = np.sqrt((BETA - Uyy) / winds.values) * Rt
else:
	Ks = np.sqrt(BETA / winds.values - F ** 2 / N2H2_term) * Rt

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0',
								    interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0',
								    interpolation='linear')

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0',
									  interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0',
									  interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index <\
				  ninio34.ninio34_index.quantile(0.75, dim='dim_0',
								 interpolation='linear'),
				  ninio34.ninio34_index >\
				  ninio34.ninio34_index.quantile(0.25, dim='dim_0',
 interpolation='linear'))
#enso during weak PoV

index_ninio_WPV = np.logical_and(index_ninio_all.values, index_SPV_upper.values)
index_ninia_WPV = np.logical_and(index_ninia_all.values, index_SPV_upper.values)
index_normal_WPV = np.logical_and(index_normal_all.values, index_SPV_upper.values)

#enso during strong PoV
index_ninio_SPV = np.logical_and(index_ninio_all.values, index_SPV_lower.values)
index_ninia_SPV = np.logical_and(index_ninia_all.values, index_SPV_lower.values)
index_normal_SPV = np.logical_and(index_normal_all.values, index_SPV_lower.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_ninio_WPV = np.nanmean(Ks[i, index_ninio_WPV, :, :], axis=0)
	var_normal_WPV = np.nanmean(Ks[i, index_normal_WPV, :, :], axis=0)
	var_ninia_WPV = np.nanmean(Ks[i, index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.nanmean(Ks[i, index_ninio_SPV, :, :], axis=0)
	var_normal_SPV = np.nanmean(Ks[i, index_normal_SPV, :, :], axis=0)
	var_ninia_SPV = np.nanmean(Ks[i, index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.nanmean(Ks[i, index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.nanmean(Ks[i, index_normal_all.values, :, :], axis=0)
	var_ninia_all = np.nanmean(Ks[i, index_ninia_all.values, :, :], axis=0)	
	tit = 'Composites S4 Ks 200hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + NAME + '200_composites_ENSO_' + month[i] +'_SPoV.png'
	plots.PlotKsEnsoCompositesPoV(var_ninio_all-var_normal_all, var_ninia_all- var_normal_all,
				    var_ninio_WPV - var_normal_WPV, var_ninia_WPV - var_normal_WPV,
				    var_ninio_SPV - var_normal_SPV, var_ninia_SPV - var_normal_SPV,
				    winds.latitude, winds.longitude, tit, filename)
for i in np.arange(1, 5):
	winds_s = winds.isel(month=range(i, i+3)).mean(dim='month')
	Umean = np.transpose(np.tile(winds_s.mean(dim='latitude'), [nlats, 1, 1]),
			     [1, 0, 2])
	Uyy = metpy.calc.second_derivative(winds_s, axis=1,
					   x=6.378e6 * np.cos(winds.latitude.values * np.pi / 180))
	if NAME == "Ks":
		Ks = np.sqrt((BETA[0, :, :, :] - Uyy) / winds_s.values) * Rt[0, :, :, :]
	else:
		Ks = np.sqrt(BETA[0, :, :, :] / winds_s.values - F[0, :, :, :] ** 2 / N2H2_term) * Rt[0, :, :, :]
	var_ninio_WPV = np.nanmean(Ks[index_ninio_WPV, :, :], axis=0)
	var_normal_WPV = np.nanmean(Ks[index_normal_WPV, :, :], axis=0)
	var_ninia_WPV = np.nanmean(Ks[index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.nanmean(Ks[index_ninio_SPV, :, :], axis=0)
	var_normal_SPV = np.nanmean(Ks[index_normal_SPV, :, :], axis=0)
	var_ninia_SPV = np.nanmean(Ks[index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.nanmean(Ks[index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.nanmean(Ks[index_normal_all.values, :, :], axis=0)
	var_ninia_all = np.nanmean(Ks[index_ninia_all.values, :, :], axis=0)	
	tit = 'Composites S4 Ks 200hPa Conditioned - SPoV - ' + seas[i]
	filename = FIG_PATH + NAME + '200_composites_ENSO_' + seas[i] +'_SPoV.png'
	plots.PlotKsEnsoCompositesPoV(var_ninio_all-var_normal_all, var_ninia_all- var_normal_all,
				    var_ninio_WPV - var_normal_WPV, var_ninia_WPV - var_normal_WPV,
				    var_ninio_SPV - var_normal_SPV, var_ninia_SPV - var_normal_SPV,
				    winds.latitude, winds.longitude, tit, filename)


