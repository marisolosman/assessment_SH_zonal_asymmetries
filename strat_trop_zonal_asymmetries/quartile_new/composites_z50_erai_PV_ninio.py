#composites of PV events conditioned on ENSO phase
import numpy as np
import xarray as xr
import os
import calendar
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT = 'hgt_erai_50.nc4'
FILE_NINIO = 'fogt/ninio34_erai_index.nc4'
FILE_PV = 'fogt/SPV_index_erai.nc4'
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV)

#search for EN years 
index_EN = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for LN years
index_LN = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

index_neutral = np.logical_and(ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'), ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'))

#PV during all phases
index_SPV_all = PV_index.SPV_mon <= PV_index.SPV_mon.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_all = PV_index.SPV_mon >= PV_index.SPV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(PV_index.SPV_mon < PV_index.SPV_mon.quantile(0.75, dim='dim_0', interpolation='linear'), PV_index.SPV_mon > PV_index.SPV_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

#PV  during EN
index_SPV_EN = np.logical_and(index_SPV_all.values, index_EN.values)
index_WPV_EN = np.logical_and(index_WPV_all.values, index_EN.values)

#PV  during EN
index_SPV_LN = np.logical_and(index_SPV_all.values, index_LN.values)
index_WPV_LN = np.logical_and(index_WPV_all.values, index_LN.values)

index_normal = np.logical_and(index_normal_all.values, index_neutral.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

hgt = xr.open_dataset(PATH_DATA + FILE_HGT)
hgt = hgt - hgt.mean(dim='longitude')
hgt['time'] = hgt.valid_time.values
hgt = hgt.sel(**{'time':slice('1981-08-01', '2018-02-28')})
hgt = hgt.sel(time=np.logical_or(hgt.time.values >= np.datetime64('2003-08-01'),
				 hgt.time.values <= np.datetime64('2002-07-31')))

for i in np.arange(0, 7):
	hgt_m = hgt.sel(time= hgt['time.month'] == list(calendar.month_abbr).index(month[i]))
	var_WPV_EN = np.mean(hgt_m.z.values[index_WPV_EN, :, :], axis=0)
	var_WPV_LN = np.mean(hgt_m.z.values[index_WPV_LN, :, :], axis=0)
	var_normal = np.mean(hgt_m.z.values[:, :, :], axis=0)
	var_SPV_EN = np.mean(hgt_m.z.values[index_SPV_EN, :, :], axis=0)
	var_SPV_LN = np.mean(hgt_m.z.values[index_SPV_LN, :, :], axis=0)	
	var_WPV_all = np.mean(hgt_m.z.values[index_WPV_all.values, :, :], axis=0)
	var_SPV_all = np.mean(hgt_m.z.values[index_SPV_all.values, :, :], axis=0)
	tit = 'Composites ERAI Z* 50hPa Conditioned - ENSO - ' + month[i]
	filename = FIG_PATH + 'z50_erai_composites_SPoV2_' + month[i] +'_ENSO.png'
	plots.PlotPoVCompositesENSO(var_SPV_all-var_normal, var_WPV_all-var_normal,
				    var_SPV_EN - var_normal, var_WPV_EN - var_normal,
				    var_SPV_LN - var_normal, var_WPV_LN - var_normal,
				    hgt.latitude, hgt.longitude,
				    tit, filename)
for i in np.arange(0, 5):
	hgt_erai = hgt.resample(time='QS-' + month[i]).mean(dim='time', skipna=True)
	mes = list(calendar.month_abbr).index(month[i])
	hgt_s = hgt_erai.sel(time=np.logical_and(hgt_erai['time.month'] ==mes,
						 hgt_erai['time.year']!=2002))
	var_WPV_EN = np.mean(hgt_s.z.values[index_WPV_EN, :, :], axis=0)
	var_WPV_LN = np.mean(hgt_s.z.values[index_WPV_LN, :, :], axis=0)	
	var_SPV_EN = np.mean(hgt_s.z.values[index_SPV_EN, :, :], axis=0)
	var_SPV_LN = np.mean(hgt_s.z.values[index_SPV_LN, :, :], axis=0)	
	var_WPV_all = np.mean(hgt_s.z.values[index_WPV_all.values, :, :], axis=0)
	var_SPV_all = np.mean(hgt_s.z.values[index_SPV_all.values, :, :], axis=0)
	var_normal = np.mean(hgt_s.z.values[:, :, :], axis=0)
	tit = 'Composites ERAI Z* 50hPa Conditioned - ENSO - ' + seas[i]
	filename = FIG_PATH + 'z50_erai_composites_SPoV2_' + seas[i] +'_ENSO_new.png'
	plots.PlotPoVCompositesENSO(var_SPV_all - var_normal, var_WPV_all - var_normal,
				    var_SPV_EN - var_normal, var_WPV_EN - var_normal,
				    var_SPV_LN - var_normal, var_WPV_LN - var_normal,
				    hgt.latitude, hgt.longitude,
				    tit, filename)
