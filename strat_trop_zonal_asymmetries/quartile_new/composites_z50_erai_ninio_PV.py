#composites of EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import calendar
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT = 'hgt_erai_50.nc4'
FILE_NINIO = 'fogt/ninio34_erai_index.nc4'
FILE_PV = 'fogt/SPV_index_erai.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT)
hgt = hgt - hgt.mean(dim='longitude')
hgt['time'] = hgt.valid_time.values
hgt = hgt.sel(**{'time':slice('1981-08-01', '2018-02-28')})
hgt = hgt.sel(time= np.logical_or(hgt.time.values >= np.datetime64("2003-08-01"),
				   hgt.time.values <= np.datetime64("2002-07-31")))

ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_mon >= PV_index.SPV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_mon <= PV_index.SPV_mon.quantile(0.25, dim='dim_0', interpolation='linear')
#search for years with normal PV
index_SPV_normal = np.logical_and(PV_index.SPV_mon > PV_index.SPV_mon.quantile(0.25, dim='dim_0', interpolation='linear'), PV_index.SPV_mon < PV_index.SPV_mon.quantile(0.75, dim='dim_0', interpolation='linear'))

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

index_normal = np.logical_and(index_normal_all.values, index_SPV_normal.values)

print(sum(index_ninio_all), sum(index_ninia_all))
print(sum(index_ninio_WPV), sum(index_ninia_WPV))
print(sum(index_ninio_SPV), sum(index_ninia_SPV))
print(sum(index_normal))

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	hgt_m = hgt.sel(time= hgt['time.month'] == list(calendar.month_abbr).index(month[i]))
	var_ninio_WPV = np.mean(hgt_m.z.values[index_ninio_WPV, :, :], axis=0)
	var_ninia_WPV = np.mean(hgt_m.z.values[index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.mean(hgt_m.z.values[index_ninio_SPV, :, :], axis=0)
	var_ninia_SPV = np.mean(hgt_m.z.values[index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.mean(hgt_m.z.values[index_ninio_all.values, :, :], axis=0)
	var_normal = np.mean(hgt_m.z.values[ :, :, :], axis=0)
	var_ninia_all = np.mean(hgt_m.z.values[index_ninia_all.values, :, :], axis=0)	
	tit = 'Composites ERAI Z* 50hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'z50_erai_composites_ENSO_' + month[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoV(var_ninio_all - var_normal, var_ninia_all - var_normal,
				    var_ninio_WPV - var_normal, var_ninia_WPV - var_normal,
				    var_ninio_SPV - var_normal, var_ninia_SPV - var_normal,
				    hgt.latitude, hgt.longitude, tit, filename)
for i in np.arange(0, 5):
	hgt_erai = hgt.resample(time='QS-' + month[i]).mean(dim='time',skipna=True)
	mes = list(calendar.month_abbr).index(month[i])
	hgt_s = hgt_erai.sel(time= np.logical_and(hgt_erai['time.month'] == mes, hgt_erai['time.year']!=2002))
	var_ninio_WPV = np.mean(hgt_s.z.values[index_ninio_WPV, :, :], axis=0)
	var_ninia_WPV = np.mean(hgt_s.z.values[index_ninia_WPV, :, :], axis=0)	
	var_ninio_SPV = np.mean(hgt_s.z.values[index_ninio_SPV, :, :], axis=0)
	var_ninia_SPV = np.mean(hgt_s.z.values[index_ninia_SPV, :, :], axis=0)	
	var_ninio_all = np.mean(hgt_s.z.values[index_ninio_all.values, :, :], axis=0)
	var_normal = np.mean(hgt_s.z.values[:, :, :], axis=0)
	var_ninia_all = np.mean(hgt_s.z.values[index_ninia_all.values, :, :], axis=0)	
	tit = 'Composites ERAI Z* 50hPa Conditioned - SPoV - ' + seas[i]
	filename = FIG_PATH + 'z50_erai_composites_ENSO_' + seas[i] +'_SPoV.png'
	plots.PlotEnsoCompositesPoV(var_ninio_all - var_normal, var_ninia_all- var_normal,
				    var_ninio_WPV - var_normal, var_ninia_WPV - var_normal,
				    var_ninio_SPV - var_normal, var_ninia_SPV - var_normal,
				    hgt.latitude, hgt.longitude, tit, filename)


