#plot composites enso phases  and pv intensity for ERAI and S4
#add significance test
import numpy as np
import datetime
import pandas as pd
import xarray as xr
from scipy import signal
from scipy import stats
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
tvalue = [stats.t.ppf(1-0.025, i) for i in tvalue]

def PlotComposites(model_m, lat_model, lon_model, observ_m, lat_observ,
		    lon_observ,	title, filename):
	proj = ccrs.Stereographic(central_longitude=-60, centra_latitude=-90)
	fig = plt.figure(1,(10, 6.7), 300)
	ax =plt.subplot(2, 1, 1, projection=proj)
	clevs = np.arange(-60, 70, 10)
	barra = plt.cm.RdBu_r
	ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
	gl = ax.gridlines(draw_labels=False, xlocs=None, ylocs=None)
	gl.n_steps = 90
	t_value = stats.t.ppf(1 - 0.025, observ_m['df'])
	M = np.ma.array(observ_m['var'], mask=np.abs(observ_m['mask']) < t_value)
	im = ax.contourf(lon_observ, lat_observ, M, clevs,
			 transform=ccrs.PlateCarree(), cmap=barra, extend='both', vmin=-60, vmax=70)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	ax.coastlines()
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	plt.title('ERA Interim', fontsize=10)
	ax =plt.subplot(2, 1, 2, projection=proj)
	ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
	gl = ax.gridlines(draw_labels=False, xlocs=None, ylocs=None)
	gl.n_steps = 90
	t_value = stats.t.ppf(1 - 0.025, model_m['df'])
	M = np.ma.array(model_m['var'], mask=np.abs(model_m['mask']) < t_value)
	im = ax.contourf(lon_model, lat_model, M, clevs,
			 transform=ccrs.PlateCarree(), cmap=barra, extend='both',
			 vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	ax.add_feature(cartopy.feature.COASTLINE)
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	plt.title('ECMWF S4', fontsize=10)
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures_paper/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt_s4 = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt_s4 = hgt_s4 - hgt_s4.mean(dim='longitude')
FILE_HGT_ERAI = 'hgt_erai_200.nc4'
FILE_NINIO_ERAI = 'fogt/ninio34_erai_index.nc4'
FILE_PV_ERAI = 'fogt/SPV_index_erai.nc4'
hgt_erai = xr.open_dataset(PATH_DATA + FILE_HGT_ERAI)
hgt_erai.time.values = hgt_erai.valid_time.values
#remove data from august 2002 to february 2003
hgt_erai = hgt_erai.sel(time=hgt_erai.time.values[np.logical_or(hgt_erai.time.values<=np.datetime64('2002-07-31'), hgt_erai.time.values>=np.datetime64('2003-03-01'))])
hgt_erai = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')})
hgt_erai = hgt_erai - hgt_erai.mean(dim='longitude')
#abrir archivo de sst y PV del erai
ninio34_erai =  xr.open_dataset(PATH_DATA + FILE_NINIO_ERAI)
PV_erai =  xr.open_dataset(PATH_DATA + FILE_PV_ERAI)

#select ninio, ninia and neutral
index_ninio_erai_upper = ninio34_erai.ninio34_index >= ninio34_erai.ninio34_index.quantile(0.75, dim='dim_0')
index_PV_erai_upper = PV_erai.SPV_mon >= PV_erai.SPV_mon.quantile(0.75, dim='dim_0')

index_ninio_erai_lower = ninio34_erai.ninio34_index <= ninio34_erai.ninio34_index.quantile(0.25, dim='dim_0')
index_PV_erai_lower = PV_erai.SPV_mon <= PV_erai.SPV_mon.quantile(0.25, dim='dim_0')

index_ninio_erai_normal = np.logical_and(ninio34_erai.ninio34_index < ninio34_erai.ninio34_index.quantile(0.75, dim='dim_0'), ninio34_erai.ninio34_index > ninio34_erai.ninio34_index.quantile(0.25, dim='dim_0'))
index_PV_erai_normal = np.logical_and(PV_erai.SPV_mon < PV_erai.SPV_mon.quantile(0.75, dim='dim_0'), PV_erai.SPV_mon > PV_erai.SPV_mon.quantile(0.25, dim='dim_0'))

#idem s4
#abrir archivo de sst y PV del erai
ninio34_s4 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_s4 =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#select ninio, ninia and neutral
index_ninio_s4_upper = ninio34_s4.ninio34_index >= ninio34_s4.ninio34_index.quantile(0.75, dim='dim_0')
index_PV_s4_upper = PV_s4.SPV_index >= PV_s4.SPV_index.quantile(0.75, dim='dim_0')

index_ninio_s4_lower = ninio34_s4.ninio34_index <= ninio34_s4.ninio34_index.quantile(0.25, dim='dim_0')
index_PV_s4_lower = PV_s4.SPV_index <= PV_s4.SPV_index.quantile(0.25, dim='dim_0')

index_ninio_s4_normal = np.logical_and(ninio34_s4.ninio34_index < ninio34_s4.ninio34_index.quantile(0.75, dim='dim_0'), ninio34_s4.ninio34_index > ninio34_s4.ninio34_index.quantile(0.25, dim='dim_0'))
index_PV_s4_normal = np.logical_and(PV_s4.SPV_index < PV_s4.SPV_index.quantile(0.75, dim='dim_0'), PV_s4.SPV_index > PV_s4.SPV_index.quantile(0.25, dim='dim_0'))

lon_s4, lat_s4 = np.meshgrid(hgt_s4.longitude.values, hgt_s4.latitude.values)
lon_erai, lat_erai = np.meshgrid(hgt_erai.longitude.values, hgt_erai.latitude.values)

#seasonal means
season = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lmonth = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec']
#loop over seasons, selec data and performs all composites
for i in np.arange(0, 5):
	hgt_erai_seas_mean = hgt_erai['z'].resample(time='QS-' + lmonth[i]).mean(dim='time',skipna=True)
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	hgt_erai_smean = hgt_erai_seas_mean.sel(time= np.logical_and(hgt_erai_seas_mean['time.month'] == mes, hgt_erai_seas_mean['time.year']!=2002))
	hgt_erai_smean= np.nan_to_num(hgt_erai_smean, np.nanmean(hgt_erai_smean))
	hgt_erai_smean = signal.detrend(hgt_erai_smean, axis=0, type='linear')
	hgt_erai_EN = np.nanmean(hgt_erai_smean[index_ninio_erai_upper.values, :, :], axis=0)
	SS_erai_EN = np.nanvar(hgt_erai_smean[index_ninio_erai_upper.values, :, :],
			       axis=0) / np.sum(index_erai_upper.values)
	hgt_erai_LN = np.nanmean(hgt_erai_smean[index_ninio_erai_lower.values, :, :], axis=0)
	SS_erai_LN = np.nanvar(hgt_erai_smean[index_ninio_erai_lower.values, :, :],
			       axis=0) / np.sum(index_erai_lower.values)
	hgt_erai_N = np.nanmean(hgt_erai_smean[index_ninio_erai_normal.values, :, :], axis=0)
	hgt_erai_WSPV = np.nanmean(hgt_erai_smean[index_PV_erai_upper.values, :, :], axis=0)
	SS_erai_WSPV = np.nanvar(hgt_erai_smean[index_PV_erai_upper.values, :, :],
			       axis=0) / np.sum(index_PV_erai_upper.values)
	hgt_erai_SSPV = np.nanmean(hgt_erai_smean[index_PV_erai_lower.values, :, :], axis=0)
	SS_erai_SSPV = np.nanvar(hgt_erai_smean[index_PV_erai_lower.values, :, :],
			       axis=0) / np.sum(index_PV_erai_lower.values)
	hgt_erai_NSPV = np.nanmean(hgt_erai_smean[index_PV_erai_normal.values, :, :], axis=0)
	hgt_s4_smean = np.nanmean(hgt_s4.z.values[i:i + 3, :, :, :], axis=0)	
	hgt_s4_smean= np.nan_to_num(hgt_s4_smean, np.nanmean(hgt_s4_smean))
	hgt_s4_smean = signal.detrend(hgt_s4_smean, axis=0, type='linear')
	hgt_s4_EN = np.nanmean(hgt_s4_smean[index_ninio_s4_upper.values, :, :], axis=0)
	SS_s4_EN = np.nanvar(hgt_s4_smean[index_ninio_s4_upper.values, :, :],
			       axis=0) / np.sum(index_s4_upper.values)
	hgt_s4_LN = np.nanmean(hgt_s4_smean[index_ninio_s4_lower.values, :, :], axis=0)
	SS_s4_LN = np.nanvar(hgt_s4_smean[index_ninio_s4_lower.values, :, :],
			       axis=0) / np.sum(index_s4_lower.values)
	hgt_s4_N = np.nanmean(hgt_s4_smean[index_ninio_s4_normal.values, :, :], axis=0)
	hgt_s4_WSPV = np.nanmean(hgt_s4_smean[index_PV_s4_upper.values, :, :], axis=0)
	SS_s4_WSPV = np.nanvar(hgt_s4_smean[index_PV_s4_upper.values, :, :],
			       axis=0) / np.sum(index_PV_s4_upper.values)
	hgt_s4_SSPV = np.nanmean(hgt_s4_smean[index_PV_s4_lower.values, :, :], axis=0)
	SS_s4_SSPV = np.nanvar(hgt_s4_smean[index_PV_s4_lower.values, :, :],
			       axis=0) / np.sum(index_PV_s4_lower.values)
	hgt_s4_NSPV = np.nanmean(hgt_s4_smean[index_PV_s4_normal.values, :, :], axis=0)
	#composites EN-N
	#title = 'Composites Z* 200hPa NINIO years - ' + season[i]
	#filename = FIG_PATH + 'Composites_NINIO_' + season[i] + '_det.png'
	#model = hgt_s4_EN - hgt_s4_N
	#obs = hgt_erai_EN - hgt_erai_N
	#PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
	#composites LN-N
	#title = 'Composites Z* 200hPa NINIA years - ' + season[i]
	#filename = FIG_PATH + 'Composites_NINIA_' + season[i] + '_det.png'
	#model = hgt_s4_LN - hgt_s4_N
	#obs = hgt_erai_LN - hgt_erai_N
	#PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
	#composites EN-LN
	title = 'Composites Z* 200hPa NINIO-NINIA years - ' + season[i]
	filename = FIG_PATH + 'Composites_ENSO_' + season[i] + '_det.png'
	model = hgt_s4_EN - hgt_s4_LN
	obs = hgt_erai_EN - hgt_erai_LN
	var_model = {'var': model, 'mask': model / np.sqrt(SS_s4_EN + SS_s4_LN),
		     'df': np.sum(index_ninio_s4_upper.values) + np.sum(index_ninio_s4_lower.values)}
	var_obs = {'var': obs, 'mask': obs / np.sqrt(SS_erai_EN + SS_erai_LN),
		     'df': np.sum(index_erai_upper.values) + np.sum(index_erai_lower.values)}
	PlotComposites(var_model, lat_s4, lon_s4, var_obs, lat_erai, lon_erai, title, filename)
	#composites SSPV
#	title = 'Composites Z* 200hPa Strong SPV years - ' + season[i]
#	filename = FIG_PATH + 'Composites_SSPV_' + season[i] + '_det.png'
#	model = hgt_s4_SSPV - hgt_s4_NSPV
#	obs = hgt_erai_SSPV - hgt_erai_NSPV
#	PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
#	#composites WSPV
#	title = 'Composites Z* 200hPa Weak SPV years - ' + season[i]
#	filename = FIG_PATH + 'Composites_WSPV_' + season[i] + '_det.png'
#	model = hgt_s4_WSPV - hgt_s4_NSPV
#	obs = hgt_erai_WSPV - hgt_erai_NSPV
#	PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
	#composites SSPV - WsPV
	title = 'Composites Z* 200hPa Strong-Weak SPV years - ' + season[i]
	filename = FIG_PATH + 'Composites_SPV_' + season[i] + '_det.png'
	model = hgt_s4_SSPV - hgt_s4_WSPV
	obs = hgt_erai_SSPV - hgt_erai_WSPV
	var_model = {'var': model, 'mask': model / np.sqrt(SS_s4_SSPV + SS_s4_WSPV),
		     'df': np.sum(index_PV_s4_upper.values) + np.sum(index_PV_s4_lower.values)}
	var_obs = {'var': obs, 'mask': obs / np.sqrt(SS_erai_SSPV + SS_erai_WSPV),
		     'df': np.sum(index_PV_erai_upper.values) + np.sum(index_PV_erai_lower.values)}
	PlotComposites(var_model, lat_s4, lon_s4, var_obs, lat_erai, lon_erai, title, filename)

