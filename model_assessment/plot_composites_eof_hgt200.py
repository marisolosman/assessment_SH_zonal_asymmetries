#plor meand and std for HGT 200 A-N for ERA Interim and S4
import numpy as np
import datetime
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
def PlotComposites(model_m, lat_model, lon_model, observ_m, lat_observ,
		    lon_observ,	title, filename):
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1,(10, 6.7), 300)
	ax =plt.subplot(2, 1, 1, projection=proj)
	clevs = np.arange(-60, 70, 10)
	barra = plt.cm.RdBu_r
	ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
	im = ax.contourf(lon_observ, lat_observ, observ_m, clevs,
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
	im = ax.contourf(lon_model, lat_model, model_m, clevs,
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

ROUTE = '~/datos/data/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc'
FILE_NINIO_S4 = 'ninio34_index.nc'
FILE_PV_S4 = 'PV_index.nc'
hgt_s4 = xr.open_dataset(ROUTE + FILE_HGT_S4)
hgt_s4 = hgt_s4 - hgt_s4.mean(dim='longitude')
FILE_HGT_ERAI = 'hgt_erai_200.nc4'
FILE_NINIO_ERAI = 'ninio34_erai_index.nc'
FILE_PV_ERAI = 'PV_index_erai.nc4'
hgt_erai = xr.open_dataset(ROUTE + FILE_HGT_ERAI)
hgt_erai.time.values = hgt_erai.valid_time.values
hgt_erai.z.values = hgt_erai.z.values / 10
hgt_erai = hgt_erai.sel(**{'time':slice('1981-08-01', '2017-02-01')})
hgt_erai = hgt_erai - hgt_erai.mean(dim='longitude')
#abrir archivo de sst y PV del erai
ninio34_erai =  xr.open_dataset(ROUTE + FILE_NINIO_ERAI)
PV_erai =  xr.open_dataset(ROUTE + FILE_PV_ERAI)
#select ninio, ninia and neutral
index_ninio_erai_upper = ninio34_erai.ninio34_mon >= ninio34_erai.ninio34_mon.quantile(0.75, dim='dim_0')
index_PV_erai_upper = PV_erai.PV_mon >= PV_erai.PV_mon.quantile(0.75, dim='dim_0')

index_ninio_erai_lower = ninio34_erai.ninio34_mon <= ninio34_erai.ninio34_mon.quantile(0.25, dim='dim_0')
index_PV_erai_lower = PV_erai.PV_mon <= PV_erai.PV_mon.quantile(0.25, dim='dim_0')

index_ninio_erai_normal = np.logical_and(ninio34_erai.ninio34_mon < ninio34_erai.ninio34_mon.quantile(0.75, dim='dim_0'), ninio34_erai.ninio34_mon > ninio34_erai.ninio34_mon.quantile(0.25, dim='dim_0'))
index_PV_erai_normal = np.logical_and(PV_erai.PV_mon < PV_erai.PV_mon.quantile(0.75, dim='dim_0'), PV_erai.PV_mon > PV_erai.PV_mon.quantile(0.25, dim='dim_0'))

#idem s4
#abrir archivo de sst y PV del erai
ninio34_s4 =  xr.open_dataset(ROUTE + FILE_NINIO_S4)
PV_s4 =  xr.open_dataset(ROUTE + FILE_PV_S4)

#select ninio, ninia and neutral
index_ninio_s4_upper = ninio34_s4.ninio34_mon >= ninio34_s4.ninio34_mon.quantile(0.75, dim='dim_0')
index_PV_s4_upper = PV_s4.PV_mon >= PV_s4.PV_mon.quantile(0.75, dim='dim_0')

index_ninio_s4_lower = ninio34_s4.ninio34_mon <= ninio34_s4.ninio34_mon.quantile(0.25, dim='dim_0')
index_PV_s4_lower = PV_s4.PV_mon <= PV_s4.PV_mon.quantile(0.25, dim='dim_0')

index_ninio_s4_normal = np.logical_and(ninio34_s4.ninio34_mon < ninio34_s4.ninio34_mon.quantile(0.75, dim='dim_0'), ninio34_s4.ninio34_mon > ninio34_s4.ninio34_mon.quantile(0.25, dim='dim_0'))
index_PV_s4_normal = np.logical_and(PV_s4.PV_mon < PV_s4.PV_mon.quantile(0.75, dim='dim_0'), PV_s4.PV_mon > PV_s4.PV_mon.quantile(0.25, dim='dim_0'))

lon_s4, lat_s4 = np.meshgrid(hgt_s4.longitude.values, hgt_s4.latitude.values)
lon_erai, lat_erai = np.meshgrid(hgt_erai.longitude.values, hgt_erai.latitude.values)
#month = ['Aug', 'Sep', 'Oct', 'Nov']
#
#for i in np.arange(0, 4):
#	title = 'HGT 200-hPa Mean (contour) and SD (shaded) - ' + month[i]
#	filename = './figures/hgt200_mu_sd_' + month[i] + '.png'
#	PlotMeanStd(hgt_s4_mean.z.values[i, :, :], hgt_s4_sd.z.values[i, :, :], lat_s4,
#		    lon_s4, hgt_erai_mean.z.values[7 + i, :, :],
#		    hgt_erai_sd.z.values[7 + i, :, :], lat_erai, lon_erai, title, filename)
#
#seasonal means
season = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lmonth = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec']
#loop over seasons, selec data and performs all composites
for i in np.arange(0, 5):
	hgt_erai_seas_mean = hgt_erai['z'].resample(time='QS-' + lmonth[i]).mean(dim='time',skipna=True)
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	hgt_erai_smean = hgt_erai_seas_mean.sel(time= (hgt_erai_seas_mean['time.month'] == mes))
	hgt_erai_EN = np.nanmean(hgt_erai_smean.values[index_ninio_erai_upper.values, :, :], axis=0)
	hgt_erai_LN = np.nanmean(hgt_erai_smean.values[index_ninio_erai_lower.values, :, :], axis=0)
	hgt_erai_N = np.nanmean(hgt_erai_smean.values[index_ninio_erai_normal.values, :, :], axis=0)
	hgt_erai_WSPV = np.nanmean(hgt_erai_smean.values[index_PV_erai_upper.values, :, :], axis=0)
	hgt_erai_SSPV = np.nanmean(hgt_erai_smean.values[index_PV_erai_lower.values, :, :], axis=0)
	hgt_erai_NSPV = np.nanmean(hgt_erai_smean.values[index_PV_erai_normal.values, :, :], axis=0)
	hgt_s4_smean = np.nanmean(hgt_s4.z.values[i:i + 3, :, :, :], axis=0)	
	hgt_s4_EN = np.nanmean(hgt_s4_smean[index_ninio_s4_upper.values, :, :], axis=0)
	hgt_s4_LN = np.nanmean(hgt_s4_smean[index_ninio_s4_lower.values, :, :], axis=0)
	hgt_s4_N = np.nanmean(hgt_s4_smean[index_ninio_s4_normal.values, :, :], axis=0)
	hgt_s4_WSPV = np.nanmean(hgt_s4_smean[index_PV_s4_upper.values, :, :], axis=0)
	hgt_s4_SSPV = np.nanmean(hgt_s4_smean[index_PV_s4_lower.values, :, :], axis=0)
	hgt_s4_NSPV = np.nanmean(hgt_s4_smean[index_PV_s4_normal.values, :, :], axis=0)
	#composites EN-N
	title = 'Composites Z* 200hPa NINIO years - ' + season[i]
	filename = './figures/Composites_NINIO_' + season[i] + '.png'
	model = hgt_s4_EN - hgt_s4_N
	obs = hgt_erai_EN - hgt_erai_N
	PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
	#composites LN-N
	title = 'Composites Z* 200hPa NINIA years - ' + season[i]
	filename = './figures/Composites_NINIA_' + season[i] + '.png'
	model = hgt_s4_LN - hgt_s4_N
	obs = hgt_erai_LN - hgt_erai_N
	PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
	#composites EN-LN
	title = 'Composites Z* 200hPa NINIO-NINIA years - ' + season[i]
	filename = './figures/Composites_ENSO_' + season[i] + '.png'
	model = hgt_s4_EN - hgt_s4_LN
	obs = hgt_erai_EN - hgt_erai_LN
	PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
	#composites SSPV
	title = 'Composites Z* 200hPa Strong SPV years - ' + season[i]
	filename = './figures/Composites_SSPV_' + season[i] + '.png'
	model = hgt_s4_SSPV - hgt_s4_NSPV
	obs = hgt_erai_SSPV - hgt_erai_NSPV
	PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
	#composites WSPV
	title = 'Composites Z* 200hPa Weak SPV years - ' + season[i]
	filename = './figures/Composites_WSPV_' + season[i] + '.png'
	model = hgt_s4_WSPV - hgt_s4_NSPV
	obs = hgt_erai_WSPV - hgt_erai_NSPV
	PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)
	#composites SSPV - WsPV
	title = 'Composites Z* 200hPa Strong-Weak SPV years - ' + season[i]
	filename = './figures/Composites_SSPV_' + season[i] + '.png'
	model = hgt_s4_SSPV - hgt_s4_WSPV
	obs = hgt_erai_SSPV - hgt_erai_WSPV
	PlotComposites(model, lat_s4, lon_s4, obs, lat_erai, lon_erai, title, filename)

