#sort el ninio events
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def PlotComposites(var_pos, var_neg, lat, lon, title, filename):
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	ax = plt.subplot(2, 1, 1, projection=proj)
	clevs = np.arange(-60, 70, 10)
	barra = plt.cm.RdBu_r
	ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
	im = ax.contourf(lon, lat, var_pos, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	ax.coastlines()
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	plt.title('+NINIO3.4 YEARS')
	
	ax1 = plt.subplot(2, 1, 2, projection=proj)
	barra = plt.cm.RdBu_r
	ax1.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
	im = ax1.contourf(lon, lat, var_neg, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	ax1.coastlines()
	ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
	ax1.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax1.xaxis.set_major_formatter(lon_formatter)
	ax1.yaxis.set_major_formatter(lat_formatter)
	plt.title('-NINIO3.4 YEARS')
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
def PlotCompDiff(var, lat, lon, title, filename):
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 3.7), 300)
	ax = plt.subplot(projection=proj)
	clevs = np.arange(-60, 70, 10)
	barra = plt.cm.RdBu_r
	ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
	im = ax.contourf(lon, lat, var, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	ax.coastlines()
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	plt.colorbar(im, orientation='horizontal')
	plt.title(title)
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()


ROUTE = '~/data/'

#abro el archivo de geopotencial

hgt_erai = xr.open_dataset(ROUTE + 'hgt200_erai.nc')
hgt_erai.time.values = hgt_erai.valid_time.values
hgt_erai.z.values = hgt_erai.z.values / 10
nlat, nlon = np.shape(hgt_erai.z.values)[1:]
hgt_erai =  hgt_erai.sel(time = (hgt_erai['time.year']!=2002))
hgt_erai = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')})
#rearange hgt values to math those sst
hgt_erai_monthly = hgt_erai.sel(time = np.logical_and(hgt_erai['time.month'] >= 8, hgt_erai['time.month'] <= 11))
hgt_aso_erai = hgt_erai.resample(time='QS-Aug').mean(dim='time', skipna='True')
hgt_son_erai = hgt_erai.resample(time='QS-Sep').mean(dim='time', skipna='True')

hgt_erai_seasonal = xr.concat([hgt_aso_erai.sel(time=(hgt_aso_erai['time.month'] == 8)), hgt_son_erai.sel(time=(hgt_son_erai['time.month'] == 9))], dim='time').dropna(dim='time')

ninio34 = xr.open_dataset(ROUTE + 'ninio34_erai_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']

#print(ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'))
#print(ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

index_monthly_upper = ninio34.ninio34_mon >= ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34.ninio34_mon <= ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))
#rearange hgt values to use the sst index
hgt_new = np.transpose(np.reshape(hgt_erai_monthly.z.values, [36, 4, nlat, nlon]), [1, 0, 2, 3])
for i in np.arange(0,4):
	var_neg = np.mean(hgt_new[i, index_monthly_lower.values, :, :], axis=0) - np.mean(hgt_new[i, index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(hgt_new[i, index_monthly_upper.values, :, :], axis=0) - np.mean(hgt_new[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites ERAI HGT 200hPa - ' + month[i]
	filename = './figures/hgt_200_erai_composites_NINIO_' + month[i] +'_3.png'
	PlotComposites(var_pos, var_neg, hgt_erai.latitude, hgt_erai.longitude, tit, filename)
	var = np.mean(hgt_new[i, index_monthly_upper.values, :, :], axis=0) - np.mean(hgt_new[i, index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites differences EN-LN Years - ' + month[i]
	filename = './figures/hgt_200_erai_composites_diff_NINIO_' + month[i] +'_3.png'
	PlotCompDiff(var, hgt_erai.latitude, hgt_erai.longitude, tit, filename)

hgt_seas_new = np.reshape(hgt_erai_seasonal.z.values, [2, 36, nlat, nlon])
for i in np.arange(0,2):
	var_neg = np.mean(hgt_seas_new[i, index_monthly_lower.values, :, :], axis=0) - np.mean(hgt_seas_new[i, index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(hgt_seas_new[i, index_monthly_upper.values, :, :], axis=0) - np.mean(hgt_seas_new[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites ERAI HGT 200hPa - ' + seas[i]
	filename = './figures/hgt_200_erai_composites_NINIO_' + seas[i] +'_4.png'
	PlotComposites(var_pos, var_neg, hgt_erai.latitude, hgt_erai.longitude, tit, filename)
	var = np.mean(hgt_seas_new[i, index_monthly_upper.values, :, :], axis=0) - np.mean(hgt_seas_new[i, index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites differences EN-LN Years - ' + seas[i]
	filename = './figures/hgt_200_erai_composites_diff_NINIO_' + seas[i] +'_4.png'
	PlotCompDiff(var, hgt_erai.latitude, hgt_erai.longitude, tit, filename)

