#composites of PV events conditioned on ninio events
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def PlotZComposites(var_pos, var_neg, lat, lon, title, filename):
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
def PlotZCompDiff(var, lat, lon, title, filename):
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


RUTA='~/data/'
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'monthly_hgt200_aug_feb.nc')
ds_zonal_mean = ds.mean('longitude')
ds = ds - ds_zonal_mean
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
PV_index = xr.open_dataset(RUTA + 'PV_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

#search for years with EN events
index_monthly_upper = ninio34.ninio34_mon >= ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear')

# compute PV composites conditioned on EN anomalies
PV_index_EN = PV_index.sel(dim_0 = index_monthly_upper.values)
ds_EN = ds.sel(realiz = index_monthly_upper.values)

index_monthly_upper = PV_index_EN.PV_mon >= PV_index_EN.PV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
index_monthly_lower = PV_index_EN.PV_mon <= PV_index_EN.PV_mon.quantile(0.25, dim='dim_0', interpolation='linear')
for i in np.arange(0,7):
	var = np.mean(ds_EN.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_EN.z.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Z* 200hPa Composites differences SPW-WPV Years - ' + month[i] + ' - EN'
	filename = './figures/z_200_composites_diff_PV_' + month[i] +'_EN.png'
	PlotZCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,5):
	var = ds_EN.isel(month=range(i, i + 3)).mean(dim='month')
	var = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z[index_monthly_upper.values, :, :], axis=0)
	tit = 'Z* 200hPa Composites differences SPV-WPV Years - ' + seas[i] + ' - EN'
	filename = './figures/z_200_composites_diff_PV_' + seas[i] +'_EN.png'
	PlotZCompDiff(var, ds.latitude, ds.longitude, tit, filename)

#search for years with LN events
index_monthly_lower = ninio34.ninio34_mon <= ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear')

# compute PV composites conditioned on LN anomalies
PV_index_EN = PV_index.sel(dim_0 = index_monthly_lower.values)
ds_EN = ds.sel(realiz = index_monthly_lower.values)

index_monthly_upper = PV_index_EN.PV_mon >= PV_index_EN.PV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
index_monthly_lower = PV_index_EN.PV_mon <= PV_index_EN.PV_mon.quantile(0.25, dim='dim_0', interpolation='linear')

for i in np.arange(0,7):
	var = np.mean(ds_EN.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_EN.z.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Z* 200hPa Composites differences SPV-WPV Years - ' + month[i] + ' - LN'
	filename = './figures/z_200_composites_diff_PV_' + month[i] +'_LN.png'
	PlotZCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,5):
	var = ds_EN.isel(month=range(i, i + 3)).mean(dim='month')
	var = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z[index_monthly_upper.values, :, :], axis=0)
	tit = 'Z* 200hPa Composites differences SPV-WPV Years - ' + seas[i] + ' - LN'
	filename = './figures/z_200_composites_diff_PV_' + seas[i] +'_LN.png'
	PlotZCompDiff(var, ds.latitude, ds.longitude, tit, filename)

