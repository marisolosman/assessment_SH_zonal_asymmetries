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

def PlotMeanStd(model_m, model_sd, lat_model, lon_model, observ_m, observ_sd,
		lat_observ, lon_observ,	title, filename):
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1,(10, 6.7), 300)
	ax =plt.subplot(2, 1, 1, projection=proj)
	clevs = np.arange(0, 160, 20)
	barra = plt.cm.get_cmap('YlOrRd')
	ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
	im = ax.contourf(lon_observ, lat_observ, observ_sd, clevs,
			 transform=ccrs.PlateCarree(), cmap=barra, extend='both',
			 vmin=0, vmax=160)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	ax.coastlines()
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	ax.contour(lon_observ, lat_observ, observ_m,
		   levels=np.arange(8000, 14000, 400), linewidths=0.4,
		   colors='k', transform=ccrs.PlateCarree())
	plt.title('ERA Interim', fontsize=10)
	ax =plt.subplot(2, 1, 2, projection=proj)
	clevs = np.arange(0, 160, 20)
	barra = plt.cm.get_cmap('YlOrRd')
	ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
	im = ax.contourf(lon_model, lat_model, model_sd, clevs,
			 transform=ccrs.PlateCarree(), cmap=barra, extend='both',
			 vmin=0, vmax=160)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	ax.add_feature(cartopy.feature.COASTLINE)
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	ax.contour(lon_model, lat_model, model_m,
		   levels=np.arange(8000, 14000, 400), linewidths=0.4,
		   colors='k', transform=ccrs.PlateCarree())
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
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/model_assessment/'
ERAI = 'hgt_erai_200.nc4'
hgt_erai = xr.open_dataset(ROUTE + ERAI)
hgt_erai.time.values = hgt_erai.valid_time.values
#remove data from august 2002 to february 2003
hgt_erai = hgt_erai.sel(time=hgt_erai.time.values[np.logical_or(hgt_erai.time.values<=np.datetime64('2002-07-31'), hgt_erai.time.values>=np.datetime64('2003-03-01'))])
hgt_erai_mean = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')}).groupby('time.month').mean(dim='time', skipna='True')
hgt_erai_sd = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')}).groupby('time.month').std(dim='time', skipna='True')


FILE = 'monthly_hgt200_aug_feb.nc4'
hgt_s4 = xr.open_dataset(ROUTE + FILE)
hgt_s4_mean = hgt_s4.mean(dim='realiz')
hgt_s4_sd = hgt_s4.std(dim='realiz')


lon_s4, lat_s4 = np.meshgrid(hgt_s4.longitude.values, hgt_s4.latitude.values)
lon_erai, lat_erai = np.meshgrid(hgt_erai.longitude.values, hgt_erai.latitude.values)
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
mm = [7, 8, 9, 10, 11, 0, 1]
for i in np.arange(0, 7):
	title = 'HGT 200-hPa Mean (contour) and SD (shaded) - ' + month[i]
	filename = FIG_PATH + 'hgt200_mu_sd_' + month[i] + '.png'
	PlotMeanStd(hgt_s4_mean.z.values[i, :, :], hgt_s4_sd.z.values[i, :, :], lat_s4,
		    lon_s4, hgt_erai_mean.z.values[mm[i], :, :],
		    hgt_erai_sd.z.values[mm[i], :, :], lat_erai, lon_erai, title, filename)

#seasonal means
season = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lmonth = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec']
for i in np.arange(0, 5):
	hgt_erai_seas_mean = hgt_erai['z'].sel(**{'time':slice('1981-08-01',
							       '2018-02-01')}).resample(time='QS-' + lmonth[i]).mean(dim='time')
	hgt_s4_smean = np.mean(hgt_s4.z.values[i:i + 3, :, :, :], axis=0)	
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	hgt_erai_smean = hgt_erai_seas_mean.sel(time= (hgt_erai_seas_mean['time.month'] == mes))
	title = 'HGT 200-hPa Mean (contour) and SD (shaded) - ' + season[i]
	filename = FIG_PATH + 'hgt200_mu_sd_' + season[i] + '.png'
	PlotMeanStd(np.mean(hgt_s4_smean, axis=0), np.std(hgt_s4_smean, axis=0), lat_s4,
		    lon_s4, np.nanmean(hgt_erai_smean, axis=0), np.nanstd(hgt_erai_smean, axis=0),
		    lat_erai, lon_erai, title, filename)



