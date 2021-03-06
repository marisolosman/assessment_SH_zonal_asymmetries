#test mean and std differences between ERA Interim and S4 HGT 200hPa
import numpy as np
import numpy.ma as ma
from scipy.stats import f
from scipy.stats import t
import datetime
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

def PlotMeanStdTest(t, F, lat_observ, lon_observ, title, filename):
	proj = ccrs.PlateCarree(central_longitude=180) #define proj: cyl eq
	fig = plt.figure(1,(10, 6.7), 300)
	ax =plt.subplot(2, 1, 1, projection=proj)
	clevs = np.arange(-9, 11, 2)
	barra = plt.cm.get_cmap('RdBu_r')
	ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree()) #set map limits
	im = ax.contourf(lon_observ, lat_observ, t, clevs,
			 transform=ccrs.PlateCarree(), cmap=barra, extend='both',
			 vmin=clevs[0], vmax=clevs[-1]) 
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	plt.colorbar(im, orientation='vertical')
	ax.coastlines()
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	plt.title('Bias in Mean', fontsize=10)
	if np.sum(np.sum(F.mask))<114000:
		ax =plt.subplot(2, 1, 2, projection=proj)
		clevs = np.arange(0, 3.5, 0.5)
		barra = plt.cm.get_cmap('YlOrRd')
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon_observ, lat_observ, F, clevs,
				 transform=ccrs.PlateCarree(), cmap=barra, extend='both',
				 vmin=0, vmax=3)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		plt.colorbar(im, orientation='vertical')
		ax.add_feature(cartopy.feature.COASTLINE)
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title('Bias in Variance', fontsize=10)
	plt.suptitle(title, fontsize=12, x=0.44, y=0.9)
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

ROUTE = '~/datos/data/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/model_assessment/'
FILE = './fogt/monthly_winds50_aug_feb.nc4'
winds_s4 = xr.open_dataset(ROUTE + FILE, chunks={'realiz':10})
winds_s4_mean = winds_s4.mean(dim='realiz')
winds_s4_sd = winds_s4.std(dim='realiz')

ERAI = 'winds_erai_50.nc4'
winds_erai = xr.open_dataset(ROUTE + ERAI)
winds_erai.time.values = winds_erai.valid_time.values
#remove data from august 2002 to february 2003
winds_erai = winds_erai.isel(time=np.logical_or(winds_erai.time.values<=np.datetime64('2002-07-31'), winds_erai.time.values>=np.datetime64('2003-03-01')))
winds_erai_mean = winds_erai.groupby('time.month').mean(dim='time')
winds_erai_sd = winds_erai.groupby('time.month').std(dim='time')

lon_s4, lat_s4 = np.meshgrid(winds_s4.longitude.values, winds_s4.latitude.values)
lon_erai, lat_erai = np.meshgrid(winds_erai.longitude.values, winds_erai.latitude.values)
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
mm = [7, 8, 9, 10, 11, 0, 1]

for i in np.arange(0, 7):
	title = 'Model Bias - ' + month[i]
	filename = FIG_PATH + 'u50_bias_' + month[i] + '.png'
	#interp model grid to obs grid and compute test
	winds_s4_mean_new = winds_s4_mean.u[i,:,:].interp_like(winds_erai_mean.u[mm[i], :, :], method='linear')
	winds_s4_sd_new = winds_s4_sd.u[i,:,:].interp_like(winds_erai_sd.u[mm[i], :, :], method='linear')
	SE= np.sqrt(np.power(winds_erai_sd.u.values[mm[i], :, :], 2)/ 36 + np.power(winds_s4_sd_new.values, 2)/ (36 * 51))
	tt = (winds_s4_mean_new.values - winds_erai_mean.u.values[mm[i], :, :])
	DF = 36 + 51 * 36 - 2
	#DF = np.power(SE, 4) / (np.power(np.power(winds_erai_sd[7 + i, :, :], 2)/ 36) / 35 + 
	#			np.power(np.power(winds_s4_sd_new, 2)/ (36 * 51)) /(36 * 51 - 1))
	t_cut = t.ppf(0.025, DF)
	df_cut = f.ppf([0.025, 0.975], 35, 36 * 51 -1)
	F = np.power (winds_erai_sd.u.values[mm[i], :, :], 2) / np.power(winds_s4_sd_new.values, 2)
	tt = ma.masked_array(tt, mask=np.logical_and((tt/SE)>t_cut, (tt/SE)<np.abs(t_cut)))
	F = ma.masked_array(F, mask=np.logical_and(F>df_cut[0], F<df_cut[1]))
	PlotMeanStdTest(tt, F, lat_erai, lon_erai, title, filename)

#seasonal means
season = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lmonth = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec']

for i in np.arange(0, 5):
	winds_erai_seas_mean = winds_erai['u'].sel(**{'time':slice('1981-08-01',
							       '2018-02-01')}).resample(time='QS-' + lmonth[i]).mean(dim='time')
	winds_s4_smean_mean = winds_s4.u[i:i + 3, :, :, :].mean(dim='month').mean(dim='realiz')
	winds_s4_smean_sd = winds_s4.u[i:i + 3, :, :, :].mean(dim='month').std(dim='realiz')
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	winds_erai_smean = winds_erai_seas_mean.sel(time= (winds_erai_seas_mean['time.month'] == mes))
	winds_erai_smean_mean = winds_erai_smean.mean(dim='time', skipna='True')
	winds_erai_smean_sd = winds_erai_smean.std(dim='time', skipna='True')
	title = 'Model Bias - ' + season[i]
	filename = FIG_PATH + 'u50_bias_' + season[i] + '.png'
	#interpolation
	winds_s4_smean_mean_new = winds_s4_smean_mean.interp_like(winds_erai_smean_mean, method='linear')
	winds_s4_smean_sd_new = winds_s4_smean_sd.interp_like(winds_erai_smean_sd, method='linear')
	SE= np.sqrt(np.power(winds_erai_smean_sd.values, 2)/ 36 + np.power(winds_s4_smean_sd_new.values, 2)/ (36 * 51))
	tt = (winds_s4_smean_mean_new.values - winds_erai_smean_mean.values)
	DF = 36 + 51 * 36 - 2
	#DF = np.power(SE, 4) / (np.power(np.power(winds_erai_sd[7 + i, :, :], 2)/ 36) / 35 + 
	#			np.power(np.power(winds_s4_sd_new, 2)/ (36 * 51)) /(36 * 51 - 1))
	t_cut = t.ppf(0.025, DF)
	df_cut = f.ppf([0.025, 0.975], 35, 36 * 51 -1)
	F = np.power (winds_erai_smean_sd.values, 2) / np.power(winds_s4_smean_sd_new.values, 2)
	tt = ma.masked_array(tt, mask=np.logical_and((tt/SE)>t_cut, (tt/SE)<np.abs(t_cut)))
	F = ma.masked_array(F, mask=np.logical_and(F>df_cut[0], F<df_cut[1]))
	PlotMeanStdTest(tt, F, lat_erai, lon_erai, title, filename)


