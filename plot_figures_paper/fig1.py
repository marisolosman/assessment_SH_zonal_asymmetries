#test mean and std differences between ERA Interim and S4 HGT 200hPa
import numpy as np
import numpy.ma as ma
from scipy.stats import f
from scipy.stats import t
import datetime
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs	
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

def PlotMeanTest(t, F, lat_observ, lon_observ, title, filename):
	mes = ['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
	theta = np.linspace(0, 2 * np.pi, 100)
	center, radius = [0.5, 0.5], 0.5
	verts = np.vstack([np.sin(theta), np.cos(theta)]).T
	circle = mpath.Path(verts * radius + center)
	proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
	xticks = [-180, -90, 0, 90, 180]
	yticks = [-90, -65, -40, -15, 0]
	month = ['September', 'October', 'November', 'December', 'January', 'February']
	clevs = np.arange(-150, 175, 25)
	barra = plt.cm.get_cmap('RdBu_r')
	fig = plt.figure(1,(8, 5.7), 300)
	for i in np.arange(6):
		ax =plt.subplot(2, 3, i + 1, projection=proj)
		ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree()) #set map limits
		gl = ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels=False, 
				 crs=ccrs.PlateCarree(), color='black', alpha=0.6,
				 linewidth=0.3, linestyle='--')
		gl.n_steps = 90	
		ax.set_boundary(circle, transform=ax.transAxes)
		ax.coastlines()
		im = ax.contourf(lon_observ, lat_observ, t[i + 1, :, :], clevs,
				transform=ccrs.PlateCarree(), cmap=barra, extend='both',
				vmin=clevs[0], vmax=clevs[-1])
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(month[i], fontsize=10)
	plt.suptitle(title, fontsize=12, x=0.52, y=0.9)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.34, 0.1, 0.35, 0.05])
	cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	cb.ax.tick_params(labelsize=8)
	plt.savefig(filename, dpi=300, bbox_inches='tight')
	plt.clf()
	plt.cla()
	plt.close()

def PlotStdTest(t, F, lat_observ, lon_observ, title, filename):
	mes = ['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
	theta = np.linspace(0, 2 * np.pi, 100)
	center, radius = [0.5, 0.5], 0.5
	verts = np.vstack([np.sin(theta), np.cos(theta)]).T
	circle = mpath.Path(verts * radius + center)
	proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
	xticks = [-180, -90, 0, 90, 180]
	yticks = [-90, -65, -40, -15, 0]
	month = ['September', 'October', 'November', 'December', 'January', 'February']
	clevs = np.arange(0, 3.5, 0.5)
	barra = plt.cm.get_cmap('YlOrRd')
	fig = plt.figure(1,(8, 5.7), 300)
	for i in np.arange(6):
		ax =plt.subplot(2, 3, i + 1, projection=proj)
		ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree()) #set map limits
		gl = ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels=False, 
			     crs=ccrs.PlateCarree(), color='black', linewidth=0.3, linestyle='--')
		gl.n_steps = 90	
		ax.set_boundary(circle, transform=ax.transAxes)
		ax.coastlines()
		im = ax.contourf(lon_observ, lat_observ, F[i + 1, :, :], clevs,
				transform=ccrs.PlateCarree(), cmap=barra, extend='max',
				vmin=clevs[0], vmax=clevs[-1])
		barra.set_over(barra(barra.N-1))
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(month[i], fontsize=10)
	plt.suptitle(title, fontsize=12, x=0.52, y=0.9)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.39, 0.1, 0.25, 0.05])
	cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	cb.ax.tick_params(labelsize=8)
	plt.savefig(filename, dpi=300, bbox_inches='tight')
	plt.clf()
	plt.cla()
	plt.close()

ROUTE = '~/datos/data/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures_paper/'

if not(os.path.isfile(ROUTE + 'datos_fig1.npz')):
	FILE = 'monthly_hgt200_aug_feb.nc4'
	ERAI = 'hgt_erai_200.nc4'
	hgt_erai = xr.open_dataset(ROUTE + ERAI)
	hgt_erai.time.values = hgt_erai.valid_time.values
	hgt_erai = xr.decode_cf(hgt_erai)
	index_time = np.logical_or(hgt_erai.time.values<=np.datetime64('2002-07-31'), hgt_erai.time.values>=np.datetime64('2003-03-01'))
	#remove data from august 2002 to february 2003
	hgt_erai = hgt_erai.sel(time=index_time)
	hgt_erai_mean = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')}).groupby('time.month').mean(dim='time')
	hgt_erai_sd = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')}).groupby('time.month').std(dim='time')

	hgt_s4 = xr.open_dataset(ROUTE + FILE)
	hgt_s4_mean = hgt_s4.mean(dim='realiz')
	hgt_s4_sd = hgt_s4.std(dim='realiz')

	lon_s4, lat_s4 = np.meshgrid(hgt_s4.longitude.values, hgt_s4.latitude.values)
	lon_erai, lat_erai = np.meshgrid(hgt_erai.longitude.values, hgt_erai.latitude.values)
	month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
	mm = [7, 8, 9, 10, 11, 0, 1]

	tt = np.ma.empty([7, lat_erai.shape[0], lat_erai.shape[1]])
	F = np.ma.empty([7, lat_erai.shape[0], lat_erai.shape[1]])

	for i in np.arange(0, 7):
		#interp model grid to obs grid and compute test
		hgt_s4_mean_new = hgt_s4_mean.z[i,:,:].interp_like(hgt_erai_mean.z[mm[i], :, :], method='linear')
		hgt_s4_sd_new = hgt_s4_sd.z[i,:,:].interp_like(hgt_erai_sd.z[mm[i], :, :], method='linear')
		SE= np.sqrt(np.power(hgt_erai_sd.z.values[mm[i], :, :], 2)/ 36 + np.power(hgt_s4_sd_new.values, 2)/ (36 * 51))
		tt[i, :, :] = (hgt_s4_mean_new.values - hgt_erai_mean.z.values[mm[i], :, :])
		tt[np.isnan(tt)] = 0
		DF = 36 + 51 * 36 - 2
		#DF = np.power(SE, 4) / (np.power(np.power(hgt_erai_sd[7 + i, :, :], 2)/ 36) / 35 + 
		#			np.power(np.power(hgt_s4_sd_new, 2)/ (36 * 51)) /(36 * 51 - 1))
		t_cut = t.ppf(0.025, DF)
		df_cut = f.ppf([0.025, 0.975], 35, 36 * 51 -1)
		F[i, :, :] = np.power (hgt_erai_sd.z.values[mm[i], :, :], 2) / np.power(hgt_s4_sd_new.values, 2)
		F[np.isnan(F)] = 1
		tt[i, :, :] = ma.masked_array(tt[i, :, :], mask=np.logical_and(np.greater(tt[i, :, :]/SE, t_cut),
				 np.less(tt[i, :, :]/SE, np.abs(t_cut))))
		F[i, :, :] = ma.masked_array(F[i, :, :], mask=np.logical_and(F[i, :, :]>df_cut[0],
					     F[i, :, :]<df_cut[1]))
	title = 'GPH 200hPa Model Bias (m)'
	filename = FIG_PATH + 'fig1.eps'
	PlotMeanTest(tt, F, lat_erai, lon_erai, title, filename)

	title = 'GPH 200hPa Variance Ratio'
	filename = FIG_PATH + 'fig2.eps'
	PlotStdTest(tt, F, lat_erai, lon_erai, title, filename)
	np.savez('datos_fig1.npz', tt=tt, F=F, lat_erai=lat_erai, lon_erai=lon_erai,
		 title=title, filename=filename)
else:
	var = np.load('datos_fig1.npz')
	title = 'GPH 200hPa Model Bias (m)'
	filename = FIG_PATH + 'fig1.eps'
	PlotMeanTest(var['tt'], var['F'], var['lat_erai'], var['lon_erai'], title, filename)

	title = 'GPH 200hPa Variance Ratio'
	filename = FIG_PATH + 'fig2.eps'
	PlotStdTest(var['tt'], var['F'], var['lat_erai'], var['lon_erai'], title, filename)

