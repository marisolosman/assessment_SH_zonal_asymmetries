#plot correlation between eof and sst for ERA Interim and S4
import numpy as np
import numpy.ma as ma
import datetime
import pandas as pd
import xarray as xr
from eofs.standard import Eof
from scipy import signal
from scipy.stats import t
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

def PlotCorrelation(correlat, lat_observ, lon_observ, title, filename):
	rcorte = t.ppf(0.025, 36-2)/np.sqrt(36 - 2 + t.ppf(0.025, 36-2) ** 2)
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1,(10, 6), 300)
	for i in range(3):
		ax =plt.subplot(3, 1, i + 1, projection=proj)
		clevs = np.linspace(-1, 1, 11)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon_observ, lat_observ, correlat[i, :, :], clevs,
				 transform=ccrs.PlateCarree(), cmap=barra, vmin=-1, vmax=1)
		ax.contour(lon_observ, lat_observ, correlat[i, :, :], [rcorte, -rcorte],
			   colors='k', linewidths=0.5, transform=ccrs.PlateCarree())
		ax.coastlines()
		ax.gridlines()#crs=proj, linewidth=0.3, linestyle='-')
		plt.title('EOF ' + str(i + 1), fontsize=10)
	plt.suptitle(title, fontsize=12)
	fig.subplots_adjust(bottom=0.17, top=0.87, hspace=0.2, wspace=0.07)
	cbar_ax = fig.add_axes([0.375, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()
ROUTE = '~/datos/data/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/model_assessment/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
hgt_s4 = xr.open_dataset(ROUTE + FILE_HGT_S4)
hgt_s4 = hgt_s4.sel(**{'latitude':slice(-20, -90)})
FILE_HGT_ERAI = 'hgt_erai_200.nc4'
hgt_erai = xr.open_dataset(ROUTE + FILE_HGT_ERAI)
hgt_erai.time.values = hgt_erai.valid_time.values
#remove data from august 2002 to february 2003
hgt_erai = hgt_erai.sel(time=hgt_erai.time.values[np.logical_or(hgt_erai.time.values<=np.datetime64('2002-07-31'), hgt_erai.time.values>=np.datetime64('2003-03-01'))])
hgt_erai = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')})
hgt_erai = hgt_erai.sel(**{'latitude':slice(-20, -90)})
lon_s4, lat_s4 = np.meshgrid(hgt_s4.longitude.values, hgt_s4.latitude.values)
lon_erai, lat_erai = np.meshgrid(hgt_erai.longitude.values, hgt_erai.latitude.values)
#open observed and forecast sst
FILE_SST_ERAI = 'fogt/sst_erai.nc4'
FILE_SST_S4 = 'fogt/sst_s4_aug_feb.nc4'
sst_erai = xr.open_dataset(ROUTE + FILE_SST_ERAI)
sst_erai.time.values = sst_erai.valid_time.values
#remove data from august 2002 to february 2003
sst_erai = sst_erai.sel(time=sst_erai.time.values[np.logical_or(sst_erai.time.values<=np.datetime64('2002-07-31'), sst_erai.time.values>=np.datetime64('2003-03-01'))])
sst_erai = sst_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')})
sst_erai = sst_erai.sel(**{'latitude':slice(20, -90)})
sst_s4 = xr.open_dataset(ROUTE + FILE_SST_S4)
sst_s4 = sst_s4.sel(**{'latitude':slice(20, -90)})
lon_s4_sst, lat_s4_sst = np.meshgrid(sst_s4.longitude.values, sst_s4.latitude.values)
lon_erai_sst, lat_erai_sst = np.meshgrid(sst_erai.longitude.values, sst_erai.latitude.values)

#seasonal means
season = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lmonth = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec']
#loop over seasons, selec data and performs all composites
NMODES = 3
tcc = np.empty([NMODES, 10])
s4_pc_sorted = np.empty([5, 36, 3])
s4_eof_sorted = np.empty([5, 3, np.shape(lat_s4)[0], np.shape(lat_s4)[1]])
s4_eval_sorted = np.empty([5, 3])
s4_realiz_pc_sorted = np.empty([5, 36 * 51, 3])
s4_realiz_eof_sorted = np.empty([5, 3, np.shape(lat_s4)[0], np.shape(lat_s4)[1]])
s4_realiz_eval_sorted = np.empty([5, 3])
pc_erai = np.empty([5, 36, 3])
corr_sst_s4 = np.empty([5, 3, np.shape(lat_s4_sst)[0], np.shape(lat_s4_sst)[1]])
corr_sst_erai = np.empty([5, 3, np.shape(lat_erai_sst)[0], np.shape(lat_erai_sst)[1]])

for i in np.arange(0, 5):
	hgt_erai_seas_mean = hgt_erai['z'].resample(time='QS-' + lmonth[i]).mean(dim='time',skipna=True)
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	hgt_erai_smean = hgt_erai_seas_mean.sel(time= np.logical_and(hgt_erai_seas_mean['time.month'] == mes, hgt_erai_seas_mean['time.year']!=2002))
	hgt_s4_smean = np.nanmean(hgt_s4.z.values[i:i + 3, :, :, :], axis=0)
	#detrend data
	hgt_s4_smean = signal.detrend(hgt_s4_smean, axis=0, type='linear')
	hgt_erai_smean = np.nan_to_num(hgt_erai_smean.values, np.nanmean(hgt_erai_smean.values))
	hgt_erai_smean = signal.detrend(hgt_erai_smean, axis=0, type='linear')
	#eof analysis obs
	# Compute anomalies by removing the time-mean
	z_mean = np.nanmean(hgt_erai_smean, axis=0)
	z_anom = hgt_erai_smean - z_mean
	# Create an EOF solver to do the EOF analysis. Square-root of cosine of
	# latitude weights are applied before the computation of EOFs.
	solver = Eof(z_anom)#, weights=wgts)
	eofs = solver.eofsAsCorrelation(neofs=5)
	exp_var = solver.varianceFraction()
	pcs = solver.pcs(npcs=5, pcscaling=1)
	pc_erai[i, :, :] = pcs[:, 0:3]
	#eof analysis model mean
	# Compute anomalies by removing the time-mean
	z_mean = np.nanmean(hgt_s4_smean, axis=0)
	#computo media del ensamble
	hgt_s4m_smean = np.mean(np.reshape(hgt_s4_smean, [36, 51, 99, 512]), axis=1)
	z_anom = hgt_s4m_smean - z_mean
	solver_s4 = Eof(z_anom)#, weights=wgts)
	eofs_s4 = solver_s4.eofsAsCorrelation(neofs=10)
	exp_var_s4 = solver_s4.varianceFraction()
	pcs_s4 = solver_s4.pcs(npcs=10, pcscaling=1)
	#eof analysis model realizations
	z_mean = np.nanmean(hgt_s4_smean, axis=0)
	z_anom = hgt_s4_smean - z_mean
	solver_s4_realiz = Eof(z_anom)#, weights=wgts)
	eofs_s4_realiz = solver_s4_realiz.eofsAsCorrelation(neofs=5)
	exp_var_s4_realiz = solver_s4_realiz.varianceFraction()
	pcs_s4_realiz = solver_s4_realiz.pcs(npcs=5, pcscaling=1)
	#compute correlation between observed and forecast PCs
	eofs = np.reshape(eofs,[5, np.shape(lat_erai)[0] * np.shape(lat_erai)[1] ] )
	#eofs_s4 = np.reshape(eofs_s4,[10, np.shape(lat_s4)[0] * np.shape(lat_s4)[1] ] )
	pcc = np.empty([NMODES, 10])
	tcc = np.empty([NMODES, 10])
	for j in np.arange(NMODES): #loop over observed modes
		for k in np.arange(10): #loop over forecast modes
			ds = xr.DataArray(eofs_s4[k, :, :], coords=[hgt_s4.latitude.values, hgt_s4.longitude.values], dims=['lat', 'lon'])
			dsi = ds.interp(lat=hgt_erai.latitude.values, lon=hgt_erai.longitude.values)
			tcc[j, k] = np.corrcoef(pcs[:, j], pcs_s4[:, k])[0, 1]
			mx = ma.masked_array(dsi.values, mask=np.isnan(dsi.values))			
			pcc[j, k] = ma.corrcoef(eofs[j, :], 
						   np.reshape(mx, [lat_erai.shape[0] * lat_erai.shape[1]]))[0, 1]
	##compute pss
	pss = tcc * pcc
	##find max pss for each mode
	print(season[i])
	for j in np.arange(3):
		index = np.nanargmax(pss[j, :])
		print(index, tcc[j, index], pcc[j, index])
		s4_eof_sorted[i, j, :, :] = np.sign(pcc[j, index]) * eofs_s4[index, :, :]
		s4_eval_sorted[i, j] = exp_var[index]
		s4_pc_sorted[i, :, j] = np.sign(pcc[j, index]) * pcs_s4[:, index]
		pss = np.delete(pss, index, 1)
		pcc = np.delete(pcc, index, 1)
		tcc = np.delete(tcc, index, 1)
		pcs_s4 = np.delete(pcs_s4, index, 1)
		eofs_s4 = np.delete(eofs_s4, index, 0)
		exp_var_s4 = np.delete(exp_var_s4, index)
	#compute correlation between observed and forecast PCs
	pcc = np.empty([NMODES, 10])

	eofs = np.reshape(eofs,[5, np.shape(lat_erai)[0] * np.shape(lat_erai)[1] ] )
	for j in np.arange(NMODES): #loop over observed modes
		for k in np.arange(5): #loop over forecast modes
			ds = xr.DataArray(eofs_s4_realiz[k, :, :], coords=[hgt_s4.latitude.values, hgt_s4.longitude.values], dims=['lat', 'lon'])
			dsi = ds.interp(lat=hgt_erai.latitude.values, lon=hgt_erai.longitude.values)
			mx = ma.masked_array(dsi.values, mask=np.isnan(dsi.values))			
			pcc[j, k] = ma.corrcoef(eofs[j, :], 
						   np.reshape(mx, [lat_erai.shape[0] * lat_erai.shape[1]]))[0, 1]
	##compute pss
	pss = np.abs(pcc)
	##find max pss for each mode
	print(season[i])
	for j in np.arange(3):
		index = np.nanargmax(pss[j, :])
		print(index, pcc[j,index])
		s4_realiz_eof_sorted[i, j, :, :] = np.sign(pcc[j, index]) * eofs_s4_realiz[index, :, :]
		s4_realiz_eval_sorted[i, j] = exp_var_s4_realiz[index]
		s4_realiz_pc_sorted[i, :, j] = np.sign(pcc[j, index]) * pcs_s4_realiz[:, index]
		pss = np.delete(pss, index, 1)
		pcc = np.delete(pcc, index, 1)
		pcs_s4_realiz = np.delete(pcs_s4_realiz, index, 1)
		eofs_s4_realiz = np.delete(eofs_s4_realiz, index, 0)
		exp_var_s4_realiz = np.delete(exp_var_s4_realiz, index)
	sst_erai_seas_mean = sst_erai['sst'].resample(time='QS-' + lmonth[i]).mean(dim='time',skipna=True)
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	sst_erai_smean = sst_erai_seas_mean.sel(time= np.logical_and(sst_erai_seas_mean['time.month'] == mes, sst_erai_seas_mean['time.year']!=2002))
	sst_erai_smean = np.nan_to_num(sst_erai_smean, np.nanmean(sst_erai_smean))
	sst_erai_smean = signal.detrend(sst_erai_smean, axis=0, type='linear')
	sst_s4_smean = np.nanmean(sst_s4.sst.values[i:i + 3, :, :, :], axis=0)
	sst_s4_smean = np.nanmean(np.reshape(sst_s4_smean, [157, 512, 36, 51]), axis=3)
	sst_s4_smean = np.nan_to_num(sst_s4_smean, np.nanmean(sst_s4_smean))

	sst_s4_smean = signal.detrend(sst_s4_smean, axis=2, type='linear')

	for j in np.arange(3): #loop over modes
		for k in range(np.shape(lat_s4_sst)[0]):
			for l in range(np.shape(lat_s4_sst)[1]):
				corr_sst_s4[i, j, k, l] = np.corrcoef(s4_pc_sorted[i, :, j], sst_s4_smean[k, l, :])[0, 1]
		for k in range(np.shape(lat_erai_sst)[0]):
			for l in range(np.shape(lat_erai_sst)[1]):
				corr_sst_erai[i, j, k, l] = np.corrcoef(pc_erai[i, :, j], sst_erai_smean[:, k, l])[0, 1]


	title = 'Correlation PCs vs SST - S4 - ' + season[i]
	filename = FIG_PATH + 'corr_pcs_S4_' + season[i] + '_det.png'
	PlotCorrelation(corr_sst_s4[i, :, :, :], lat_s4_sst, lon_s4_sst, title, filename)
	title = 'Correlation PCs vs SST - ERAI - ' + season[i]
	filename = FIG_PATH + 'corr_pcs_erai_' + season[i] + '_det.png'
	PlotCorrelation(corr_sst_erai[i, :, :, :], lat_erai_sst, lon_erai_sst, title, filename)

#put all eofs and pcs together
s4_realiz_pc_sorted = np.reshape(s4_realiz_pc_sorted, [5, 36, 51, 3])
eofs = np.reshape(eofs,[5, np.shape(lat_erai)[0] * np.shape(lat_erai)[1] ] )

