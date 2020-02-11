#plor meand and std for HGT 200 A-N for ERA Interim and S4
import numpy as np
import numpy.ma as ma
import datetime
import pandas as pd
import xarray as xr
from eofs.standard import Eof
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

def PlotEOF(eof, lat_observ, lon_observ, title, filename):
	proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)

	fig = plt.figure(1,(18, 7), 300)
	for i in range(3):
		ax =plt.subplot(1, 3, i + 1, projection=proj)
		clevs = np.linspace(-1, 1, 11)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, -20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon_observ, lat_observ, eof[i, :, :], clevs,
				 transform=ccrs.PlateCarree(), cmap=barra, vmin=-1, vmax=1)
		ax.coastlines()
		ax.gridlines()#crs=proj, linewidth=0.3, linestyle='-')
		plt.title('EOF ' + str(i + 1), fontsize=10)
	plt.suptitle(title, fontsize=12)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.07)
	cbar_ax = fig.add_axes([0.375, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()
def PlotPCs(obs_pcs, em_pcs, members_pcs, years, title, filename):
	fig1 = plt.figure(figsize = (15, 21), dpi = 300)  #fig size in inches
	for j in np.arange(3): #loop over Pcs
		ax = plt.subplot(3, 1, j + 1)
		ax.plot(years, obs_pcs[:, j], color = 'k', linewidth=1.5, label = 'Obs')
		ax.plot(years, em_pcs[:, j], color = 'b', linewidth=1.5, label = 'S4')
		ax.fill_between(years, np.min(members_pcs[:, :, j], 1), np.max(members_pcs[:, :, j], 1),
				facecolor='#089FFF', alpha=0.5)
		ax.axhline(y=0, xmin=0, linestyle='--', linewidth=0.8, color='b', alpha=0.4)
		ax.set_xlim((years[0] - 1, years[-1] + 1))
		ax.set_ylim((-4, 4))
		plt.title('PC '+ str(j + 1), fontsize=10)
	plt.suptitle(title ,fontsize=12, x=0.52, y=0.93)
	fig1.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')
	plt.clf()
	plt.cla()
	plt.close()
def PlotScree(var_exp, ntimes, title, filename):
	error_eval = np.sqrt(2 / ntimes) * var_exp * 100
	plt.figure
	plt.errorbar(np.arange(1,11), var_exp[0:10] * 100, error_eval[0:10],
		     color='b', linewidth=1.5)
	plt.xlabel('Modes')
	plt.ylabel('percentage variance (%)')
	plt.title(title)
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape', papertype='A4')
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
for i in np.arange(0, 5):
	hgt_erai_seas_mean = hgt_erai['z'].resample(time='QS-' + lmonth[i]).mean(dim='time',skipna=True)
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	hgt_erai_smean = hgt_erai_seas_mean.sel(time= np.logical_and(hgt_erai_seas_mean['time.month'] == mes, hgt_erai_seas_mean['time.year']!=2002))
	hgt_s4_smean = np.nanmean(hgt_s4.z.values[i:i + 3, :, :, :], axis=0)
	#eof analysis obs
	# Compute anomalies by removing the time-mean
	z_mean = np.nanmean(hgt_erai_smean.values, axis=0)
	z_anom = hgt_erai_smean.values - z_mean
	# Create an EOF solver to do the EOF analysis. Square-root of cosine of
	# latitude weights are applied before the computation of EOFs.
	solver = Eof(z_anom)#, weights=wgts)
	eofs = solver.eofsAsCorrelation(neofs=5)
	exp_var = solver.varianceFraction()
	pcs = solver.pcs(npcs=5, pcscaling=1)
	pc_erai[i, :, :] = pcs[:, 0:3]
	title = 'Observed HGT 200hPa EOFs - ' + season[i]
	filename = FIG_PATH + 'obs_eof_' + season[i] + '.png'
	PlotEOF(eofs[0:3, :, :], lat_erai, lon_erai, title, filename)
	filename = FIG_PATH + 'obs_scree_' + season[i] + '.png'
	ttle = 'Variance Explained by Observed modes - ' + season[i]
	PlotScree(exp_var, 36, title, filename)
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
	title = 'S4 HGT 200hPa EOFs - ' + season[i]
	filename = FIG_PATH + 'S4_eof_' + season[i] + '.png'
	PlotEOF(eofs_s4[0:3, :, :], lat_s4, lon_s4, title, filename)
	filename = FIG_PATH + 'S4_scree_' + season[i] + '.png'
	ttle = 'Variance Explained by S4 modes - ' + season[i]
	PlotScree(exp_var_s4, 36*51, title, filename)
	#eof analysis model realizations
	z_mean = np.nanmean(hgt_s4_smean, axis=0)
	z_anom = hgt_s4_smean - z_mean
	solver_s4_realiz = Eof(z_anom)#, weights=wgts)
	eofs_s4_realiz = solver_s4_realiz.eofsAsCorrelation(neofs=5)
	exp_var_s4_realiz = solver_s4_realiz.varianceFraction()
	pcs_s4_realiz = solver_s4_realiz.pcs(npcs=5, pcscaling=1)
	title = 'S4 realizations HGT 200hPa EOFs - ' + season[i]
	filename = FIG_PATH + 'S4_realiz_eof_' + season[i] + '.png'
	PlotEOF(eofs_s4_realiz[0:3, :, :], lat_s4, lon_s4, title, filename)
	filename = FIG_PATH + 'S4_realiz_scree_' + season[i] + '.png'
	ttle = 'Variance Explained by S4 realizations modes - ' + season[i]
	PlotScree(exp_var_s4_realiz, 36*51, title, filename)
	#compute correlation between observed and forecasted PCs
	eofs = np.reshape(eofs,[5, np.shape(lat_erai)[0] * np.shape(lat_erai)[1] ] )
	#eofs_s4 = np.reshape(eofs_s4,[10, np.shape(lat_s4)[0] * np.shape(lat_s4)[1] ] )
	pcc = np.empty([NMODES, 10])
	tcc = np.empty([NMODES, 10])
	for j in np.arange(NMODES): #loop over observed modes
		for k in np.arange(10): #loop over forecasted modes
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
	#compute correlation between observed and forecasted PCs
	pcc = np.empty([NMODES, 10])

	eofs = np.reshape(eofs,[5, np.shape(lat_erai)[0] * np.shape(lat_erai)[1] ] )
	for j in np.arange(NMODES): #loop over observed modes
		for k in np.arange(5): #loop over forecasted modes
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

#put all eofs and pcs together
s4_realiz_pc_sorted = np.reshape(s4_realiz_pc_sorted, [5, 36, 51, 3])

for i in range(5):
	title = 'S4 realizations HGT 200hPa EOFs - ' + season[i]
	filename = FIG_PATH + 'sorted_S4_realiz_eof_' + season[i] + '.png'
	PlotEOF(s4_realiz_eof_sorted[i, :, :, :], lat_s4, lon_s4, title, filename)
	title = 'S4 HGT 200hPa EOFs - ' + season[i]
	filename = FIG_PATH + 'sorted_S4_eof_' + season[i] + '.png'
	PlotEOF(s4_eof_sorted[i, :, :, :], lat_s4, lon_s4, title, filename)
	title = 'Observed and S4 HGT 200hPa PCs - ' + season[i]
	filename = FIG_PATH + 'PCs_' + season[i] + '.png'
	PlotPCs(pc_erai[i, :, :], s4_pc_sorted[i, :, :], s4_realiz_pc_sorted[i, :, :, :], np.arange(1981, 2017), title, filename)

