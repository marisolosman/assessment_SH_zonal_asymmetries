import numpy as np
import datetime
import pandas as pd
import xarray as xr
from eofs.standard import Eof
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '~/datos/data/' #este spanglish no te lo robo amiga
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/model_assessment/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_PV_ERAI = 'fogt/SPV_index_erai.nc4'
SPV_erai =  xr.open_dataset(RUTA + FILE_PV_ERAI)
hgt_s4 = xr.open_dataset(RUTA + FILE_HGT_S4)
hgt_s4 = hgt_s4.sel(**{'latitude':slice(-20, -90)})
FILE_HGT_ERAI = 'hgt_erai_200.nc4'
hgt_erai = xr.open_dataset(RUTA + FILE_HGT_ERAI)
hgt_erai.time.values = hgt_erai.valid_time.values
#remove data from august 2002 to february 2003
hgt_erai = hgt_erai.sel(time=hgt_erai.time.values[np.logical_or(hgt_erai.time.values<=np.datetime64('2002-07-31'), hgt_erai.time.values>=np.datetime64('2003-03-01'))])
hgt_erai = hgt_erai.sel(**{'time':slice('1981-08-01', '2018-02-01')})
hgt_erai = hgt_erai.sel(**{'latitude':slice(-20, -90)})

#select SSPV WSPV and NSPV dates
index_SPV_erai_upper = SPV_erai.SPV_mon >= SPV_erai.SPV_mon.quantile(0.75, dim='dim_0')
index_SPV_erai_lower = SPV_erai.SPV_mon <= SPV_erai.SPV_mon.quantile(0.25, dim='dim_0')
index_SPV_erai_normal = np.logical_and(SPV_erai.SPV_mon < SPV_erai.SPV_mon.quantile(0.75, dim='dim_0'), SPV_erai.SPV_mon > SPV_erai.SPV_mon.quantile(0.25, dim='dim_0'))

index_SSPV_erai_upper = SPV_erai.SPV_mon >= SPV_erai.SPV_mon.quantile(0.90, dim='dim_0')
index_SSPV_erai_lower = SPV_erai.SPV_mon <= SPV_erai.SPV_mon.quantile(0.10, dim='dim_0')
index_SSPV_erai_normal = np.logical_and(SPV_erai.SPV_mon < SPV_erai.SPV_mon.quantile(0.90, dim='dim_0'), SPV_erai.SPV_mon > SPV_erai.SPV_mon.quantile(0.10, dim='dim_0'))

#idem s4
SPV_s4 =  xr.open_dataset(RUTA + FILE_PV_S4)
index_SPV_s4_upper = SPV_s4.SPV_index >= SPV_s4.SPV_index.quantile(0.75, dim='dim_0')
index_SPV_s4_lower = SPV_s4.SPV_index <= SPV_s4.SPV_index.quantile(0.25, dim='dim_0')
index_SPV_s4_normal = np.logical_and(SPV_s4.SPV_index < SPV_s4.SPV_index.quantile(0.75, dim='dim_0'), SPV_s4.SPV_index > SPV_s4.SPV_index.quantile(0.25, dim='dim_0'))

index_SSPV_s4_upper = SPV_s4.SPV_index >= SPV_s4.SPV_index.quantile(0.90, dim='dim_0')
index_SSPV_s4_lower = SPV_s4.SPV_index <= SPV_s4.SPV_index.quantile(0.10, dim='dim_0')
index_SSPV_s4_normal = np.logical_and(SPV_s4.SPV_index < SPV_s4.SPV_index.quantile(0.90, dim='dim_0'), SPV_s4.SPV_index > SPV_s4.SPV_index.quantile(0.10, dim='dim_0'))

season = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lmonth = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec']
#loop over seasons, selec data and performs boxplot
SAM_erai = np.zeros([5, index_SPV_erai_upper.values.shape[0]])
SAM_s4 = np.zeros([5, index_SPV_s4_upper.values.shape[0]])
#should be the other way around
SAM_erai_SSPV = np.zeros([5, index_SPV_erai_upper.values.sum()])
SAM_erai_sSSPV = np.zeros([5, index_SSPV_erai_upper.values.sum()])
SAM_erai_WSPV = np.zeros([5, index_SPV_erai_lower.values.sum()])
SAM_erai_sWSPV = np.zeros([5, index_SSPV_erai_lower.values.sum()])
SAM_s4_SSPV = np.zeros([5, index_SPV_s4_upper.values.sum()])
SAM_s4_sSSPV = np.zeros([5, index_SSPV_s4_upper.values.sum()])
SAM_s4_WSPV = np.zeros([5, index_SPV_s4_lower.values.sum()])
SAM_s4_sWSPV = np.zeros([5, index_SSPV_s4_lower.values.sum()])
eof_erai = np.zeros([5, hgt_erai.latitude.values.shape[0]])
eof_s4 = np.zeros([5, hgt_s4.latitude.values.shape[0]])
sign_s4 = np.array([1, 1, 1, 1, -1])
sign_erai = np.array([1, -1, -1, -1, 1])
for i in np.arange(0, 5):
	aux = hgt_erai['z'].resample(time='QS-' + lmonth[i]).mean(dim='time',skipna=True)
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	aux = aux.sel(time= np.logical_and(aux['time.month'] == mes, aux['time.year']!=2002))
	X_zm = aux.mean(dim='longitude')
	X_an = X_zm - X_zm.mean(dim='time')
	solver = Eof(X_an.values)
	pcs = solver.pcs(npcs=1, pcscaling=1)
	eof_erai[i, :] = solver.eofs(neofs=1)[0,:]
	SAM_erai [i, :] = sign_erai[i] * pcs[:, 0]
	SAM_erai_SSPV[i, :] = sign_erai[i] * pcs[index_SPV_erai_upper, 0]
	SAM_erai_WSPV[i, :] = sign_erai[i] * pcs[index_SPV_erai_lower, 0] 
	SAM_erai_sSSPV[i, :] = sign_erai[i] * pcs[index_SSPV_erai_upper, 0]
	SAM_erai_sWSPV[i, :] = sign_erai[i] * pcs[index_SSPV_erai_lower, 0]
	hgt_s4_smean = np.nanmean(np.nanmean(hgt_s4.z.values[i:i + 3, :, :, :], axis=0), axis=2)
	hgt_s4_smean = hgt_s4_smean - np.nanmean(hgt_s4_smean, axis=0)
        # Create an EOF solver to do the EOF analysis. Square-root of cosine of
        # latitude weights are applied before the computation of EOFs.
	solver = Eof(hgt_s4_smean)#, weights=wgts)
	pcs = solver.pcs(npcs=1, pcscaling=1)
	eof_s4[i, :] = solver.eofs(neofs=1)[0,:]
	SAM_s4 [i, :] = sign_s4[i] * pcs[:, 0]
	SAM_s4_SSPV[i, :] = sign_s4[i] * pcs[index_SPV_s4_upper, 0]
	SAM_s4_WSPV[i, :] = sign_s4[i] * pcs[index_SPV_s4_lower, 0]
	SAM_s4_sSSPV[i, :] = sign_s4[i] * pcs[index_SSPV_s4_upper, 0]
	SAM_s4_sWSPV[i, :] = sign_s4[i] * pcs[index_SSPV_s4_lower, 0]

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
# plot violin plot
axes[0].boxplot(SAM_erai_SSPV.T, positions=np.arange(1.20, 8.8,1.6))
axes[0].boxplot(SAM_erai_sSSPV.T, positions=np.arange(1.8,9.3,1.6))

axes[0].set_title('Zonal mean SAM index at 200hPa - ERAI')

# plot box plot
axes[1].violinplot(SAM_s4_SSPV.T,positions=np.arange(1.20, 8.8, 1.6),
		   showmeans=False,
		   showmedians=True)
axes[1].violinplot(SAM_s4_sSSPV.T, positions=np.arange(1.8, 9.3, 1.6),
		   showmeans=False, showmedians=True)

axes[1].set_title('Zonal mean SAM index at 200hPa - S4')
for ax in axes:
	# adding horizontal grid lines
	ax.yaxis.grid(True)
	#ax.set_xticks(range(5))
	#ax.set_xtickslabel(season)
	ax.set_xlim([0.5, 9])
	ax.set_ylim([-4.8, 4.8])
	ax.set_xlabel('Seasons')
	ax.set_ylabel('SAM index')
# add x-tick labels
plt.setp(axes, xticks=np.arange(1.5, 9, 1.6), xticklabels=season)

# add x-tick labels
plt.savefig(FIG_PATH + 'boxplot_SAM_WSPV.png', res=400)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
# plot violin plot
axes[0].boxplot(SAM_erai.T, positions=np.arange(1.20, 8.8,1.6))
axes[0].boxplot(SAM_erai_SSPV.T, positions=np.arange(1.8,9.3,1.6))

axes[0].set_title('Zonal mean SAM index at 200hPa - ERAI')

# plot box plot
axes[1].violinplot(SAM_s4.T,positions=np.arange(1.20, 8.8, 1.6),
		   showmeans=False,
		   showmedians=True)
axes[1].violinplot(SAM_s4_SSPV.T, positions=np.arange(1.8, 9.3, 1.6),
		   showmeans=False, showmedians=True)

axes[1].set_title('Zonal mean SAM index at 200hPa - S4')
for ax in axes:
	# adding horizontal grid lines
	ax.yaxis.grid(True)
	#ax.set_xticks(range(5))
	#ax.set_xtickslabel(season)
	ax.set_xlim([0.5, 9])
	ax.set_ylim([-4.8, 4.8])
	ax.set_xlabel('Seasons')
	ax.set_ylabel('SAM index')
# add x-tick labels
plt.setp(axes, xticks=np.arange(1.5, 9, 1.6), xticklabels=season)

# add x-tick labels
plt.savefig(FIG_PATH + 'boxplot_SAM_SPV.png', res=400)
limite = np.max(np.abs(SPV_erai.SPV_mon.values))

fig = plt.figure(figsize=(12, 12), dpi=300)
for i in range(5):
	results = sm.OLS(SAM_erai[i, :], sm.add_constant(SPV_erai.SPV_mon.values)).fit()
	x_pred = np.linspace(-limite, limite, 1000)
	y_pred = x_pred * results.params[1] + results.params[0]
	print(results.summary())
	x_pred2 = sm.add_constant(x_pred)
	sdev, lower, upper = wls_prediction_std(results, exog=x_pred2, alpha=0.05)
	interv = results.conf_int(alpha=0.05)
	plt.subplot(3, 2, i + 1)
	plt.scatter(SPV_erai.SPV_mon.values, SAM_erai[i,:])
	plt.plot(x_pred, y_pred)
	plt.fill_between(x_pred, lower, upper, color='#888888', alpha=0.1)
	#plt.xlim([-6, 6])
	plt.ylim([-4, 4])
	plt.xlabel('PoV Index')
	plt.ylabel('Zonal SAM index')
	plt.title(season[i] + '- correlation: ' + '{:03.2f}'.format(np.corrcoef(SPV_erai.SPV_mon.values, SAM_erai[i, :]) [0,1]))
plt.suptitle('Regression SAM vs PoV - ERAI')
plt.subplots_adjust(wspace=0.3, hspace=0.30)
plt.savefig(FIG_PATH + 'regression_sam_pov.png', res=400, box_inches='tight'
)
limite = np.max(np.abs(SPV_s4.SPV_index.values))
#scatter plot between PoV and SAM
fig = plt.figure(figsize=(14, 12), dpi=300)
for i in range(5):
	results = sm.OLS(SAM_s4[i, :], sm.add_constant(SPV_s4.SPV_index.values)).fit()
	x_pred = np.linspace(-limite, limite, 1000)
	y_pred = x_pred * results.params[1] + results.params[0]
	print(results.summary())
	x_pred2 = sm.add_constant(x_pred)
	sdev, lower, upper = wls_prediction_std(results, exog=x_pred2, alpha=0.05)
	interv = results.conf_int(alpha=0.05)
	plt.subplot(3, 2, i + 1)
	plt.scatter(SPV_s4.SPV_index.values, SAM_s4[i,:])
	plt.plot(x_pred, y_pred)
	plt.fill_between(x_pred, lower, upper, color='#888888', alpha=0.1)
	plt.xlim([-limite, limite])
	plt.ylim([-4, 4])
	plt.xlabel('PoV Index')
	plt.ylabel('Zonal SAM index')
	plt.title(season[i] + '- correlation: ' + '{:03.2f}'.format(np.corrcoef(SPV_s4.SPV_index.values, SAM_s4[i, :])[0,1]))
plt.subplots_adjust(wspace=0.3, hspace=0.30)
plt.suptitle('Regression SAM vs PoV - S4')
plt.savefig(FIG_PATH + 'regression_sam_pov_s4.png', res=400,bbox_inches='tight' )


