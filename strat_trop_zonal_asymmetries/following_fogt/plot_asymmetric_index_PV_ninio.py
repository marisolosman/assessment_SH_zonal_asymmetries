#plot asymmetric index for ninio events conditioned on PV
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_MONTHLY_ASYMMETRIES = 'monthly_asymmetric_index_SPoV_enso_polar.nc4'
asymmetric_index = xr.open_dataset(PATH_DATA_2 + FILE_MONTHLY_ASYMMETRIES)

years = np.concatenate([np.arange(1981, 2002), np.arange(2003, 2018)])

fig1 = plt.figure(figsize = (20, 21), dpi = 300)  #fig size in inches
filename = FIG_PATH + 'monthly_asymmetric_index_PV_LN-all.jpg'
for j in np.arange(7): #loop over Months
	indexall = np.reshape(asymmetric_index.asymm_all.values[j, :], [36, 51])
	indexEN = np.reshape(asymmetric_index.asymm_EN.values[j, :], [36, 51])
	indexLN = np.reshape(asymmetric_index.asymm_LN.values[j, :], [36, 51])
	ax = plt.subplot(4, 2, j + 1)
	ax.plot(years, np.mean(indexall, 1)/1e3, color = 'gray', linewidth=5, label = 'All Years')
	ax.plot(years, np.max(indexall, 1)/1e3, color = 'silver',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexall, 1)/1e3, color = 'silver',
		 linewidth=3, linestyle='--', alpha=0.55)

	#ax.plot(years, np.mean(indexEN, 1)/1e3, color = '#67a9cf', linewidth=5, label = 'EN')
	#ax.plot(years, np.max(indexEN, 1)/1e3, color = '#67a9cf',
	#	 linewidth=3, linestyle='--', alpha=0.55)
	#ax.plot(years, np.min(indexEN, 1)/1e3, color = '#67a9cf',
	#	 linewidth=3, linestyle='--', alpha=0.55)

	ax.plot(years, np.mean(indexLN, 1)/1e3, color = '#ef8a62', linewidth=5, label = 'LN')
	ax.plot(years, np.max(indexLN, 1)/1e3, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexLN, 1)/1e3, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)

	#ax.fill_between(years, np.min(indexEN, 1) / 1e3, np.max(indexEN, 1) / 1e3,
	#		facecolor='#67a9cf', alpha=0.55)
	ax.fill_between(years, np.min(indexLN, 1) / 1e3, np.max(indexLN, 1) / 1e3,
			facecolor='#ef8a62', alpha=0.55, linewidth=2)
	ax.fill_between(years, np.min(indexall, 1) / 1e3, np.max(indexall, 1) / 1e3,
			facecolor='silver', alpha=0.55)
	ax.set_xlim((years[0] - 1, years[-1] + 1))
	ax.set_ylim((-8.5, 8.5))
	plt.title(asymmetric_index.month.values[j], fontsize=10)
ax.legend(bbox_to_anchor=(1.5, 0.8), loc='lower left', borderaxespad=0, fontsize=14, ncol=2)
# it will place the legend on the outer right-hand side of the last axes
plt.suptitle('Asymmetric Index SPoV-WPoV' ,fontsize=12, x=0.52, y=0.93)
fig1.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')
plt.clf()
plt.cla()
plt.close()


FILE_SEASONAL_ASYMMETRIES = 'seasonal_asymmetric_index_SPoV_enso_polar.nc4'
asymmetric_index = xr.open_dataset(PATH_DATA_2 + FILE_SEASONAL_ASYMMETRIES)
fig1 = plt.figure(figsize = (20, 21), dpi = 300)  #fig size in inches
filename = FIG_PATH + 'seasonal_asymmetric_index_PV_LN-all.jpg'
for j in np.arange(5): #loop over Seasons
	indexall = np.reshape(asymmetric_index.asymm_all.values[j, :], [36, 51])
	indexEN = np.reshape(asymmetric_index.asymm_EN.values[j, :], [36, 51])
	indexLN = np.reshape(asymmetric_index.asymm_LN.values[j, :], [36, 51])
	ax = plt.subplot(4, 2, j + 1)
	ax.plot(years, np.mean(indexall, 1)/1e3, color = 'gray', linewidth=5, label = 'All years')
	ax.plot(years, np.max(indexall, 1)/1e3, color = 'silver',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexall, 1)/1e3, color = 'silver',
		 linewidth=3, linestyle='--', alpha=0.55)
	#ax.plot(years, np.mean(indexEN, 1)/1e3, color = '#67a9cf', linewidth=5, label = 'EN')
	#ax.plot(years, np.max(indexEN, 1)/1e3, color = '#67a9cf',
	#	 linewidth=3, linestyle='--', alpha=0.55)
	#ax.plot(years, np.min(indexEN, 1)/1e3, color = '#67a9cf',
	#	 linewidth=3, linestyle='--', alpha=0.55)

	ax.plot(years, np.mean(indexLN, 1)/1e3, color = '#ef8a62', linewidth=5, label = 'LN')
	ax.plot(years, np.max(indexLN, 1)/1e3, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexLN, 1)/1e3, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)

	ax.fill_between(years, np.min(indexall, 1) / 1e3, np.max(indexall, 1) / 1e3,
			facecolor='silver', alpha=0.55)
	#ax.fill_between(years, np.min(indexEN, 1) / 1e3, np.max(indexEN, 1) / 1e3,
	#		facecolor='#67a9cf', alpha=0.55)
	ax.fill_between(years, np.min(indexLN, 1) / 1e3, np.max(indexLN, 1) / 1e3,
			facecolor='#ef8a62', alpha=0.55, linewidth=2)

	ax.set_xlim((years[0] - 1, years[-1] + 1))
	ax.set_ylim((-7, 7))
	plt.title(asymmetric_index.seas.values[j], fontsize=10)
ax.legend(bbox_to_anchor=(1.5, 0.8), loc='lower left', borderaxespad=0, fontsize=14, ncol=2)
plt.suptitle('Asymmetric Index SPoV-WPoV' ,fontsize=12, x=0.52, y=0.93)
fig1.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')
plt.clf()
plt.cla()
plt.close()


