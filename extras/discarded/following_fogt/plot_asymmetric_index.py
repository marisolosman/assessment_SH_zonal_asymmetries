#plot asymmetric index for ninio events conditioned on PV
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_MONTHLY_ASYMMETRIES = 'monthly_asimmetric_index_enso_SPoV_polar.nc4'
asymmetric_index = xr.open_dataset(PATH_DATA_2 + FILE_MONTHLY_ASYMMETRIES)

years = np.concatenate([np.arange(1981, 2002), np.arange(2003, 2018)])

fig1 = plt.figure(figsize = (20, 21), dpi = 300)  #fig size in inches
filename = FIG_PATH + 'monthly_asymmetric_index_ninio_polar_strong_weak.jpg'
for j in np.arange(7): #loop over Months
	indexall = np.reshape(asymmetric_index.asimm_ninio_all.values[j, :], [36, 51])
	indexSPV = np.reshape(asymmetric_index.asimm_ninio_SPV.values[j, :], [36, 51])
	indexWPV = np.reshape(asymmetric_index.asimm_ninio_WPV.values[j, :], [36, 51])
	ax = plt.subplot(4, 2, j + 1)
	#ax.plot(years, np.mean(indexall, 1)/1e0, color = 'gray', linewidth=5, label = 'All PoV')
	#ax.plot(years, np.max(indexall, 1)/1e0, color = 'silver',
	#	 linewidth=3, linestyle='--', alpha=0.55)
	#ax.plot(years, np.min(indexall, 1)/1e0, color = 'silver',
	#	 linewidth=3, linestyle='--', alpha=0.55)

	ax.plot(years, np.mean(indexWPV, 1)/1e0, color = '#67a9cf', linewidth=5, label = 'Weak PoV')
	ax.plot(years, np.max(indexWPV, 1)/1e0, color = '#67a9cf',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexWPV, 1)/1e0, color = '#67a9cf',
		 linewidth=3, linestyle='--', alpha=0.55)

	ax.plot(years, np.mean(indexSPV, 1)/1e0, color = '#ef8a62', linewidth=5, label = 'Strong PoV')
	ax.plot(years, np.max(indexSPV, 1)/1e0, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexSPV, 1)/1e0, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)

	ax.fill_between(years, np.min(indexWPV, 1) / 1e0, np.max(indexWPV, 1) / 1e0,
			facecolor='#67a9cf', alpha=0.55)
	ax.fill_between(years, np.min(indexSPV, 1) / 1e0, np.max(indexSPV, 1) / 1e0,
			facecolor='#ef8a62', alpha=0.55, linewidth=2)
	#ax.fill_between(years, np.min(indexall, 1) / 1e0, np.max(indexall, 1) / 1e0,
	#		facecolor='silver', alpha=0.55)
	ax.set_xlim((years[0] - 1, years[-1] + 1))
	ax.set_ylim((-2, 2))
	plt.title(asymmetric_index.month.values[j], fontsize=10)
ax.legend(bbox_to_anchor=(1.5, 0.8), loc='lower left', borderaxespad=0, fontsize=14, ncol=2)
# it will place the legend on the outer right-hand side of the last axes
plt.suptitle('Asymmetric Index Ninio' ,fontsize=12, x=0.52, y=0.93)
fig1.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')
plt.clf()
plt.cla()
plt.close()

fig1 = plt.figure(figsize = (20, 21), dpi = 300)  #fig size in inches
filename = FIG_PATH + 'monthly_asymmetric_index_ninia_polar_all_strong.jpg'
for j in np.arange(7): #loop over Months
	indexSPV = np.reshape(asymmetric_index.asimm_ninia_SPV.values[j, :], [36, 51])
	indexWPV = np.reshape(asymmetric_index.asimm_ninia_WPV.values[j, :], [36, 51])
	indexall = np.reshape(asymmetric_index.asimm_ninia_all.values[j, :], [36, 51])

	ax = plt.subplot(4, 2, j + 1)
	ax.plot(years, np.mean(indexall, 1)/1e0, color = 'gray', linewidth=5, label = 'All PoV')
	ax.plot(years, np.max(indexall, 1)/1e0, color = 'silver',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexall, 1)/1e0, color = 'silver',
		 linewidth=3, linestyle='--', alpha=0.55)

	#ax.plot(years, np.mean(indexWPV, 1)/1e0, color = '#67a9cf', linewidth=5, label = 'Weak PoV')
	#ax.plot(years, np.max(indexWPV, 1)/1e0, color = '#67a9cf',
	#	 linewidth=3, linestyle='--', alpha=0.55)
	#ax.plot(years, np.min(indexWPV, 1)/1e0, color = '#67a9cf',
	#	 linewidth=3, linestyle='--', alpha=0.55)

	ax.plot(years, np.mean(indexSPV, 1)/1e0, color = '#ef8a62', linewidth=5, label = 'Strong PoV')
	ax.plot(years, np.max(indexSPV, 1)/1e0, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexSPV, 1)/1e0, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)

	ax.fill_between(years, np.min(indexall, 1) / 1e0, np.max(indexall, 1) / 1e0,
			facecolor='silver', alpha=0.55)
	#ax.fill_between(years, np.min(indexWPV, 1) / 1e0, np.max(indexWPV, 1) / 1e0,
	#		facecolor='#67a9cf', alpha=0.55)
	ax.fill_between(years, np.min(indexSPV, 1) / 1e0, np.max(indexSPV, 1) / 1e0,
			facecolor='#ef8a62', alpha=0.55, linewidth=2)

	ax.set_xlim((years[0] - 1, years[-1] + 1))
	ax.set_ylim((-2, 2))
	plt.title(asymmetric_index.month.values[j], fontsize=10)
ax.legend(bbox_to_anchor=(1.5, 0.8), loc='lower left', borderaxespad=0, fontsize=14, ncol=2)
plt.suptitle('Asymmetric Index Ninia' ,fontsize=12, x=0.52, y=0.93)
fig1.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')
plt.clf()
plt.cla()
plt.close()

for i in np.arange(7):
	index_ninia_SPV = np.reshape(asymmetric_index.asimm_ninia_SPV.values[j, :], [36, 51])
	index_ninia_WPV = np.reshape(asymmetric_index.asimm_ninia_WPV.values[j, :], [36, 51])
	index_ninio_SPV = np.reshape(asymmetric_index.asimm_ninio_SPV.values[j, :], [36, 51])
	index_ninio_WPV = np.reshape(asymmetric_index.asimm_ninio_WPV.values[j, :], [36, 51])

	fig1 = plt.figure(figsize = (20, 21), dpi = 300)  #fig size in inches
	filename = FIG_PATH + 'monthly_asymmetric_index_' + asymmetric_index.month.values[i] +  '_polar.jpg'
	ax = plt.subplot(2, 3, 1)
	ax.plot(years, np.mean(index_ninio_WPV, 1)/1, color = 'k', linewidth=1.5, label = 'Weak PoV')
	ax.fill_between(years, np.min(index_ninio_WPV, 1) / 1, np.max(index_ninio_WPV, 1) / 1,
			facecolor='darkgray', alpha=0.35)
	ax.set_xlim((years[0] - 1, years[-1] + 1))
	plt.title("E Years - Weak PoV")
	ax = plt.subplot(2, 3, 2)
	ax.plot(years, np.mean(index_ninio_SPV, 1)/1, color = 'k', linewidth=1.5, label = 'Weak PoV')
	ax.fill_between(years, np.min(index_ninio_SPV, 1) / 1, np.max(index_ninio_SPV, 1) / 1,
			facecolor='darkgray', alpha=0.35)
	ax.set_xlim((years[0] - 1, years[-1] + 1))
	plt.title("E Years - Strong PoV")
	ax = plt.subplot(2, 3, 3)
	heatmap = (index_ninio_WPV - index_ninio_SPV) / 1
	im = plt.imshow(heatmap, cmap=plt.cm.RdBu, interpolation='bilinear')
	plt.colorbar(im)

	ax = plt.subplot(2, 3, 4)
	ax.plot(years, np.mean(index_ninio_SPV, 1)/1e0, color = 'k', linewidth=1.5, label = 'Strong PoV')
	ax.fill_between(years, np.min(index_ninio_SPV, 1) / 1e0, np.max(index_ninio_SPV, 1) / 1,
			facecolor='darkgray', alpha=0.35)
	ax.set_xlim((years[0] - 1, years[-1] + 1))
	plt.title("E Years - Strong PoV")
	ax = plt.subplot(2, 3, 5)
	ax.plot(years, np.mean(index_ninia_SPV, 1)/1e0, color = 'k', linewidth=1.5, label = 'Strong PoV')
	ax.fill_between(years, np.min(index_ninia_SPV, 1) / 1e0, np.max(index_ninia_SPV, 1) / 1,
			facecolor='darkgray', alpha=0.35)
	ax.set_xlim((years[0] - 1, years[-1] + 1))
	plt.title("L Years - Strong PoV")
	ax = plt.subplot(2, 3, 6)
	heatmap = index_ninio_SPV - index_ninia_SPV
	im = plt.imshow(heatmap, cmap=plt.cm.RdBu, interpolation='bilinear')
	plt.colorbar(im)
	plt.suptitle(asymmetric_index.month.values[j] ,fontsize=12, x=0.52, y=0.93)
#	fig1.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')
	plt.clf()
	plt.cla()
plt.close()




FILE_SEASONAL_ASYMMETRIES = 'seasonal_asimmetric_index_enso_SPoV_polar.nc4'
asymmetric_index = xr.open_dataset(PATH_DATA_2 + FILE_SEASONAL_ASYMMETRIES)
fig1 = plt.figure(figsize = (20, 21), dpi = 300)  #fig size in inches
filename = FIG_PATH + 'seasonal_asymmetric_index_ninio_polar_strong_weak.jpg'
for j in np.arange(5): #loop over Seasons
	indexall = np.reshape(asymmetric_index.asimm_ninio_all.values[j, :], [36, 51])
	indexSPV = np.reshape(asymmetric_index.asimm_ninio_SPV.values[j, :], [36, 51])
	indexWPV = np.reshape(asymmetric_index.asimm_ninio_WPV.values[j, :], [36, 51])
	ax = plt.subplot(4, 2, j + 1)
	#ax.plot(years, np.mean(indexall, 1)/1e0, color = 'gray', linewidth=5, label = 'All PoV')
	#ax.plot(years, np.max(indexall, 1)/1e0, color = 'silver',
	#	 linewidth=3, linestyle='--', alpha=0.55)
	#ax.plot(years, np.min(indexall, 1)/1e0, color = 'silver',
	#	 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.mean(indexWPV, 1)/1e0, color = '#67a9cf', linewidth=5, label = 'Weak PoV')
	ax.plot(years, np.max(indexWPV, 1)/1e0, color = '#67a9cf',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexWPV, 1)/1e0, color = '#67a9cf',
		 linewidth=3, linestyle='--', alpha=0.55)

	ax.plot(years, np.mean(indexSPV, 1)/1e0, color = '#ef8a62', linewidth=5, label = 'Strong PoV')
	ax.plot(years, np.max(indexSPV, 1)/1e0, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexSPV, 1)/1e0, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)

	#ax.fill_between(years, np.min(indexall, 1) / 1e0, np.max(indexall, 1) / 1e0,
	#		facecolor='silver', alpha=0.55)
	ax.fill_between(years, np.min(indexWPV, 1) / 1e0, np.max(indexWPV, 1) / 1e0,
			facecolor='#67a9cf', alpha=0.55)
	ax.fill_between(years, np.min(indexSPV, 1) / 1e0, np.max(indexSPV, 1) / 1e0,
			facecolor='#ef8a62', alpha=0.55, linewidth=2)

	ax.set_xlim((years[0] - 1, years[-1] + 1))
	ax.set_ylim((-2, 2))
	plt.title(asymmetric_index.seas.values[j], fontsize=10)
ax.legend(bbox_to_anchor=(1.5, 0.8), loc='lower left', borderaxespad=0, fontsize=14, ncol=2)
plt.suptitle('Asymmetric Index Ninio' ,fontsize=12, x=0.52, y=0.93)
fig1.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')
plt.clf()
plt.cla()
plt.close()



fig1 = plt.figure(figsize = (20, 21), dpi = 300)  #fig size in inches
filename = FIG_PATH + 'seasonal_asymmetric_index_ninia_polar_strong_weak.jpg'
for j in np.arange(5): #loop over Seasons
	indexall = np.reshape(asymmetric_index.asimm_ninia_all.values[j, :], [36, 51])
	indexSPV = np.reshape(asymmetric_index.asimm_ninia_SPV.values[j, :], [36, 51])
	indexWPV = np.reshape(asymmetric_index.asimm_ninia_WPV.values[j, :], [36, 51])
	ax = plt.subplot(4, 2, j + 1)
	#ax.plot(years, np.mean(indexall, 1)/1e0, color = 'gray', linewidth=5, label = 'All PoV')
	#ax.plot(years, np.max(indexall, 1)/1e0, color = 'silver',
	#	 linewidth=3, linestyle='--', alpha=0.55)
	#ax.plot(years, np.min(indexall, 1)/1e0, color = 'silver',
	#	 linewidth=3, linestyle='--', alpha=0.55)

	ax.plot(years, np.mean(indexSPV, 1)/1e0, color = '#ef8a62', linewidth=5, label = 'Strong PoV')
	ax.plot(years, np.max(indexSPV, 1)/1e0, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexSPV, 1)/1e0, color = '#ef8a62',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.mean(indexWPV, 1)/1e0, color = '#67a9cf', linewidth=5, label = 'Weak PoV')
	ax.plot(years, np.max(indexWPV, 1)/1e0, color = '#67a9cf',
		 linewidth=3, linestyle='--', alpha=0.55)
	ax.plot(years, np.min(indexWPV, 1)/1e0, color = '#67a9cf',
		 linewidth=3, linestyle='--', alpha=0.55)

	#ax.fill_between(years, np.min(indexall, 1) / 1e0, np.max(indexall, 1) / 1e0,
	#		facecolor='silver', alpha=0.55)
	ax.fill_between(years, np.min(indexWPV, 1) / 1e0, np.max(indexWPV, 1) / 1e0,
			facecolor='#67a9cf', alpha=0.55)
	ax.fill_between(years, np.min(indexSPV, 1) / 1e0, np.max(indexSPV, 1) / 1e0,
			facecolor='#ef8a62', alpha=0.55, linewidth=2)
	ax.set_xlim((years[0] - 1, years[-1] + 1))
	ax.set_ylim((-2, 2))
	plt.title(asymmetric_index.seas.values[j], fontsize=10)
ax.legend(bbox_to_anchor=(1.5, 0.8), loc='lower left', borderaxespad=0, fontsize=14, ncol=2)
plt.suptitle('Asymmetric Index Ninia' ,fontsize=12, x=0.52, y=0.93)
fig1.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')
plt.clf()
plt.cla()
plt.close()

