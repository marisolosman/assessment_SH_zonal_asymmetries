import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.path as mpath
from matplotlib.transforms import offset_copy
def PlotCompositesdivPF(var, lat, lon, title, filename):
	var_keys = [*var]
	tit = ['Weak SPoV (All ENSO)', 'Strong SpoV (All ENSO)',
		'Ninio (All PoV)', 'Ninia (All PoV)',
		'Ninio (Weak PoV)', 'Ninia (Weak PoV)',
		'Ninio (Strong PoV)', 'Ninia (Strong PoV)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 8.9), 300)
	clevs = np.arange(-0.75, 1.25, 0.5)
	barra = plt.cm.RdBu_r
	for i in range(8):
		ax = plt.subplot(4, 2, i + 1, projection=proj)
		ax.set_extent([0, 359, -90, -18], crs=ccrs.PlateCarree())
		lat = lat[lat < 0]
		lat1 = lat[lat < -20]
		M = var[var_keys[i * 3]] * 1e6
		M = M[lat < -20, :]
		im = ax.contourf(lon, lat1, M, clevs, transform=ccrs.PlateCarree(),
				 vmin=clevs[0], vmax=clevs[-1], cmap=barra, alpha=0.5,
				 extend='both')
		lat1 = lat[lat < 0]
		u = var[var_keys[i * 3 + 1]][lat1 < -30, :]
		v = var[var_keys[i * 3 + 2]][lat1 < -30, :]
		lat1 = lat1[lat1 < -30]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.0007
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		im2 = ax.quiver(lon[2:-1:10], lat1[2:-1:8], u_mask[2:-1:8, 2:-1:10],
			  v_mask[2:-1:8, 2:-1:10], transform=ccrs.PlateCarree(),
			  angles='xy', scale=10)
		#barra.set_under(barra(0))
		#barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	#fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.25, top=0.85, hspace=0.1, wspace=0.05)
	plt.quiverkey(im2, 0, 1.3, 0.5, "0.5 m2/s2", coordinates='axes', color='k')
	#fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotEnsoCompositesPoVU(var, lat, lon, title, filename):
	var_keys = [*var]
	tit = ['Ninio (All PoV)', 'Ninia (All PoV)',
		'Ninio (Weak PoV)', 'Ninia (Weak PoV)',
		'Ninio (Strong PoV)', 'Ninia (Strong PoV)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-7, 8, 2)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, var[var_keys[i]], clevs, transform=ccrs.PlateCarree(),
				cmap=barra, extend='both', vmin=-7, vmax=7)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
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

def PlotCompositesRWSChiWENSOPV(var_1, var_2, var_3, var_4, var_5, var_6, var_7, var_8,
			       u_1, u_2, u_3, u_4, u_5, u_6, u_7, u_8, v_1, v_2, v_3, v_4,
			       v_5, v_6, v_7, v_8, uu1, uu2, uu3, uu4, uu5, uu6, uu7, uu8,
			       lat, lon, title, filename):
	var = [var_1, var_5, var_2, var_6, var_3, var_7, var_4, var_8]
	tit = ['Weak-Strong SPV cond Ninio', 'Weak-Strong SPV cond Ninia',
		   'Ninio-Ninia cond Weak SPV', 'Ninio-Ninia cond Strong SPV',
		   'Ninio & Weak PoV', 'Ninia & Weak PoV',
		   'Ninio & Strong PoV', 'Ninia & Strong PoV']
	uchi = [u_1, u_5, u_2, u_6, u_3, u_7, u_4, u_8]
	vchi = [v_1, v_5, v_2, v_6, v_3, v_7, v_4, v_8]
	uu = [uu1, uu2, uu3, uu4, uu5, uu6, uu7, uu8]
	clevs = np.arange(-1.2, 1.6, 0.8)
	clevs2 = [30, 40]
	barra = plt.cm.RdBu_r
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 7), 300)
	for i in range(8):
		ax = plt.subplot(4, 2, i + 1, projection=proj)
		ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, var[i] / 10, clevs,
				transform=ccrs.PlateCarree(), cmap=barra, extend='both',
				vmin=-1.2, vmax=1.2)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		M = np.sqrt(np.add(np.power(uchi[i], 2), np.power(vchi[i], 2))) < 0.3
		#mask array
		u_mask = ma.array(uchi[i], mask = M)
		v_mask = ma.array(vchi[i], mask = M)
		im2 = ax.quiver(lon[0:-1:10], lat[0:-1:8], u_mask[0:-1:8, 0:-1:10],
				v_mask[0:-1:8, 0:-1:10], transform=ccrs.PlateCarree(),
				angles='xy', scale=40)
		ax.contour(lon, lat, uu[i], clevs2, linewidths=0.8, colors='green',
			   transform=ccrs.PlateCarree())
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.15, wspace=0.05)
	cbar_ax = fig.add_axes([0.34, 0.1, 0.25, 0.05])
	plt.quiverkey(im2, 0.5, -0.8, 1, "1 m/s", coordinates='axes', color='k')
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotCompositesdivPF(var50, var200, lat, lon, title, filename):
	var_keys = [*var50]
#	tit = ['Weak-Strong SPV \n cond Ninio', 'Weak-Strong SPV \n cond Ninia',
#		   'Ninio-Ninia \n cond Weak SPV', 'Ninio-Ninia \n cond Strong SPV',
#		   'Ninio & \n Weak SPV', 'Ninia & \n Strong SPV',
#		   'Ninio & \n Strong SPV', 'Ninia & \n Weak SPV']
	tit = ['Weak SPV', 'Strong SPV', 'Ninio', 'Ninia',
		   'Ninio & \n Weak SPV', 'Ninia & \n Strong SPV',
		   'Ninio & \n Strong SPV', 'Ninia & \n Weak SPV']
	lev = ['50 hPa', '200 hPa']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig, ax = plt.subplots(nrows=8, ncols=2, figsize=(5, 8), subplot_kw={'projection': proj})
	clevs = np.arange(-1.25, 1.75, 0.5)
	barra = plt.cm.RdBu_r
	for i in range(8):
		#ax = plt.subplot(8, 2, i * 2+ 1, projection=proj)
		ax[i, 0].set_extent([0, 359, -90, -18], crs=ccrs.PlateCarree())
		lat = lat[lat < 0]
		lat1 = lat[lat < -20]
		M = var50[var_keys[i * 3]] * 1e6
		M = M[lat < -20, :]
		ax[i, 0].contourf(lon, lat1, M, clevs, transform=ccrs.PlateCarree(),
				 vmin=clevs[0], vmax=clevs[-1], cmap=barra, alpha=0.5,
				 extend='both')
		lat1 = lat[lat < 0]
		u = var50[var_keys[i * 3 + 1]][lat1 < -30, :]
		v = var50[var_keys[i * 3 + 2]][lat1 < -30, :]
		lat1 = lat1[lat1 < -30]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.0007
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		ax[i, 0].quiver(lon[2:-1:10], lat1[2:-1:8], u_mask[2:-1:8, 2:-1:10],
			  v_mask[2:-1:8, 2:-1:10], transform=ccrs.PlateCarree(),
			  angles='xy', scale=25)
		ax[i, 0].coastlines(linewidth=0.4)
		#ax[i, 0].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax[i, 0].gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax[i, 0].xaxis.set_major_formatter(lon_formatter)
		ax[i, 0].yaxis.set_major_formatter(lat_formatter)

		ax[i, 1].set_extent([0, 359, -90, -18], crs=ccrs.PlateCarree())
		lat = lat[lat < 0]
		lat1 = lat[lat < -20]
		M = var200[var_keys[i * 3]] * 1e6
		M = M[lat < -20, :]
		im = ax[i, 1].contourf(lon, lat1, M, clevs, transform=ccrs.PlateCarree(),
				 vmin=clevs[0], vmax=clevs[-1], cmap=barra, alpha=0.5,
				 extend='both')
		lat1 = lat[lat < 0]
		u = var200[var_keys[i * 3 + 1]][lat1 < -30, :]
		v = var200[var_keys[i * 3 + 2]][lat1 < -30, :]
		lat1 = lat1[lat1 < -30]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.0007
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		im2 = ax[i, 1].quiver(lon[2:-1:10], lat1[2:-1:8], u_mask[2:-1:8, 2:-1:10],
			  v_mask[2:-1:8, 2:-1:10], transform=ccrs.PlateCarree(),
			  angles='xy', scale=25)
		ax[i, 1].coastlines(linewidth=0.4)
		#ax[i, 1].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax[i, 1].gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax[i, 1].xaxis.set_major_formatter(lon_formatter)
		ax[i, 1].yaxis.set_major_formatter(lat_formatter)
		if i==0:
			plt.quiverkey(im2, 0.9, 1.1, 0.5, r'$0.5 m^2/s^2$', labelpos='N',labelsep=0.03,
				      coordinates='axes', fontproperties={'size':7}, color='k')

	pad = 5
	for bx, col in zip(ax[0], lev):
		bx.annotate(col, xy=(0.5, 1), xytext=(0, pad),
			    xycoords='axes fraction', textcoords='offset points',
			    size='large', ha='center', va='baseline')
	for bx, row in zip(ax[:, 0], tit):
		bx.annotate(row, xy=(0, 0.5), xytext=(-bx.yaxis.labelpad-pad, 0),
			    xycoords='axes fraction', rotation=90, textcoords='offset points',
			    size='large', ha='center', va='center', fontsize=6)
	fig.subplots_adjust(top=0.85, hspace=0.35)#, wspace=0.03)
	plt.suptitle(title, fontsize=10, x=0.49, y=0.92)
	plt.savefig(filename, dpi=500, bbox_inches='tight', papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotCompositesdivPFx4(var50, var200, lat, lon, title, filename):
	var_keys = [*var50]
#	tit = ['Weak-Strong SPV \n cond Ninio', 'Weak-Strong SPV \n cond Ninia',
#		   'Ninio-Ninia \n cond Weak SPV', 'Ninio-Ninia \n cond Strong SPV',
#		   'Ninio & \n Weak SPV', 'Ninia & \n Strong SPV',
#		   'Ninio & \n Strong SPV', 'Ninia & \n Weak SPV']
	tit = ['Ninio & \n Weak SPV', 'Ninia & \n Strong SPV',
		   'Ninio & \n Strong SPV', 'Ninia & \n Weak SPV']
	lev = ['50 hPa', '200 hPa']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(5, 3.6), subplot_kw={'projection': proj})
	clevs = np.arange(-1.25, 1.75, 0.5)
	barra = plt.cm.RdBu_r
	for i in range(4):
		#ax = plt.subplot(8, 2, i * 2+ 1, projection=proj)
		ax[i, 0].set_extent([0, 359, -90, -18], crs=ccrs.PlateCarree())
		lat = lat[lat < 0]
		lat1 = lat[lat < -20]
		M = var50[var_keys[i * 3]] * 1e6
		M = M[lat < -20, :]
		ax[i, 0].contourf(lon, lat1, M, clevs, transform=ccrs.PlateCarree(),
				 vmin=clevs[0], vmax=clevs[-1], cmap=barra, alpha=0.5,
				 extend='both')
		lat1 = lat[lat < 0]
		u = var50[var_keys[i * 3 + 1]][lat1 < -30, :]
		v = var50[var_keys[i * 3 + 2]][lat1 < -30, :]
		lat1 = lat1[lat1 < -30]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.0007
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		ax[i, 0].quiver(lon[2:-1:10], lat1[2:-1:8], u_mask[2:-1:8, 2:-1:10],
			  v_mask[2:-1:8, 2:-1:10], transform=ccrs.PlateCarree(),
			  angles='xy', scale=25)
		ax[i, 0].coastlines(linewidth=0.4)
		#ax[i, 0].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax[i, 0].gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax[i, 0].xaxis.set_major_formatter(lon_formatter)
		ax[i, 0].yaxis.set_major_formatter(lat_formatter)

		ax[i, 1].set_extent([0, 359, -90, -18], crs=ccrs.PlateCarree())
		lat = lat[lat < 0]
		lat1 = lat[lat < -20]
		M = var200[var_keys[i * 3]] * 1e6
		M = M[lat < -20, :]
		im = ax[i, 1].contourf(lon, lat1, M, clevs, transform=ccrs.PlateCarree(),
				 vmin=clevs[0], vmax=clevs[-1], cmap=barra, alpha=0.5,
				 extend='both')
		lat1 = lat[lat < 0]
		u = var200[var_keys[i * 3 + 1]][lat1 < -30, :]
		v = var200[var_keys[i * 3 + 2]][lat1 < -30, :]
		lat1 = lat1[lat1 < -30]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.0007
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		im2 = ax[i, 1].quiver(lon[2:-1:10], lat1[2:-1:8], u_mask[2:-1:8, 2:-1:10],
			  v_mask[2:-1:8, 2:-1:10], transform=ccrs.PlateCarree(),
			  angles='xy', scale=25)
		ax[i, 1].coastlines(linewidth=0.4)
		#ax[i, 1].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax[i, 1].gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax[i, 1].xaxis.set_major_formatter(lon_formatter)
		ax[i, 1].yaxis.set_major_formatter(lat_formatter)
		if i==0:
			plt.quiverkey(im2, 0.9, 1.1, 0.5, r'$0.5 m^2/s^2$', labelpos='N',labelsep=0.03,
				      coordinates='axes', fontproperties={'size':7}, color='k')
	pad = 5
	for bx, col in zip(ax[0], lev):
		bx.annotate(col, xy=(0.5, 1), xytext=(0, pad),
			    xycoords='axes fraction', textcoords='offset points',
			    size='large', ha='center', va='baseline', fontsize=9)
	for bx, row in zip(ax[:, 0], tit):
		bx.annotate(row, xy=(0, 0.5), xytext=(-bx.yaxis.labelpad-pad, 0),
			    xycoords='axes fraction', rotation=90, textcoords='offset points',
			    size='large', ha='center', va='center', fontsize=6)
	fig.subplots_adjust(top=0.8, hspace=0.15)#, wspace=0.03)
	plt.suptitle(title, fontsize=10, x=0.49, y=0.92)
	plt.savefig(filename, dpi=500, bbox_inches='tight', papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

#def PlotCompositesdivPF(var50, var200, lat, lon, title, filename):
#	var_keys = [*var50]
#	tit = ['Weak-Strong SPV \n cond Ninio', 'Weak-Strong SPV \n cond Ninia',
#		   'Ninio-Ninia \n cond Weak SPV', 'Ninio-Ninia \n cond Strong SPV',
#		   'Ninio & \n Weak SPV', 'Ninia & \n Weak SPV',
#		   'Ninio & \n Strong SPV', 'Ninia & \n Strong SPV']
#	lev = ['50 hPa', '200 hPa']
#	proj = ccrs.PlateCarree(central_longitude=180)
#	fig, ax = plt.subplots(nrows=8, ncols=2, figsize=(5, 6.8), subplot_kw={'projection': proj})
#	clevs = np.arange(-0.75, 1.25, 0.5)
#	barra = plt.cm.RdBu_r
#	for i in range(8):
#		#ax = plt.subplot(8, 2, i * 2+ 1, projection=proj)
#		ax[i, 0].set_extent([0, 359, -90, -18], crs=ccrs.PlateCarree())
#		lat = lat[lat < 0]
#		lat1 = lat[lat < -20]
#		M = var50[var_keys[i * 3]] * 1e5
#		M = M[lat < -20, :]
#		ax[i, 0].contourf(lon, lat1, M, clevs, transform=ccrs.PlateCarree(),
#				 vmin=clevs[0], vmax=clevs[-1], cmap=barra, alpha=0.5,
#				 extend='both')
#		lat1 = lat[lat < 0]
#		u = var50[var_keys[i * 3 + 1]][lat1 < -30, :]
#		v = var50[var_keys[i * 3 + 2]][lat1 < -30, :]
#		lat1 = lat1[lat1 < -30]
#		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.0007
#		#mask array
#		u_mask = ma.array(u,mask = M)
#		v_mask = ma.array(v,mask = M)
#		ax[i, 0].quiver(lon[2:-1:10], lat1[2:-1:8], u_mask[2:-1:8, 2:-1:10],
#			  v_mask[2:-1:8, 2:-1:10], transform=ccrs.PlateCarree(),
#			  angles='xy', scale=25)
#		ax[i, 0].coastlines(linewidth=0.4)
#		#ax[i, 0].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
#		ax[i, 0].gridlines(crs=proj, linewidth=0.3, linestyle='-')
#		lon_formatter = LongitudeFormatter(zero_direction_label=True)
#		lat_formatter = LatitudeFormatter()
#		ax[i, 0].xaxis.set_major_formatter(lon_formatter)
#		ax[i, 0].yaxis.set_major_formatter(lat_formatter)
#
#		ax[i, 1].set_extent([0, 359, -90, -18], crs=ccrs.PlateCarree())
#		lat = lat[lat < 0]
#		lat1 = lat[lat < -20]
#		M = var200[var_keys[i * 3]] * 1e5
#		M = M[lat < -20, :]
#		im = ax[i, 1].contourf(lon, lat1, M, clevs, transform=ccrs.PlateCarree(),
#				 vmin=clevs[0], vmax=clevs[-1], cmap=barra, alpha=0.5,
#				 extend='both')
#		lat1 = lat[lat < 0]
#		u = var200[var_keys[i * 3 + 1]][lat1 < -30, :]
#		v = var200[var_keys[i * 3 + 2]][lat1 < -30, :]
#		lat1 = lat1[lat1 < -30]
#		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.0007
#		#mask array
#		u_mask = ma.array(u,mask = M)
#		v_mask = ma.array(v,mask = M)
#		im2 = ax[i, 1].quiver(lon[2:-1:10], lat1[2:-1:8], u_mask[2:-1:8, 2:-1:10],
#			  v_mask[2:-1:8, 2:-1:10], transform=ccrs.PlateCarree(),
#			  angles='xy', scale=25)
#		ax[i, 1].coastlines(linewidth=0.4)
#		#ax[i, 1].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
#		ax[i, 1].gridlines(crs=proj, linewidth=0.3, linestyle='-')
#		lon_formatter = LongitudeFormatter(zero_direction_label=True)
#		lat_formatter = LatitudeFormatter()
#		ax[i, 1].xaxis.set_major_formatter(lon_formatter)
#		ax[i, 1].yaxis.set_major_formatter(lat_formatter)
#		if i==0:
#			plt.quiverkey(im2, 0.9, 1.1, 0.5, r'$0.5 m^2/s^2$', labelpos='N',labelsep=0.03,
#				      coordinates='axes', fontproperties={'size':7}, color='k')
#
#	pad = 5
#	for bx, col in zip(ax[0], lev):
#		bx.annotate(col, xy=(0.5, 1), xytext=(0, pad),
#			    xycoords='axes fraction', textcoords='offset points',
#			    size='large', ha='center', va='baseline')
#	for bx, row in zip(ax[:, 0], tit):
#		bx.annotate(row, xy=(0, 0.5), xytext=(-bx.yaxis.labelpad-pad, 0),
#			    xycoords='axes fraction', rotation=90, textcoords='offset points',
#			    size='large', ha='center', va='center', fontsize=6)
#	fig.subplots_adjust(top=0.85)#, hspace=0.3, wspace=0.03)
#	plt.suptitle(title, fontsize=10, x=0.49, y=0.92)
#	plt.savefig(filename, dpi=500, bbox_inches='tight', papertype='A4')
#	plt.clf()
#	plt.cla()
#	plt.close()
#
