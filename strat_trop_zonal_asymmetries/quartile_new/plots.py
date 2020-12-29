import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.path as mpath

def PlotEnsoCompositesPoV(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Ninio (All PoV)', 'Ninia (All PoV)',
		'Ninio (Weak PoV)', 'Ninia (Weak PoV)',
		'Ninio (Strong PoV)', 'Ninia (Strong PoV)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
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

def PlotEnsoCompositesPoVex(var1, var2, var3, var4, lat, lon, title, filename):
	vars = [var1, var2, var3, var4]
	tit = ['Ninio (Neutral PoV)', 'Ninia (Neutral PoV)',
		'Ninio (Strong + Weak PoV)', 'Ninia (Strong + Weak PoV)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(4):
		ax = plt.subplot(2, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
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

def PlotKsEnsoCompositesPoV(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Ninio (All PoV) - Normal', 'Ninia (All PoV) - Normal',
		'Ninio (Weak PoV) - Normal', 'Ninia (Weak PoV) - Normal',
		'Ninio (Strong PoV) - Normal', 'Ninia (Strong PoV) - Normal']
	proj = ccrs.PlateCarree(central_longitude=180)
	clevs = [3] #np.arange(1, 5, 1)
	barra = plt.cm.OrRd_r
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
		im = ax.contour(lon, lat, vars[i], levels=clevs, transform=ccrs.PlateCarree())
#				cmap=barra, extend='both', vmin=clevs[0], vmax=clevs[-1])
		#barra.set_under(barra(0))
		#barra.set_over(barra(barra.N-1))
		ax.coastlines()
		#plt.clabel(im, inline=1, fontsize=10)
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	#fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	#cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	#fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotPoVCompositesENSO(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Strong SPoV (All)', 'Weak SPoV (All)',
		'Strong SPoV (Ninio)', 'Weak SPoV (Ninio)',
		'Strong SPoV (Ninia)', 'Weak SPoV (Ninia)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
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
def PlotPoVCompositesENSOex(var1, var2, var3, var4, lat, lon, title, filename):
	vars = [var1, var2, var3, var4]
	tit = ['Strong SPoV (Neutral ENSO)', 'Weak SPoV (Neutral)',
		'Strong SPoV (EN + LN)', 'Weak SPoV (EN + LN)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(4):
		ax = plt.subplot(2, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
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

def PlotEnsoCompositesPoVZPF(var, lat, lon, title, filename):
	var_keys = [*var]
	tit = ['Ninio (All PoV)', 'Ninia (All PoV)',
		'Ninio (Weak PoV)', 'Ninia (Weak PoV)',
		'Ninio (Strong PoV)', 'Ninia (Strong PoV)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		#im = ax.contour(lon, lat, var[var_keys[i * 3]], clevs, transform=ccrs.PlateCarree(),
		#		cmap=barra, extend='both', vmin=-60, vmax=60)
		lat1 = lat[lat < 0]
		u = var[var_keys[i * 3 + 1]][lat1 < -18, :]
		v = var[var_keys[i * 3 + 2]][lat1 < -18, :]
		lat1 = lat1[lat1 < -18]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.00000007
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6],
			  v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(),
			  angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3,
			  headaxislength=2.8)
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
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	#cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	#fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotEnsoCompositesPoVdivPF(var, lat, lon, title, filename):
	var_keys = [*var]
	tit = ['Ninio (All PoV)', 'Ninia (All PoV)',
		'Ninio (Weak PoV)', 'Ninia (Weak PoV)',
		'Ninio (Strong PoV)', 'Ninia (Strong PoV)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	clevs = np.arange(-0.75, 1.25, 0.5)
	barra = plt.cm.RdBu_r
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
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
		im2 = ax.quiver(lon[2:-1:10], lat1[2:-1:10], u_mask[2:-1:10, 2:-1:10],
			  v_mask[2:-1:10, 2:-1:10], transform=ccrs.PlateCarree(),
			  angles='xy', scale=20)
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
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()
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
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
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

def PlotPoVCompositesENSOU(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Strong SPoV (All)', 'Weak SPoV (All)',
		'Strong SPoV (Ninio)', 'Weak SPoV (Ninio)',
		'Strong SPoV (Ninia)', 'Weak SPoV (Ninia)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-7, 8, 1)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
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

def PlotPoVCompositesENSOZPF(var, lat, lon, title, filename):
	var_keys = [*var]
	print(var_keys)
	tit = ['Strong SPoV', 'Strong SPoV (El Ninio)', 'Strong SPoV (La Ninia)',
	       'Weak SPoV', 'Weak SPoV (El Ninio)', 'Weak SPoV (La Ninia)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	order = [1, 3, 5, 2, 4, 6]
	for i in range(6):
		ax = plt.subplot(3, 2, order[i], projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contour(lon, lat, var[var_keys[i * 3]], clevs, transform=ccrs.PlateCarree(),
				cmap=barra, extend='both', vmin=-60, vmax=60)
		lat1 = lat[lat < 0]
		u = var[var_keys[i * 3 + 1]][lat1 < -18, :]
		v = var[var_keys[i * 3 + 2]][lat1 < -18, :]
		lat1 = lat1[lat1 < -18]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.00000007
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6],
			  v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(),
			  angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3,
			  headaxislength=2.8)
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

def PlotPoVCompositesENSOdivPF(var, lat, lon, title, filename):
	var_keys = [*var]
	print(var_keys)
	tit = ['Strong SPoV', 'Strong SPoV (El Ninio)', 'Strong SPoV (La Ninia)',
	       'Weak SPoV', 'Weak SPoV (El Ninio)', 'Weak SPoV (La Ninia)']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	order = [1, 3, 5, 2, 4, 6]
	clevs = np.arange(-0.3, 0.6, 0.3)

	for i in range(6):
		ax = plt.subplot(3, 2, order[i], projection=proj)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, -20], crs=ccrs.PlateCarree())		
		im = ax.contourf(lon, lat, var[var_keys[i * 3]] * 1e7, clevs, linewidths=0.7,
				transform=ccrs.PlateCarree(), vmin=clevs[0], vmax=clevs[-1],
				cmap=barra, alpha=0.5, extend='both')

		lat1 = lat[lat < 0]
		u = var[var_keys[i * 3 + 1]][lat1 < -18, :]
		v = var[var_keys[i * 3 + 2]][lat1 < -18, :]
		lat1 = lat1[lat1 < -18]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.0007
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6],
			  v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(),
			  angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3,
			  headaxislength=2.8)
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
#	fig.subplots_adjust(right=0.8)
#	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
#	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
#	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotCompositesRWSChiWPVENSO(var_1, var_2, var_3, var_4, var_5, var_6,
			       u_1, u_2, u_3, u_4, u_5, u_6, v_1, v_2, v_3, v_4,
			       v_5, v_6, lat, lon, title, filename):
	var = [var_1, var_4, var_2, var_5, var_3, var_6]
	uchi = [u_1, u_4, u_2, u_5, u_3, u_6]
	vchi = [v_1, v_4, v_2, v_5, v_3, v_6]
	tit = ['Strong SPoV', 'Weak SPoV', 'Strong SPoV (El Ninio)', 'Weak SPoV (El Ninio)',
	       'Strong SPoV (La Ninia)', 'Weak SPoV (La Ninia)']
	clevs = np.arange(-1.2, 1.6, 0.8)
	barra = plt.cm.RdBu_r
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
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
		ax.quiver(lon[0:-1:6], lat[0:-1:6], u_mask[0:-1:6, 0:-1:6], v_mask[0:-1:6, 0:-1:6],
			  transform=ccrs.PlateCarree(), angles='xy')#, headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
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
def PlotCompositesRWSChiWENSOPV(var_1, var_2, var_3, var_4, var_5, var_6,
			       u_1, u_2, u_3, u_4, u_5, u_6, v_1, v_2, v_3, v_4,
			       v_5, v_6, lat, lon, title, filename):
	var = [var_1, var_4, var_2, var_5, var_3, var_6]
	uchi = [u_1, u_4, u_2, u_5, u_3, u_6]
	vchi = [v_1, v_4, v_2, v_5, v_3, v_6]
	tit = ['El Ninio', 'La Ninia', 'El Ninio (Weak SPoV)', 'La Ninia (Weak SPoV)',
	       'El Ninio (Strong SPoV)', 'La Ninia (Strong SPoV)']
	clevs = np.arange(-1.2, 1.6, 0.8)
	barra = plt.cm.RdBu_r
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
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
		ax.quiver(lon[0:-1:6], lat[0:-1:6], u_mask[0:-1:6, 0:-1:6], v_mask[0:-1:6, 0:-1:6],
			  transform=ccrs.PlateCarree(), angles='xy')#, headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
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

#========================================================================================================
def PlotEnsoCompositesPoVDiff(var1, var2, var3, lat, lon, title, filename):
	vars = [var1, var2, var3]
	tit = ['Ninio - Ninia (All PoV)',
		'Ninio - Ninia (In phase SPoV)',
		'Ninio - Ninia (Out of phase SPoV)']
	proj = ccrs.PlateCarree(central_longitude=180)
	barra = plt.cm.RdBu_r
	clevs = np.arange (-60, 70, 10)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(3):
		ax = plt.subplot(3, 1, i + 1, projection=proj)
		ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], levels=clevs, transform=ccrs.PlateCarree(),
				cmap=barra, extend='both', vmin=clevs[0], vmax=clevs[-1])
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

def PlotEnsoComposites(var_pos, var_neg, lat, lon, title, filename):
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

def PlotCompPlumbDiff(var, u, v, lat, lon, title, filename):
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 3.7), 300)
	ax = plt.subplot(projection=proj)
	clevs = np.arange(-60, 70, 10)
	barra = plt.cm.RdBu_r
	ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
	im = ax.contourf(lon, lat, var, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	lat1 = lat[lat < 0]
	u = u[lat1 < -18, :]
	v = v[lat1 < -18, :]
	lat1 = lat[lat < -18]
	Q50=np.percentile(np.sqrt(np.add(np.power(u,2),np.power(v,2))), 1) 
	M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.08
	#mask array
	u_mask = ma.array(u,mask = M)
	v_mask = ma.array(v,mask = M)
	ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.5, headlength=1.8, width=1.5e-3, headaxislength=2.5)
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

def PlotCompositesPlumb(var_pos, var_neg, u_pos, u_neg, v_pos, v_neg, lat,
			lon, title, filename):
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	ax = plt.subplot(2, 1, 1, projection=proj)
	clevs = np.arange(-60, 70, 10)
	barra = plt.cm.RdBu_r
	ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
	im = ax.contourf(lon, lat, var_pos, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	lat1 = lat[lat < 0]
	u_pos = u_pos[lat1 < -18, :]
	v_pos = v_pos[lat1 < -18, :]
	lat1 = lat[lat < -18]
	Q50=np.percentile(np.sqrt(np.add(np.power(u_pos,2),np.power(v_pos,2))), 1) 
	M = np.sqrt(np.add(np.power(u_pos,2),np.power(v_pos,2))) < 0.07
	#mask array
	u_mask = ma.array(u_pos,mask = M)
	v_mask = ma.array(v_pos,mask = M)
	ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
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
	ax1.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
	im = ax1.contourf(lon, lat, var_neg, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	lat1 = lat[lat < 0]
	u_neg = u_neg[lat1 < -18, :]
	v_neg = v_neg[lat1 < -18, :]
	lat1 = lat[lat < -18]
	Q50=np.percentile(np.sqrt(np.add(np.power(u_neg,2),np.power(v_neg,2))), 1) 
	M = np.sqrt(np.add(np.power(u_neg,2),np.power(v_neg,2))) < 0.07
	#mask array
	u_mask = ma.array(u_neg,mask = M)
	v_mask = ma.array(v_neg,mask = M)
	ax1.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
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

def PlotCompositesPlumbPV(var_pos, var_neg, u_pos, u_neg, v_pos, v_neg, lat,
			lon, title, filename):
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	ax = plt.subplot(2, 1, 1, projection=proj)
	clevs = np.arange(-60, 70, 10)
	barra = plt.cm.RdBu_r
	ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
	im = ax.contourf(lon, lat, var_pos, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	lat1 = lat[lat < 0]
	u_pos = u_pos[lat1 < -18, :]
	v_pos = v_pos[lat1 < -18, :]
	lat1 = lat[lat < -18]
	Q50=np.percentile(np.sqrt(np.add(np.power(u_pos,2),np.power(v_pos,2))), 1) 
	M = np.sqrt(np.add(np.power(u_pos,2),np.power(v_pos,2))) < 0.07
	#mask array
	u_mask = ma.array(u_pos,mask = M)
	v_mask = ma.array(v_pos,mask = M)
	ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
	ax.coastlines()
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	plt.title('Strong PV YEARS')

	ax1 = plt.subplot(2, 1, 2, projection=proj)
	barra = plt.cm.RdBu_r
	ax1.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
	im = ax1.contourf(lon, lat, var_neg, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	lat1 = lat[lat < 0]
	u_neg = u_neg[lat1 < -18, :]
	v_neg = v_neg[lat1 < -18, :]
	lat1 = lat[lat < -18]
	Q50=np.percentile(np.sqrt(np.add(np.power(u_neg,2),np.power(v_neg,2))), 1) 
	M = np.sqrt(np.add(np.power(u_neg,2),np.power(v_neg,2))) < 0.07
	#mask array
	u_mask = ma.array(u_neg,mask = M)
	v_mask = ma.array(v_neg,mask = M)
	ax1.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
	ax1.coastlines()
	ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
	ax1.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax1.xaxis.set_major_formatter(lon_formatter)
	ax1.yaxis.set_major_formatter(lat_formatter)
	plt.title('Weak PV YEARS')
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

def PlotCompDivPlumbDiff(var,u, v, lat, lon, title, filename):
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 3.7), 300)
	ax = plt.subplot(projection=proj)
	clevs = np.arange(-2.4, 4, 1.6)
	barra = plt.cm.RdBu_r
	ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
	im = ax.contourf(lon, lat, var / 10, clevs,
			transform=ccrs.PlateCarree(), cmap=barra, extend='both',
			vmin=-2.4, vmax=2.4)
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	M = np.sqrt(np.add(np.power(u, 2),np.power(v, 2))) < 0.5
	#mask array
	u_mask = ma.array(u, mask = M)
	v_mask = ma.array(v, mask = M)
	ax.quiver(lon[0:-1:8], lat[0:-1:8], u_mask[0:-1:8, 0:-1:8], v_mask[0:-1:8, 0:-1:8], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
	ax.coastlines()
	ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
	ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
	lon_formatter = LongitudeFormatter(zero_direction_label=True)
	lat_formatter = LatitudeFormatter()
	ax.xaxis.set_major_formatter(lon_formatter)
	ax.yaxis.set_major_formatter(lat_formatter)
	plt.colorbar(im, orientation='horizontal', shrink=0.25)
	plt.title(title)
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotEnsoCompositesSAM(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Ninio - Normal : All SAM', 'Ninia - Normal : All SAM',
		'Ninio - Normal : Weak SAM', 'Ninia - Normal : Weak SAM',
		'Ninio - Normal : Strong SAM', 'Ninia - Normal: Strong SAM']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
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

def PlotPoVCompositesDiffENSO(var1, var2, var3, lat, lon, title, filename):
	vars = [var1, var2, var3]
	tit = ['All', 'El Ninio', 'La Ninia']
	proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
	fig = plt.figure(1, (10, 6.7), 300)
	center, radius = [0.5, 0.5], 0.5
	theta = np.linspace(0, 2 * np.pi, 100)
	verts = np.vstack([np.sin(theta), np.cos(theta)]).T
	circle = mpath.Path(verts * radius + center)
	for i in range(3):
		ax = plt.subplot(3, 1, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, -20], crs=ccrs.PlateCarree())
		ax.set_boundary(circle, transform=ax.transAxes)
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
				cmap=barra, extend='both', vmin=-60, vmax=60)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		#ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		#ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		#lon_formatter = LongitudeFormatter(zero_direction_label=True)
		#lat_formatter = LatitudeFormatter()
		#ax.xaxis.set_major_formatter(lon_formatter)
		#ax.yaxis.set_major_formatter(lat_formatter)
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



