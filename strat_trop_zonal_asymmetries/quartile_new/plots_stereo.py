import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.path as mpath
from cartopy.util import add_cyclic_point

def PlotEnsoCompositesPoV(var1, var2, var3, var4, var5, var6, var7, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Ninio (All PoV)', 'Ninia (All PoV)',
		'Ninio (Weak PoV)', 'Ninia (Weak PoV)',
		'Ninio (Strong PoV)', 'Ninia (Strong PoV)']
	proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
	fig = plt.figure(1, (6, 10), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
		theta = np.linspace(0, 2 * np.pi, 100)
		center, radius = [0.5, 0.5], 0.5
		verts = np.vstack([np.sin(theta), np.cos(theta)]).T
		circle = mpath.Path(verts * radius + center)
		ax.set_boundary(circle, transform=ax.transAxes)
#		ax.contour(lon, lat, var7, np.arange(-200, 250, 50),
#			    transform=ccrs.PlateCarree(), colors='black')
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
				cmap=barra, extend='both', vmin=-60, vmax=60)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.85, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.22, 0.1, 0.5, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()
def PlotEnsoCompositesPoVDiff(var1, var2, var3, lat, lon, title, filename):
	vars = [var1, var2, var3]
	tit = ['Ninio - Ninia (All PoV)',
		'Ninio - Ninia (Weak SPoV)',
		'Ninio - Ninia (Strong SPoV']
	proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
	fig = plt.figure(1, (6, 10), 300)
	for i in range(3):
		ax = plt.subplot(3, 1, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
		theta = np.linspace(0, 2 * np.pi, 100)
		center, radius = [0.5, 0.5], 0.5
		verts = np.vstack([np.sin(theta), np.cos(theta)]).T
		circle = mpath.Path(verts * radius + center)
		ax.set_boundary(circle, transform=ax.transAxes)
#		ax.contour(lon, lat, var7, np.arange(-200, 250, 50),
#			    transform=ccrs.PlateCarree(), colors='black')
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
				cmap=barra, extend='both', vmin=-60, vmax=60)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.85, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.22, 0.1, 0.5, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotPoVCompositesENSO(var1, var2, var3, var4, var5, var6, var7, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Strong SPoV (All)', 'Weak SPoV (All)',
		'Strong SPoV (Ninio)', 'Weak SPoV (Ninio)',
		'Strong SPoV (Ninia)', 'Weak SPoV (Ninia)']
	proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
	fig = plt.figure(1, (6, 10), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
		theta = np.linspace(0, 2 * np.pi, 100)
		center, radius = [0.5, 0.5], 0.5
		verts = np.vstack([np.sin(theta), np.cos(theta)]).T
		circle = mpath.Path(verts * radius + center)
		ax.set_boundary(circle, transform=ax.transAxes)
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
				cmap=barra, extend='both', vmin=-60, vmax=60)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		#ax.contour(lon, lat, var7, np.arange(-200, 240, 40),
		#	    transform=ccrs.PlateCarree(), colors='black')
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.85, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.05])
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
	proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
	theta = np.linspace(0, 2 * np.pi, 100)
	center, radius = [0.5, 0.5], 0.5
	verts = np.vstack([np.sin(theta), np.cos(theta)]).T
	circle = mpath.Path(verts * radius + center)
	fig = plt.figure(1, (6, 10), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
		ax.set_boundary(circle, transform=ax.transAxes)
		im = ax.contourf(lon, lat, var[var_keys[i * 3]], clevs, transform=ccrs.PlateCarree(),
				cmap=barra, extend='both', vmin=-60, vmax=60)
		lat1 = lat[lat < 0]
		u = var[var_keys[i * 3 + 1]][lat1 < -18, :]
		v = var[var_keys[i * 3 + 2]][lat1 < -18, :]
		lat1 = lat1[lat1 < -18]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.07
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
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.85, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.25, 0.1, 0.5, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

