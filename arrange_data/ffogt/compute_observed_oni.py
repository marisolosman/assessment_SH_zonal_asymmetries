#this code takes observed sst over ninio34 data and compute oni index to define el ninio years
import numpy as np
from scipy import signal
import xarray as xr
import pandas as pd
FILE = '/home/users/vg140344/datos/data/nino34.long.data'
data = open(FILE, 'r')
data.readline()
MISSING_VALUE = -99.99
sst = []
for line in data.readlines():
	s = str.split(line)
	if s[0]==str(MISSING_VALUE): break
	sst.append([float(s[i]) for i in range(len(s))])
data.close()

sst = np.array(sst)

years = sst[:, 0]

sst = sst[:, 1:]

#select the 1981-2018 period to compute oni

sst = sst[np.logical_and(years>=1981, years<=2018), :]

years = years[np.logical_and(years>=1981, years<=2018)]

#detrend data and compute 3month running mean

sst = np.reshape(sst, sst.shape[0] * sst.shape[1])
#compute mean to add after detrend it
sst_mean = np.mean(sst)
sst_det = signal.detrend(sst, type='linear') + sst_mean

#remove seasonal cycle
sst_det = np.reshape(sst_det, (np.int(np.size(sst_det)/12), 12))

sst_det = sst_det - np.mean(sst_det, 0)
sst_det = np.reshape(sst_det, sst_det.shape[0] * sst_det.shape[1])
MODE = 'same'

oni = np.convolve(sst_det, np.ones((3,))/3, mode=MODE)

#check if 5 consecutives 3month mean have values <> to +-0.5 to define ninio and ninia years
ninio = np.flatnonzero(np.convolve(oni>0.5, np.ones(5, dtype=int), 'valid')>=5)
ninia = np.flatnonzero(np.convolve(oni<-0.5, np.ones(5, dtype=int), 'valid')>=5)
ninio_dates = np.unique(np.array([np.arange(ninio[i], ninio[i] + 5) for i in range(np.size(ninio))]).ravel())
ninia_dates = np.unique(np.array([np.arange(ninia[i], ninia[i] + 5) for i in range(np.size(ninia))]).ravel())

enso = np.zeros_like(oni)

times = pd.date_range('1981-01-01', freq='M', periods=np.size(oni), tz='UTC')

enso[ninio_dates] = 1
enso[ninia_dates] = -1

enso = xr.DataArray(enso, coords=[times], dims=['time'])

enso.to_netcdf('/home/users/vg140344/datos/data/fogt/observed_enso.nc4')


