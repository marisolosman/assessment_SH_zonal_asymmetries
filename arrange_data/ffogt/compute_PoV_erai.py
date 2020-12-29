#compute PoV index for ERAI data
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import eofdata
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '~/datos/data/'
FILE = 'hgt_erai_50.nc4'
hgt = xr.open_dataset(RUTA + FILE)
hgt = hgt.sel(**{'latitude': slice(-60, -90)}).mean(dim=['longitude', 'latitude'])
hgt.to_netcdf(RUTA + 'fogt/SPV_erai.nc4')

