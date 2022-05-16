#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 11:27:31 2021

@author: magnudry
"""

#We expect that the response in the ocean current follows the frequency of the wind stress
#Thus, need to get frequency from tau

import numpy as np
#import scipy.integrate as integrate
#import scipy.special as special
from scipy.fft import fft, ifft
from scipy import interpolate
from scipy.ndimage import gaussian_filter
import xarray as xr
import matplotlib.pyplot as plt
#from solution_lin import *
import glob
import datetime as dt
import cftime
#%%
o = cftime.datetime(2000,1,1,calendar="noleap")
print(o)
#%%
udir = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/U_comp/"
vdir = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/V_comp/"
#make a list of the first n files in the directory
ufiles = glob.glob(udir + "/*.nc")[:20] 
vfiles = glob.glob(vdir + "/*.nc")[:20]
tau_arr = np.zeros(20)
time_arr = np.zeros(20)
#try only exracting a single point at first to create a tau-timeseries in that point
#say (900,900)
for i in range(20):
    ufile = xr.open_dataset(ufiles[i],use_cftime=True)
    vfile = xr.open_dataset(vfiles[i],use_cftime=True)
    taux = ufile["sozotaux"][0,900,900].data
    tauy = vfile["sometauy"][0,900,900].data
    tau_arr[i] = np.sqrt(taux**2 + tauy**2)
    time_arr[i] = (ufile["time_centered"][0].data - o).total_seconds() 
    ufile.close()
    vfile.close()
#%%
#print(tau_arr)
#print(time_arr)
f_tau = np.abs(fft(tau_arr))
#print(f_tau)
fr_arr = 1/time_arr
fig, ax = plt.subplots(1,1,figsize = [10,8])
ax.plot(fr_arr,f_tau)
#%%