#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 13:10:50 2022

@author: magnudry
"""

#Try a numpy approach - xarray too slow, and data small enough to be loaded into memory
#insert all necessary data into a numpy ndarray and save as numpy file - wich may then be loaded into memory once per session

#plan: make u-velociy file, v-vel file, x-tau file, y-tau file and time-file

import numpy as np
import xarray as xr
import cftime
import glob
from solution_lin import *

udir = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/U_comp/"
vdir = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/V_comp/"
ufiles = glob.glob(udir + "/*.nc")
vfiles = glob.glob(vdir + "/*.nc") 
#%%
udata = np.zeros((((73,75,361,316)))) #fill with thinned datasets
vdata = np.zeros((((73,75,361,316))))
txdata = np.zeros(((73,361,316)))
tydata = np.zeros(((73,361,316)))
timedata = np.zeros(73) #fill with seconds from 1.1.2000
o = cftime.datetime(2000,1,1,calendar="noleap")
#%%
for i in range(73):
    #open datasets
    print(i+1)
    ufile = xr.open_dataset(ufiles[i],use_cftime=True)
    ufile = ufile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
    vfile = xr.open_dataset(vfiles[i],use_cftime=True)
    vfile = vfile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
    #fill arrays
    udata[i] = ufile["vozocrtx"].values
    vdata[i] = vfile["vomecrty"].values
    txdata[i] = ufile["sozotaux"].values
    tydata[i] = vfile["sometauy"].values
    timedata[i] = (ufile["time_centered"].data - o).total_seconds()
    ufile.close()
    vfile.close()
#%%
print(timedata)
#%%
np.save("/uio/lagringshotell/geofag/students/metos/magnudry/Master/timedata.npy",timedata)
#%%
np.save("/uio/lagringshotell/geofag/students/metos/magnudry/Master/txdata.npy",txdata)
#%%
np.save("/uio/lagringshotell/geofag/students/metos/magnudry/Master/tydata.npy",tydata)
#%%
np.save("/uio/lagringshotell/geofag/students/metos/magnudry/Master/udata.npy",udata)
#%%
np.save("/uio/lagringshotell/geofag/students/metos/magnudry/Master/vdata.npy",vdata)