#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 13:25:31 2022

@author: magnudry
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
#from scipy import interpolate
from scipy.ndimage import gaussian_filter
import glob
import xesmf as xe
from dir_solve_ar4 import *
#%%
#get the reduced files for regridding
info_file = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/regridder_files/deg_info.nc"
x_file = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/regridder_files/red_set_x.nc"
y_file = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/regridder_files/red_set_y.nc"

deg_info = xr.open_dataset(info_file)
red_set_x = xr.open_dataset(x_file)
red_set_y = xr.open_dataset(y_file) 
#%%
#create regridder objects
regrx = xe.Regridder(red_set_x, deg_info, "bilinear")
regry = xe.Regridder(red_set_y, deg_info, "bilinear")
#%%
#directory for the averaged 1-day files 
fdir = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4"
#file containing fonh-set
grid_meta = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/fonh_sig10.nc"
#%%
#open grid dataset
grid = xr.open_dataset(grid_meta)
print(grid)
#%%
#contour level list - fonh-vals
levs = np.array([3.43,3.6,3.85,4,4.3,4.8,5.2,6]) * 1.0e-8
#nlevs = np.linspace(4.6e-8,2.9e-7,10)
#%%
#make contours
contsdat, dat = make_m_contours(grid, levs)
#%%
print(len(contsdat[0]))

#%%
folder = "/uio/hume/student-u71/magnudry/Documents/static_Ar4/"
#%%
#m_set(contsdat, grid, folder)
#%%
folderdyn = "/uio/hume/student-u71/magnudry/Documents/dyn_Ar4/"
cont_dir = "/uio/hume/student-u71/magnudry/Documents/static_Ar4/"
#%%
file_list = sorted(glob.glob(fdir + "/*.nc"))
cont_files = glob.glob(cont_dir + "/*.nc")
#%%
print(file_list[0])
#print(cont_files)
#%%
conts_data_v(cont_dir, fdir, grid,regrx,regry,folderdyn)
#%%
st_fold = "/uio/hume/student-u71/magnudry/Documents/static_Ar4/"
dy_fold = "/uio/hume/student-u71/magnudry/Documents/dynamic_Ar4/"
#%%
#start with contour 0, and see
st0 = st_fold + "cont0.nc"
dyn0 = dy_fold + "cont1/"
times = np.arange(0,1217) * 24*3600
#%%
check = xr.open_dataset(dyn0 + "cont1t_step1.nc")
print(check["ocean_time"])
#%%
L0, H_m0 = normalisers(st0)
t0,u0,s0 = circ_solve(st0,dyn0,times,24*3600,L0,H_m0)
#%%
print(s0)

#%%
desc = "skbda"
merg0 = file_merger(t0, u0, s0, desc)
#%%
print(merg0)
#%%
st_file = "/uio/hume/student-u71/magnudry/Documents/static_Ar4/cont29.nc"
ds = xr.open_dataset(st_file)
print(ds["contx"],ds["conty"],ds["fonh"])
#%%
day_arr = np.arange(0,1217) 
#%%
plotter(t0,u0,s0,day_arr, "EB", "3.43e-8")
#%%
from natsort import natsorted, ns
dy_fold = "/uio/hume/student-u71/magnudry/Documents/dynamic_Ar4/cont0/"
dy_list = natsorted(glob.glob(dy_fold + "/*.nc"))
#%%
print(dy_list[1218])
