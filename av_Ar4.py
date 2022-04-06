#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:45:48 2022

@author: magnudry
"""

import numpy as np
import xarray as xr
import xesmf as xe
import glob
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
#%%
#create list of files
files = glob.glob(fdir + "/*.nc")
#%%
#loop over files and add to new dataset
first = True
for i in range(len(files)):
    print(i)
    temp_file = xr.open_dataset(files[i])
    temp_file = temp_file.isel(ocean_time = 0)
    #prep files for regridding
    temp_x, temp_y = grid_prep(temp_file)
    #regrid the sets - must define the regridder objects in file where function is run
    x = regrx(temp_x)
    y = regry(temp_y)
    #merge regridded sets
    file = xr.merge([x,y])
    if first:
        usum = file["ubar"].data
        vsum = file["vbar"].data
        txsum = file["sustr"].data
        tysum = file["svstr"].data
        first = False
    else:
        usum += file["ubar"].data
        vsum += file["vbar"].data
        txsum += file["sustr"].data
        tysum += file["svstr"].data
#%%
#get means from the sum-arrays
meanu = usum / 608
meanv = vsum / 608
meantx = txsum / 608
meanty = tysum / 608
#%%
print(meanu.shape)
#%%
mean_vals = xr.Dataset(data_vars = dict(mubar = (["eta_rho","xi_rho"], meanu),
                                        mvbar = (["eta_rho","xi_rho"], meanv),
                                        msustr = (["eta_rho","xi_rho"], meantx),
                                        msvstr = (["eta_rho","xi_rho"], meanty)),
                       coords = dict(lon_r = (["eta_rho","xi_rho"], deg_info.coords["lon_r"].data),
                                     lat_r = (["eta_rho","xi_rho"], deg_info.coords["lat_r"].data)),
                       attrs = dict(description = "time averaged velocity and stress fields"))
#%%
print(mean_vals)
#%%
folder = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/"
name = "A4_averaged_fields.nc"
mean_vals.to_netcdf(("{}{}".format(folder,name)))