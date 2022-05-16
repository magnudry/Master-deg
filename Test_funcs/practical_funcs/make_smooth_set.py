#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 11:13:17 2022

@author: magnudry
"""

from solution_lin import *
import xarray as xr
zgrid = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Fra_PEI/creg12.l75/CREG12.L75-REF08_mesh_zgr.nc"
folder = "/uio/hume/student-u71/magnudry/Documents/"
#%%
grid = xr.open_dataset(zgrid)
grid = grid.isel(t = 0)
#%%
#print(grid["mbathy"][1000,1000])
#%%
grid = f_cor(grid)
grid = depth_grid(grid)
#%%
print(grid.depth)
#%%
sig = 10
grid = sm_depth_grid(grid,sig)
#%%
print(grid)
#%%
depth_grid = xr.merge([grid["depth"].to_dataset(name="depth"),grid["sm_depth"].to_dataset(name="sm_depth")])
#%%
print(depth_grid)
#%%
depth_grid.to_netcdf(("{}{}".format(folder,"smooth_grid_karen.nc")))
#%%
grid = grid.isel(x = slice(0,1580,5), y = slice(0,1801,5))
#%%
#print(grid["mbathy"][0,0])
#%%
grid = fonh_grid(grid)
#%%
print(grid["fonh"][100,110])
#%%
grid.to_netcdf(("{}{}".format(folder,"fonh_thinned_sig15.nc")))
