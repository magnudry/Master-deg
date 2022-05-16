#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 11:23:22 2021

@author: magnudry
"""
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import cftime
file = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/regridded_dataset_2003.nc"
base = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/IBCAO_v4_1_400m.nc"
mesh_hgr = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Fra_PEI/creg12.l75/CREG12.L75-REF08_mesh_hgr.nc"
mesh_zgr = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Fra_PEI/creg12.l75/CREG12.L75-REF08_mesh_zgr.nc"
mask = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Fra_PEI/creg12.l75/CREG12.L75-REF08_mask.nc"
U_0 = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/U_comp/CREG12.L75-REF08_y2000m01d05.5d_gridU.nc"

ufile = xr.open_dataset(U_0,use_cftime=True)
ufile = ufile.isel(time_counter = 0, x = slice(0,1580,5), y = slice(0,1801,5))
print(ufile["time_centered"].data)
#uvel = ufile["vozocrtx"].data
#print(uvel.shape)
#print(uvel)
"""
print(ufile["time_centered"].data[0])
print(dt.datetime(2000,1,1))


t = ufile["time_centered"].data[0]
o = cftime.datetime(2000,1,1,calendar="noleap")
print((t-o).total_seconds())
#delta = dt.timedelta(days=1)
#print(delta.total_seconds())
"""