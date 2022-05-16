#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 10:39:47 2021

@author: magnudry
"""

import xarray as xr
import numpy as np

zfile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Fra_PEI/creg12.l75/CREG12.L75-REF08_mesh_zgr.nc"
zgrid = xr.open_dataset(zfile)
grdthn = zgrid.isel(t = 0, x = slice(0,1580,10), y = slice(0,1801,10))

ufile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/U_comp/CREG12.L75-REF08_y2000m01d05.5d_gridU.nc"
udata = xr.open_dataset(ufile)
udatathn = udata.isel(time_counter = 0, x = slice(0,1580,10), y = slice(0,1801,10))
#print(grdthn)
#print(udatathn)

udatathn = udatathn.assign(u_bot = udatathn["vozocrtx"].isel(depthu = grdthn["mbathy"]-1))

print(udatathn["u_bot"].data[105,105] == 0)