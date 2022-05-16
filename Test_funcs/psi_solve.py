#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 11:01:37 2021

@author: magnudry
"""

#take in contour and grid as parameters
#contour points need not correspond to grid points
#interpolate values of surrounding points - weighted by distance

from solution_lin import *
import xarray as xr
import numpy as np
from scipy import interpolate

zfile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Fra_PEI/creg12.l75/CREG12.L75-REF08_mesh_zgr.nc"
zgrid = xr.open_dataset(zfile)
grdthn = zgrid.isel(t = 0, x = slice(0,1580,10), y = slice(0,1801,10))
grdthn = f_cor(grdthn)
grdthn = fonh_grid(grdthn)
#print(grdthn)
xg = grdthn["x"].data
yg = grdthn["y"].data
#interpolation goes first
#need a file
ufile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/U_comp/CREG12.L75-REF08_y2000m01d05.5d_gridU.nc"
vfile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/V_comp/CREG12.L75-REF08_y2000m01d05.5d_gridV.nc"
udata = xr.open_dataset(ufile)
vdata = xr.open_dataset(vfile)
#print(udata["sozotaux"].attrs)

#interpolate surface stress first :)
udatathn = udata.isel(time_counter = 0, x = slice(0,1580,10), y = slice(0,1801,10))
udatathn["sozotaux"] = udatathn["sozotaux"].fillna(0.0)
taux = udatathn["sozotaux"].data
fx = interpolate.interp2d(xg,yg,taux, kind = "linear")

vdatathn = vdata.isel(time_counter = 0, x = slice(0,1580,10), y = slice(0,1801,10))
#print(vdatathn)
vdatathn["sometauy"] = vdatathn["sometauy"].fillna(0.0)
tauy = vdatathn["sometauy"].data
fy = interpolate.interp2d(xg, yg, tauy, kind = "linear")

contours = make_contours(grdthn, 2.355e-6)
c3 = contours[3]

conttaux = np.zeros(len(c3))
conttauy = np.zeros(len(c3))
for i in range(len(c3)):
    conttaux[i] = fx(c3[i][0],c3[i][1])
    conttauy[i] = fy(c3[i][0],c3[i][1])
print(conttaux)
print(conttauy)

