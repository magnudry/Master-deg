#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 10:19:44 2021

@author: magnudry
"""

from solution_lin import *
import xarray as xr
import numpy as np
from scipy import interpolate
#see psi_solve for comments
zfile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Fra_PEI/creg12.l75/CREG12.L75-REF08_mesh_zgr.nc"
zgrid = xr.open_dataset(zfile)
grdthn = zgrid.isel(t = 0, x = slice(0,1580,10), y = slice(0,1801,10))
grdthn = f_cor(grdthn)
grdthn = depth_grid(grdthn)
grdthn = fonh_grid(grdthn)
xg = grdthn["x"].data
yg = grdthn["y"].data
ufile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/U_comp/CREG12.L75-REF08_y2000m01d05.5d_gridU.nc"
vfile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/V_comp/CREG12.L75-REF08_y2000m01d05.5d_gridV.nc"
udata = xr.open_dataset(ufile)
vdata = xr.open_dataset(vfile)
udatathn = udata.isel(time_counter = 0, x = slice(0,1580,10), y = slice(0,1801,10))
vdatathn = vdata.isel(time_counter = 0, x = slice(0,1580,10), y = slice(0,1801,10))
#print(vdatathn) #vomecrty
#print(udatathn["vozocrtx"].isel(x = 100, y = 150))

#fill all nans with 0.0 (even though bathymetry has value 0.0)
udatathn["vozocrtx"] = udatathn["vozocrtx"].fillna(0.0)
uvel = udatathn["vozocrtx"].data
vdatathn["vomecrty"] = vdatathn["vomecrty"].fillna(0.0)
vvel = vdatathn["vomecrty"].data
#create contours and pick one for Eurasian basin :)
contours = make_contours(grdthn, 2.355e-6)
c3 = contours[3]
#interpolation of each dataarray in uvel
contu = np.zeros((len(c3),len(uvel[0])))
contv = np.zeros((len(c3),len(vvel[0])))
for i in range(len(uvel)): #each depth layer
    fxu = interpolate.interp2d(xg,yg,uvel[i], kind = "linear")
    fyv = interpolate.interp2d(xg,yg,vvel[i], kind = "linear")
    for j in range(len(contu)): #each contour point
        contu[j,i] = fxu(c3[j][0], c3[j][1])
        contv[j,i] = fyv(c3[j][0], c3[j][1])

print(contu[0])
print(contv[0])

#functionality for obtaining transport along contour
#first, need to interpolate bathymetry 
fbath = interpolate.interp2d(xg, yg, grdthn["depth"].data, kind = "linear")
#depth-averaged velocities along contour - use these for u_0 in analysis
meanu = np.true_divide(contu.sum(1),(contu!=0).sum(1))
meanv = np.true_divide(contv.sum(1),(contv!=0).sum(1))
#empty transport arrays
tranu = np.zeros(len(meanu))
tranv = np.zeros(len(meanv))

for k in range(len(meanu)):
    bathymetry = fbath(c3[k][0],c3[k][1])
    tranu[k] = meanu[k] * bathymetry
    tranv[k] = meanv[k] * bathymetry


    