#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:59:50 2021

@author: magnudry
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

file = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Fra_PEI/creg12.l75/CREG12.L75-REF08_mesh_zgr.nc"
data = xr.open_dataset(file)
#print(data["mbathy"])


X = data["x"][::10].data
Y = data["y"][::10].data
#print(X)
topo = data["mbathy"].isel(t = 0, x = X, y = Y).data
#print(topo)

cs = plt.contour(X,Y,topo, levels = [60])
dat0 = cs.allsegs
newdat = []
#print(len(dat0))
for i in range(len(dat0[0])):
    if len(dat0[0][i]) < 50:
        continue
    newdat.append(dat0[0][i])
print(len(newdat))
"""
for el in newdat:
    plt.plot(el[:,0],el[:,1])
plt.show()


#plt.plot(dat0[:,0],dat0[:,1])
#plt.show()

"""

