#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 14:13:45 2022

@author: magnudry
"""

from dir_solve_ar4 import *
import glob
import xarray as xr
import numpy as np
#%%
#static folder
st_fold = "/uio/hume/student-u71/magnudry/Documents/static_Ar4/"
#dynamic folder
dy_fold = "/uio/hume/student-u71/magnudry/Documents/dynamic_Ar4/"
#dynamic folder with vorticity data
dy_v_fold = "/uio/hume/student-u71/magnudry/Documents/dyn_vort_Ar4/"
#full time series - no vorticity fluxes
#%%
#start with contour 0, and see
st = st_fold + "cont12.nc"
dyn = dy_fold + "cont12/"
bas = "MB"
basin = "Makarov Basin/"
fonh = "4e-8"
#times for non arctic basins
#times = np.arange(0,2626) * 24*3600
#times for arctic basins, due to stress glitch - remove 700 first time steps
times = np.arange(0,1926) * 24*3600
R = 2e-4 #ocean-bathymetry
name = bas + fonh + str(R)
#%%
L, H_m = normalisers(st)
t,u,s = circ_solve(st,dyn,times,24*3600,L,H_m,R)
#%%
#non arctic
#day_arr = np.arange(0,2626)
#arctic
day_arr = np.arange(0,1926)
#%%
figure = plotter(t,u,s,day_arr, bas, fonh, R)
#%%
imf = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_figures/" + basin
figure.savefig(("{}{}{}".format(imf,name,".pdf")))
#%%
#merge and save the data
outf = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output/" + basin
data_m = file_merger(t,u,s,name)
data_m.to_netcdf(outf + name + ".nc")
#%%
"""
#final three years - vorticity fluxes included
#%%
#start with contour 0, and see
st = st_fold + "cont20.nc"
dyn = dy_v_fold + "cont20/"
bas = "NB"
basin = "Norwegian Basin/"
fonh = "4.8e-8"
times = np.arange(0,1044) * 24*3600
R = 2e-4 #ocean-bathymetry
name = bas + fonh + str(R)
#%%
L, H_m = normalisers(st)
t,u,v,s = circ_solve_v(st,dyn,times,24*3600,L,H_m,R)
#%%
day_arr = np.arange(0,1044)
#%%
figure = plotter_v(t,u,v,s,day_arr, bas, fonh, R)
#%%
imf = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_figures_v/" + basin
#name = "NB5.2e-8R2.5.pdf"
figure.savefig(("{}{}{}".format(imf,name,".pdf")))
#%%
#merge and save
outf = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output_v/"
data_m = file_merger_v(t,u,v,s,name)
data_m.to_netcdf(outf + name + ".nc")
"""
#%%
#Automate!
contourlist = natsorted(glob.glob("/uio/hume/student-u71/magnudry/Documents/static_Ar4" + "/*.nc"))
dynfolder = "/uio/hume/student-u71/magnudry/Documents/dynamic_Ar4/"
R_vals = np.array([0.5,1,1.5,2,2.5,5,10,15]) * 1e-4
times = np.arange(0,1926) * 24*3600
#%%
#print(contourlist[9])
dyn = dynfolder + "cont" + str(3) 
exlist = natsorted(glob.glob(dyn + "/*.nc"))
print(exlist[1])
print(R_vals[4])
#%%
for x in range(len(contourlist)):
    print(x)
    st = contourlist[x]
    dyn = dynfolder + "cont" + str(x)
    L, H_m = normalisers(st)
    for y in range(len(R_vals)):
        R = R_vals[y]
        t,u,s = circ_solve(st,dyn,times,24*3600,L,H_m,R)
        #merge and save
        name = "cont" + str(x) + "R" + str(R)
        outf = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output_allconts/"
        data_m = file_merger(t,u,s,name)
        data_m.to_netcdf(outf + name + ".nc")