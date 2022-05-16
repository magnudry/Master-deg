#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 11:25:58 2022

@author: magnudry
"""
import xarray as xr
import glob
from solution_lin import *
#%%
#test merge
ufile = xr.open_dataset("/uio/lagringshotell/geofag/students/metos/magnudry/Master/U_comp/CREG12.L75-REF08_y2000m01d05.5d_gridU.nc")
ufile = ufile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
ufile = ufile.reset_coords()
taux_arr = xr.merge([ufile["sozotaux"],ufile["utau_atmoce"],ufile["utau_iceoce"]])
print(taux_arr)
#%%
#transfer the necessary dataArrays to local disk for (hopefully) easier access
#udir = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/U_comp/"
vdir = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/V_comp/"
#ufiles = glob.glob(udir + "/*.nc") 
vfiles = glob.glob(vdir + "/*.nc")
folder = "/uio/hume/student-u71/magnudry/Documents/"
for i in range(73):
    print(i+1)
    """
    ufile = xr.open_dataset(ufiles[i])
    ufile = ufile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
    ufile = ufile.reset_coords()
    """
    vfile = xr.open_dataset(vfiles[i])
    vfile = vfile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
    vfile = vfile.reset_coords()
    #u_array = ufile["vozocrtx"]
    #v_array = vfile["vomecrty"]
    #taux_arr = xr.merge([ufile["sozotaux"],ufile["utau_atmoce"],ufile["utau_iceoce"]])
    #taux_arr = ufile["sozotaux"]
    tauy_arr = xr.merge([vfile["sometauy"],vfile["vtau_atmoce"],vfile["vtau_iceoce"]])
    #tauy_arr = vfile["sometauy"]
    #uname = ("{}{}{}".format("uvel",str(i),".nc"))
    #vname = ("{}{}{}".format("vvel",str(i),".nc"))
    #txname = ("{}{}{}".format("taux",str(i),".nc"))
    tyname = ("{}{}{}".format("tauy",str(i),".nc"))
    #u_array.to_netcdf(("{}{}{}".format(folder,"U-values/",uname)))
    #v_array.to_netcdf(("{}{}{}".format(folder,"V-values/",vname)))
    #taux_arr.to_netcdf(("{}{}{}".format(folder,"Taux/",txname)))
    tauy_arr.to_netcdf(("{}{}{}".format(folder,"Tauy/",tyname)))
    #ufile.close()
    vfile.close()
#%%
hgfile = "/uio/hume/student-u71/magnudry/Documents/CREG12.L75-REF08_mesh_hgr.nc"
hfile = xr.open_dataset(hgfile)
hfile = hfile.isel(t = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
folder = "/uio/hume/student-u71/magnudry/Documents/"
hfile.to_netcdf(("{}{}".format(folder,"mesh_hgrid.nc")))