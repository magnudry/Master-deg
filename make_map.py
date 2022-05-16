# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 11:51:59 2022

@author: dyrmo
"""
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cmocean
mpfile = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/fonh_sig10.nc"
A4means = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4means/A4means_regrid.nc"
mp = xr.open_dataset(mpfile)
means = xr.open_dataset(A4means)
#mp = mp.isel(x = slice(0,14550,10), y = slice(0,14550,10))
#%%
print(mp["lat_rho"].isel(eta_rho = 642)[1000:1100].max())
#%%
mpmin = mp.where(mp["lon_rho"] < 0, drop = True)
mpmax = mp.where(mp["lon_rho"] > 0, drop = True)
#%%
print(mpmin)
#print(mpmax)
#%%
plt.pcolormesh(mp["xi_rho"],mp["eta_rho"],mp["h"])
#%%
#The fonh levels array
levs = np.array([3.43,3.6,3.85,4,4.3,4.8,5.2,6]) * 1.0e-8
#%%
#make the map
fig = plt.figure(figsize=[15,15])
ax = plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo())
ax.set_extent([-180,180,60,90], crs=ccrs.PlateCarree())
ax.pcolormesh(mp["lon_rho"],mp["lat_rho"],mp["h"], transform = ccrs.PlateCarree(), cmap = cmocean.cm.deep)
ax.contour(mpmin["lon_rho"],mpmin["lat_rho"],mpmin["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.contour(mpmax["lon_rho"],mpmax["lat_rho"],mpmax["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax.coastlines()
fig.text(0.45,0.7,"Canadian Basin",backgroundcolor = "white",c = "k",size = "x-large",rotation = 45)
fig.text(0.57,0.62,"Makarov Basin",backgroundcolor = "white",c = "k",size = "x-large",rotation = 45)
fig.text(0.63,0.59,"Eurasian Basin",backgroundcolor = "white",c = "k",size = "x-large",rotation = 45)
fig.text(0.55,0.33,"Greenland Basin",backgroundcolor = "white",c = "k",size = "x-large",rotation = 45)
fig.text(0.58,0.27,"Norwegian Basin",backgroundcolor = "white",c = "k",size = "x-large",rotation = 45)

#%%
fig.savefig("/uio/lagringshotell/geofag/students/metos/magnudry/Master/overviewmap1.pdf")
#%%
plt.pcolormesh(mp["xi_rho"],mp["eta_rho"],mp["h"])
plt.contour(mp["xi_rho"],mp["eta_rho"],mp["fonh"],colors = "k", levels = levs)
#%%
#Norwegian Basin Map
fig = plt.figure(figsize=(8,8))
ax = plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo())
ax.set_extent([-10,12,62,73], crs=ccrs.PlateCarree())
ax.pcolormesh(mp["lon_rho"],mp["lat_rho"],mp["h"], transform = ccrs.PlateCarree(), cmap = cmocean.cm.deep)
ax.contour(mpmin["lon_rho"],mpmin["lat_rho"],mpmin["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.contour(mpmax["lon_rho"],mpmax["lat_rho"],mpmax["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax.coastlines()
fig.text(0.4,0.44, "2",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 12)
fig.text(0.4,0.36, "7",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 12)
fig.text(0.4,0.29, "14",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
fig.text(0.59,0.65, "15",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
fig.text(0.68,0.68,"20",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
fig.text(0.74,0.71, "25",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
#fig.set_tight_layout(True)
#%%
fig.savefig("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Maps/Norwegian_basin.png")
#%%
#Greenland Basin Map
fig = plt.figure(figsize= (8,8))
ax = plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo())
ax.set_extent([-16,12,71.5,79], crs=ccrs.PlateCarree())
ax.pcolormesh(mp["lon_rho"],mp["lat_rho"],mp["h"], transform = ccrs.PlateCarree(), cmap = cmocean.cm.deep)
ax.contour(mpmin["lon_rho"],mpmin["lat_rho"],mpmin["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.contour(mpmax["lon_rho"],mpmax["lat_rho"],mpmax["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax.coastlines()
fig.text(0.5,0.4,"8",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 12)
fig.text(0.5,0.35,"16",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
fig.text(0.5,0.28,"21",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
fig.text(0.55,0.25,"26",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
#%%
fig.savefig("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Maps/Greenland_basin.png")
#%%
#Joint Makarov and Eurasian Basin Map
fig = plt.figure(figsize=(8,10))
ax = plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo())
ax.set_extent([-15,160,80,90], crs=ccrs.PlateCarree())
ax.pcolormesh(mp["lon_rho"],mp["lat_rho"],mp["h"], transform = ccrs.PlateCarree(), cmap = cmocean.cm.deep)
ax.contour(mpmin["lon_rho"],mpmin["lat_rho"],mpmin["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.contour(mpmax["lon_rho"],mpmax["lat_rho"],mpmax["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax.coastlines()
#EB
fig.text(0.53,0.63,"0",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 11)
fig.text(0.6,0.69,"1",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 11)
fig.text(0.53,0.425,"3",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 11)
fig.text(0.655,0.735,"4",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 11)
fig.text(0.575,0.45,"9",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 11)
fig.text(0.68,0.59,"17",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 11)
fig.text(0.775,0.69,"22",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 11)
fig.text(0.792,0.73,"27",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 11)
#MB
fig.text(0.31,0.62,"5",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 11)
fig.text(0.32,0.58,"12",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 11)
fig.text(0.27,0.56,"18",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 11)
fig.text(0.23,0.54,"23",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 11)
#%%
fig.savefig("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Maps/Eurasian_basin.png")
#%%
#Canadian Basin Map
fig = plt.figure(figsize=(8,8))
ax = plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo())
ax.set_extent([-160,-100,71,90], crs=ccrs.PlateCarree())
ax.pcolormesh(mp["lon_rho"],mp["lat_rho"],mp["h"], transform = ccrs.PlateCarree(), cmap = cmocean.cm.deep)
ax.contour(mpmin["lon_rho"],mpmin["lat_rho"],mpmin["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.contour(mpmax["lon_rho"],mpmax["lat_rho"],mpmax["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax.coastlines()
fig.text(0.65,0.44,"6",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.2","color":"k"},c = "white",size = 12)
fig.text(0.64,0.36,"13",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
fig.text(0.64,0.32,"19",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
fig.text(0.62,0.25,"24",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
fig.text(0.73,0.35,"28",backgroundcolor = "k",bbox= {"boxstyle":"circle, pad=0.15","color":"k"},c = "white",size = 12)
#%%
fig.savefig("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Maps/Canadian_basin.png")
#%%
#Map with mean velocities
fig = plt.figure(figsize=[15,15])
ax = plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo())
ax.set_extent([-180,180,60,90], crs=ccrs.PlateCarree())
ax.pcolormesh(mp["lon_rho"],mp["lat_rho"],mp["h"], transform = ccrs.PlateCarree(), cmap = cmocean.cm.deep)
ax.contour(mpmin["lon_rho"],mpmin["lat_rho"],mpmin["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.contour(mpmax["lon_rho"],mpmax["lat_rho"],mpmax["fonh"],colors = "k", transform = ccrs.PlateCarree(), levels = levs, linewidths = 1.0)
ax.quiver(means["lon_r"].isel(xi_rho = slice(0,1602,20), eta_rho = slice(0,1202,20)),
          means["lat_r"].isel(xi_rho = slice(0,1602,20), eta_rho = slice(0,1202,20)),
          means["ubar"].isel(xi_rho = slice(0,1602,20), eta_rho = slice(0,1202,20)),
          means["vbar"].isel(xi_rho = slice(0,1602,20), eta_rho = slice(0,1202,20)), transform = ccrs.PlateCarree())
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax.coastlines()