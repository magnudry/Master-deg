#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:53:04 2022

@author: magnudry
"""

from dir_solve_ar4 import *
from scipy.ndimage import gaussian_filter1d
folder = "/uio/hume/student-u71/magnudry/Documents/dyn_vort_Ar4/"
cfolder = "/uio/hume/student-u71/magnudry/Documents/static_Ar4/"
#conts = ["cont6","cont13","cont19","cont24","cont28"] #Canadian Basin
#cor_R = np.array([2.5,10,5,5,2.5]) * 1e-4 #Canadian Basin
conts = ["cont2","cont7","cont14","cont15","cont20","cont25"] #Norwegian Basin
cor_R = np.array([5,1.5,1.5,1.5,2.5,2.5]) * 1e-4 #Norwegian Basin
#%%
for i in range(6):
    print(i)
    vel_dir = folder + conts[i] + "/"
    cont_file = cfolder + conts[i] + ".nc"
    fin = term_solver(cont_file,vel_dir,24*3600,cor_R[i])
    fin.to_netcdf("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_vort_arr/" + conts[i] + ".nc")
#%%
#plot the terms
contour = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_vort_arr/cont25.nc"
ds = xr.open_dataset(contour)
days = np.arange(0,1043)
#%%
#Smooth all with a gaussian kernel of a month
term1 = gaussian_filter1d(ds["term1"], 5)
term2 = gaussian_filter1d(ds["term2"], 5)
term3 = gaussian_filter1d(ds["term3"], 5)
term4 = gaussian_filter1d(ds["term4"], 5)
res = term3-term1-term2-term4
#%%
fig,ax = plt.subplots(1,1,figsize=(8,6))
ax.plot(days,-1*term1,lw = 1,label = "Time deriv.")
ax.plot(days,-1*term2,lw = 1,label = "Vort. fluxes")
ax.plot(days,term3,lw = 1,label = "Surf. stress")
ax.plot(days,-1*term4,lw = 1,label = "Bott. stress") 
ax.plot(days,res, label = "Residual")
ax.set_xlabel("Days", fontsize = 20)
ax.set_ylabel("Circulation [m$^{2}$/s$^{2}$]", fontsize = 20)
ax.tick_params(labelsize = 15)
#ax.plot(days,ds["term3"]-ds["term1"]-ds["term2"]-ds["term4"])
ax.legend(fontsize=13)
ax.grid()
fig.tight_layout(rect = (0,0,1,1))

#%%
#loop over all contours in Norwegian and canadian basins
contours = natsorted(glob.glob("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_vort_arr/" + "/*.nc"))
print(contours[-1][75:-3])
#%%
days = np.arange(0,1043)
for i in range(len(contours)):
    contour = contours[i]
    ds = xr.open_dataset(contour)
    #Smooth all with a gaussian kernel of a month
    term1 = gaussian_filter1d(ds["term1"], 10)   
    erm2 = gaussian_filter1d(ds["term2"], 10)
    term3 = gaussian_filter1d(ds["term3"], 10)
    term4 = gaussian_filter1d(ds["term4"], 10)
    res = term3-term1-term2-term4
    fig,ax = plt.subplots(1,1,figsize=(8,6))
    ax.plot(days,-1*term1,ls = "-",lw = 1,label = "Time deriv.") 
    ax.plot(days,-1*term2,ls = "--",lw = 1,label = "Vort. fluxes")
    ax.plot(days,term3,ls = "-.",lw = 1,label = "Surf. stress")
    ax.plot(days,-1*term4,ls = ":",lw = 1.5,label = "Bott. stress") 
    ax.plot(days,-1*res, label = "Residual")
    ax.set_xlabel("Days", fontsize = 20)
    ax.set_ylabel("Circulation [m$^{2}$/s$^{2}$]", fontsize = 20)
    ax.tick_params(labelsize = 15)
    ax.legend(fontsize=13)
    ax.grid()
    fig.tight_layout(rect = (0,0,1,1))
    fig.savefig("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_vort_figs/" + contour[75:-3] + ".pdf")