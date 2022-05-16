#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 12:42:05 2022

@author: magnudry
"""

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from dir_solve_ar4 import *
from natsort import natsorted, ns
import glob
import cmocean
from scipy import signal
from scipy import stats
#from mpl_axes_aligner import align

files = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output/"
files_v = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output_v/"

imf = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_figures/"
imf_v = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_figures_v/"
#%%
bas = "CB"
basin = "Makarov Basin"
fonh = "3.85e-8"
#times for non arctic basins
#times = np.arange(0,2626) * 24*3600
#times for arctic basins, due to stress glitch - remove 700 first time steps
#times = np.arange(0,1926) * 24*3600
#times for vorticity incl
#times = np.arange(0,1044) * 24*3600
#%%
#open file
dyn_file = xr.open_dataset(files_v + basin + "/")
#%%
#lol make script
files_a = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output_vinit/"
baslist = natsorted(glob.glob(files_a + "/*.nc"))
#days = np.arange(0,1044)
corr_arr = np.zeros(240)
#coh_arr = np.zeros((240,2))
#%%
#print(coh_arr)
print(baslist[7][79:-3])
#%%
for x in range(len(baslist)):
    ds = xr.open_dataset(baslist[x])
    crcorr = stats.pearsonr((-1*ds["tmv"]),(-1*ds["u_circ"]))
    corr_arr[x] = crcorr[0]
    """
    crcoh = signal.coherence((-1*ds["tau_circ"]),(-1*ds["u_circ"]))
    cohims = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_coherence/"
    fig,ax = plt.subplots(1,1,figsize = (10,8))
    ax.plot(crcoh[0],crcoh[1], "o-")
    ax.set_title("Coherence " + str(baslist[x][82:-3]))
    ax.set_xscale("log")
    ax.set_xlabel("frequency [1/day]")
    
    fig.savefig(cohims + str(baslist[x][82:-3]) + ".pdf")
    """            
    #figure = plotter(ds["tau_circ"],ds["u_circ"],ds["stress"], days, str(baslist[x][82:-3]))
    #imf = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_figures_allconts/"
    #figure.savefig(("{}{}{}".format(imf,str(baslist[x][82:-3]),".pdf")))
#%%
ncarr = corr_arr.reshape(30,8)
#nbas = np.array(baslist).reshape(30,8)
#%%
print(ncarr[16])
#print(nbas[1])
#%%
#make plots of condensed data for a single contour - coherence, percentages or correlation of time series 
initfiles = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output_vinit/"
#%%
initlist = natsorted(glob.glob(initfiles + "/*.nc"))
#%%
print(initlist[202])
#%%
ds = xr.open_dataset("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output/Norwegian Basin/NB5.2e-8R0.00025.nc")
#ds = xr.open_dataset(initlist[202])
#pcdiff = (-1*ds["tau_circ"] + ds["v_circ"])/(-1*ds["u_circ"])
coh = signal.coherence((-1*ds["tau_circ"]),(-1*ds["u_circ"]))
print(len(coh[1]))
corr = stats.pearsonr((-1*ds["tau_circ"]),(-1*ds["u_circ"]))
psds = signal.periodogram(-1*ds["stress"])
psdu = signal.periodogram((-1*ds["u_circ"]))
psdt = signal.periodogram((-1*ds["tau_circ"]))
#%%
#print(corr[0])
plt.plot(coh[0],coh[1],"o-")
#plt.plot(psdu[0],psdu[1])
#plt.plot(psdt[0],psdt[1])
#plt.plot(psds[0],psds[1])
plt.xscale("log")
#%%
fig,ax = plt.subplots(1,1,figsize = (10,8))
for i in range(8):
    ds = xr.open_dataset(initlist[i])
    pcdiff = (-1*ds["tau_circ"] + ds["v_circ"])/(-1*ds["u_circ"])
    coh = signal.coherence((-1*ds["tau_circ"] + ds["v_circ"]),(-1*ds["u_circ"]))
    corr = stats.pearsonr((-1*ds["tau_circ"] + ds["v_circ"]),(-1*ds["u_circ"]))
#%%
#average the coherence data for each contour, with the best R
folder = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output_allconts/"
file_list_nb = ["cont2R0.0005.nc","cont7R0.00015000000000000001.nc","cont14R0.00015000000000000001.nc","cont15R0.00015000000000000001.nc","cont20R0.00025.nc","cont25R0.00025.nc"]
file_list_gb = ["cont8R0.0001.nc","cont16R0.00025.nc","cont21R0.0005.nc","cont26R0.00025.nc"]
file_list_eb = ["cont0R5e-05.nc","cont1R0.0002.nc","cont3R0.00025.nc","cont4R0.00025.nc","cont9R0.00025.nc","cont17R0.00025.nc","cont22R0.0005.nc","cont27R0.00025.nc"]
file_list_mb = ["cont5R0.0015.nc","cont12R5e-05.nc","cont18R5e-05.nc","cont23R5e-05.nc"]
file_list_cb = ["cont6R0.00025.nc","cont13R0.001.nc","cont19R0.0005.nc","cont24R0.0005.nc","cont28R0.00025.nc"]
#%%
file = file_list_cb[0]
days = np.arange(0,1926)
ds = xr.open_dataset(folder + file)
fig = plotter(ds["tau_circ"],ds["u_circ"],ds["stress"],days)
#fig.get_fontsize()
#%%
#make individual plot functions for each basin
def plotter_nb(list_nb):
    fig,ax = plt.subplots(6,2)
#%%
files = file_list_mb
days = np.arange(0,1926)
#print(str(files[2][:-3]))
imgf = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_red_figs/"
for y in range(len(files)):
    ds = xr.open_dataset(folder + files[y])
    figure = plotter(ds["tau_circ"],ds["u_circ"],ds["stress"],days)
    figure.savefig(imgf + str(files[y][:-3]) + ".pdf")
    #print(-1*np.mean(ds["stress"]))
    
#%%
files = file_list_cb
coh_array = np.zeros(129)
for y in range(len(files)):
    ds = xr.open_dataset(folder + files[y])
    coh = signal.coherence((-1*ds["tau_circ"]),(-1*ds["u_circ"]))
    coh_array += np.array(coh[1])
coh_array_m = coh_array / len(files)
#%%
print(coh_array_m)
#%%
plt.plot(coh[0],coh_array_m)
plt.xscale("log")
#%%
fig,ax = plt.subplots(1,1,figsize = (8,6))
ax.plot(coh[0],coh_array_m, "o-")
ax.set_ylabel("Averaged coherence magnitude", fontsize = 20)
ax.set_xscale("log")
ax.set_ylim(0,1)
ax.tick_params(labelsize = 15)
ax.set_xlabel("frequency [1/day]", fontsize = 20)
fig.set_tight_layout(True)
#%%
cohims = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_coherence/"
fig.savefig(cohims + "AC_Canadian_Basin_tl" + ".pdf")
#%%
#scatter plots of the time series
test_file = xr.open_dataset(folder + file_list_cb[2])
plt.scatter(-1*test_file["tau_circ"],-1*test_file["u_circ"],marker="+")
plt.grid()
#%%
#find the indexes where tau_circ is larger than zero
dslz = test_file.where(test_file["tau_circ"] < 0, drop = True)
dssz = test_file.where(test_file["tau_circ"] > 0, drop = True)
#%%
#find linear regression
resdslz = stats.linregress((-1*dslz["tau_circ"]),(-1*dslz["u_circ"]))
resdssz = stats.linregress((-1*dssz["tau_circ"]),(-1*dssz["u_circ"]))
#%%
#find correlation
cordslz = stats.pearsonr((-1*dslz["tau_circ"]),(-1*dslz["u_circ"]))
cordssz = stats.pearsonr((-1*dssz["tau_circ"]),(-1*dssz["u_circ"]))
cor = stats.pearsonr((-1*test_file["tau_circ"]),(-1*test_file["u_circ"]))
#%%
contourdat = file_list_cb[5]
file = xr.open_dataset(folder + contourdat)
#find the indexes where tau_circ is larger than zero
dslz = file.where(file["tau_circ"] < 0, drop = True)
dssz = file.where(file["tau_circ"] > 0, drop = True)
#find linear regression
resdslz = stats.linregress((-1*dslz["tau_circ"]),(-1*dslz["u_circ"]))
resdssz = stats.linregress((-1*dssz["tau_circ"]),(-1*dssz["u_circ"]))
fig,ax = plt.subplots(1,1,figsize=(8,6))
ax.scatter(-1*file["tau_circ"],-1*file["u_circ"],marker="+",alpha = 0.3,color = "lime")
#ax.set_title(str(contourdat[:-3]),fontsize = "xx-large")
ax.tick_params(labelsize = 15)
ax.set_xlabel("Linear model $u_{m}$ [m/s]",fontsize = 20)
ax.set_ylabel("ROMS model $u_{m}$ [m/s]",fontsize = 20)
ax.plot(-1*dslz["tau_circ"], resdslz.slope * -1*dslz["tau_circ"] + resdslz.intercept, label = "Positive values")
ax.plot(-1*dssz["tau_circ"], resdssz.slope * -1*dssz["tau_circ"] + resdssz.intercept, label = "Negative values")
ax.grid()
fig.tight_layout (rect = (0,0,1,1))
#ax.legend()
#%%
fig.savefig("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_regression_figs/" + str(contourdat[:-3]) + ".pdf")
#%%
print(resdslz.slope, resdslz.intercept)
print(resdssz.slope, resdssz.intercept)
#%%
#Automate for each basin
reg_arr_nb = np.zeros((len(file_list_nb),2))
reg_arr_gb = np.zeros((len(file_list_gb),2))
reg_arr_eb = np.zeros((len(file_list_eb),2))
reg_arr_mb = np.zeros((len(file_list_mb),2))
reg_arr_cb = np.zeros((len(file_list_cb),2))
#%%
regs = reg_arr_cb
files = file_list_cb
for i in range(len(files)):
    ds = xr.open_dataset(folder + files[i])
    #find the indexes where tau_circ is larger and smaller than zero
    dslz = ds.where(ds["tau_circ"] < 0, drop = True)
    dssz = ds.where(ds["tau_circ"] > 0, drop = True)
    #find linear regression
    resdslz = stats.linregress((-1*dslz["tau_circ"]),(-1*dslz["u_circ"]))   
    if len(dssz["u_circ"]) != 0:
        resdssz = stats.linregress((-1*dssz["tau_circ"]),(-1*dssz["u_circ"]))
        regs[i] = np.array([resdslz.slope,resdssz.slope])
    else:
        regs[i][0] = resdslz.slope
#%%
print(regs)
#%%
np.save("/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_regression/reg_arr_cb.npy",regs)