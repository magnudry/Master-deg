#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 11:16:02 2022

@author: magnudry
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
#from scipy import interpolate
from scipy.ndimage import gaussian_filter
import glob
import xesmf as xe
from dir_solve_ar4 import *
#%%
#get the reduced files for regridding
info_file = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/regridder_files/deg_info.nc"
x_file = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/regridder_files/red_set_x.nc"
y_file = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/regridder_files/red_set_y.nc"

deg_info = xr.open_dataset(info_file)
red_set_x = xr.open_dataset(x_file)
red_set_y = xr.open_dataset(y_file) 
#%%
#create regridder objects
regrx = xe.Regridder(red_set_x, deg_info, "bilinear")
regry = xe.Regridder(red_set_y, deg_info, "bilinear")
#%%
#directory for the averaged 1-day files 
fdir = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4"
#file containing fonh-set
grid_meta = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/fonh_sig10.nc"
#%%
#open grid dataset
grid = xr.open_dataset(grid_meta)
#%%
filepath = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4/ocean_avg_day3409.nc"
vfile = xr.open_dataset(filepath)
#%%
print(vfile["lon_rho"][1,1])
#%%
xcomp,ycomp = grid_prep(vfile)
#%%
rx = regrx(xcomp)
ry = regry(ycomp)
both = xr.merge([rx,ry])
#%%
print(both)
print(grid)
#%%
vort = vorticity(both, grid)
#%%
print(vort)
#%%
#alternative variant of conts data that also includes the vorticity fluxes for a contour
def conts_data_v(cdir, vdir, grid, regrx, regry, folder):
    vfiles = sorted(glob.glob(vdir + "/*.nc")) #list of files w/ velocity data
    cfiles = glob.glob(cdir + "/*.nc") #list of files w/ contour data
    grid_info = xr.open_dataset(grid)
    #iterate in time
    for i in range(len(vfiles)):
        nr = i + 2799
        print("timestep",nr)
        #open NetCDFs 
        temp_file = xr.open_dataset(vfiles[i])
        temp_file = temp_file.isel(ocean_time = 0)
        #prep files for regridding
        temp_x, temp_y = grid_prep(temp_file)
        #regrid the sets - must define the regridder objects in file where function is run
        x = regrx(temp_x)
        y = regry(temp_y)
        #merge regridded sets
        file = xr.merge([x,y])
        #make vorticity da
        vort = vorticity(file, grid_info)
        #interpolate to all contours
        for j in range(len(cfiles)):
            print("contour",j)
            #possibly necessary to open and close each static set in j-loop
            cont_d = xr.open_dataset(cfiles[j])
            #make tau
            conttx, contty = cont_tau(cont_d["contx"],cont_d["conty"],file)
            #make vel
            contu, contv = cont_vel(cont_d["contx"],cont_d["conty"],file)
            #make vorticity
            contz = cont_vort(cont_d["contx"],cont_d["conty"],vort)
            #make dataset
            curr_set = xr.Dataset(data_vars = dict(
                conttx = conttx,
                contty = contty,
                contu = contu,
                contv = contv,
                contz = contz
                )) 
            curr_set = curr_set.assign_coords({"t_step" : nr, "contour" : j}) #zero indexed
            curr_set = curr_set.expand_dims(["t_step","contour"])
            #write to (local) disk
            curr_set.to_netcdf(("{}{}{}{}{}{}{}{}".format(folder,"/cont",str(j),"/cont",str(j),"t_step",str(nr),".nc")))
    return

def circ_solve_v(cont_file, vel_dir, t_arr, delta_t, L, H_m, R):
    #open static dataset
    cont_data = xr.open_dataset(cont_file)
    #list of velocity files
    vel_arr = glob.glob(vel_dir + "/*.nc")
    first = True
    for i in range(len(t_arr)):
        t = t_arr[i]
        vel = xr.open_dataset(vel_arr[i])
        vel = vel.isel(contour = 0, t_step = 0)
        norm = np.exp(-1 * R * t / H_m.data) / L.data
        stress_da = (vel["conttx"][1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                     vel["contty"][1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        tau_da = (1 / (rho_o * cont_data["depth"][1:])) * np.exp(R * t / cont_data["depth"][1:]) * stress_da * delta_t
        vort_da = np.exp(R * t / cont_data["depth"][1:]) * vel["contz"][1:] * (vel["contu"][1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:] -
                                                                               vel["contv"][1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:]) * delta_t        
        u_da = np.exp(R * t / cont_data["depth"][1:]) * (vel["contu"][1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                                                         vel["contv"][1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        if first:
            stress_set = stress_da.sum(dim = "nz") / L.data
            #include option to initialise tau-sum with u-data
            tau_iter = u_da.sum(dim = "nz")
            tau_circ = u_da.sum(dim = "nz") * norm
            #tau_iter = tau_da.sum(dim = "nz") #must make two separate tau_arrays to be able to multiply in norm correctly
            #tau_circ = tau_da.sum(dim = "nz") * 
            #could hard code the initial vorticity contribution to zero if tau_circ is initialised with u_0
            vort_iter = vort_da.sum(dim = "nz")
            vort_circ = vort_da.sum(dim = "nz")
            u_circ = u_da.sum(dim = "nz") * norm
            first = False
        else:
            stress_set = xr.concat([stress_set,stress_da.sum(dim = "nz") / L.data], dim = "t_step")
            #tau_iter = xr.concat([tau_iter, tau_da.sum(dim = "nz")], dim = "t_step")
            tau_iter += tau_da.sum(dim = "nz")
            tau_circ = xr.concat([tau_circ, tau_iter * norm], dim = "t_step")
            #tau_circ = xr.concat([tau_circ, tau_iter.sum(dim = "t_step") * norm], dim = "t_step")
            vort_iter += vort_da.sum(dim = "nz")
            vort_circ = xr.concat([vort_circ, vort_iter * norm], dim = "t_step")
            u_circ = xr.concat([u_circ, u_da.sum(dim = "nz") * norm], dim = "t_step")
    #tau_circ = tau_circ.cumsum(dim = "t_step")
    return tau_circ, u_circ, vort_circ, stress_set
#%%
#make the contour we want to examine
contdat = make_contours(grid, 5e-8)
#%%
print(contdat[1])
#%%
contd, vels = cont_data(contdat[1], grid, fdir, 608, regrx, regry)
#%%
print(vels)
#%%
L, H_m = normalisers(contd)
#%%
print(vels["conttx"])
#%%
times = np.arange(0,608,1) * 3600*24
days = np.arange(0,608,1)
#print(len(times))
#%%
tau_circ, u_circ, stress = circ_solve(contd, vels, times, 3600*24, L, H_m)
#%%
print(u_circ)
#%%
fig,axs = plt.subplots(2,1, sharex = True, figsize = (8,10))
fig.suptitle("Circulation NB, 4.8e-8, sigma10")
axs[0].plot(days,-1*tau_circ, label = "Calculated")
#axs[0].set_ylim(-0.02,0.01)
axs[0].set_title("Mean contour velocity [m/s]")
axs[0].grid(True)
axs[0].plot(days,-1*u_circ, label = "Model velocity data")
axs[0].legend()
#axs[1].set_ylim(-0.01,0.02)
axs[1].plot(days,-1*stress)
axs[1].grid(True)
axs[1].set_title("Contour average surface stress [Pa]")
axs[1].set_xlabel("days since 01.01.2012")
#%%
folder = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_figures/"
name = "NB4.8e-8sig10.pdf"
fig.savefig(("{}{}".format(folder,name)))

#%%
desc = "Norwegian Basin, f/H = 4.8e-8 (ms)-1, smoothing kernel 10 points"
merged = file_merger(tau_circ, u_circ, stress, contd, vels, desc)
#%%
print(merged)
#%%
folder = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Arctic4_output/"
name = "NB4.8e-8sig10"
merged_to_disk(merged,folder,name)