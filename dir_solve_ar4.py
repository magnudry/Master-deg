#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 17:24:35 2022

@author: magnudry
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.ndimage import gaussian_filter
import glob
from natsort import natsorted, ns

### GLOBAL CONSTANTS ###
#Drag coefficients
cdi = 0.0055 #ice-ocean
cda = 0.00125 #air-ocean
R = 3.5e-4 #ocean-bathymetry - need to find one that makes sense :)
#Densities
rho_a = 1.25 #air
rho_o = 1027.5 #ocean
#Frequencies
omega = (2 * np.pi) / 86400 #Earth's angular frequecy
#practicalities
th_fac = 5 #thinning factor (for slicing the grid)

### FUNCTIONS ###

#function for smoothing the depth array with a gaussian filter
def sm_depth_grid(grid, sig):
    sm_depth = gaussian_filter(grid["h"].data, sigma = sig) #arbitrary choice of sigma???
    da = xr.DataArray(sm_depth, dims = ["eta_rho","xi_rho"])
    grid["sm_h"] = da
    return grid

#function for adding an f/H dataarray to the dataset
def fonh_grid(grid):
    grid = grid.assign(fonh = grid["f"] / grid["sm_h"])
    #convert all zeros to nan, to avoid along-shore contour (takes a while)
    #grid = grid.where(grid["fonh"] > 0.0) fuck this
    return grid

#function for making contours for a given fonh value (level)
def make_contours(grid, level):
    cs = plt.contour(grid["xi_rho"],grid["eta_rho"],grid["fonh"], levels = [level])
    dat = cs.allsegs
    #exclude small contours (less than 50 points)
    newdat = []
    for i in range(len(dat[0])):
        if len(dat[0][i]) < 50:
            continue
        newdat.append(dat[0][i]) #put approved arrays into newdat
    return newdat

#function for making contours for multiple fonh values
def make_m_contours(grid, level_arr):
    cs = plt.contour(grid["xi_rho"],grid["eta_rho"],grid["fonh"], levels = level_arr)
    dat = cs.allsegs
    #exclude small contours (less than 50 points) and contours from South Atlantic 
    newdat = []
    for h in range(len(dat)):
        for i in range(len(dat[h])):
            if len(dat[h][i]) < 50 or dat[h][i][0][0] < 350:
                continue
            newdat.append(dat[h][i]) #put approved arrays into newdat
    return newdat, dat

#function for creating tau-arrays about a single contour
def cont_tau(contx, conty, file):
    txgrid = file["sustr"] 
    tygrid = file["svstr"] 
    conttx = txgrid.interp(xi_rho = contx, eta_rho = conty)
    contty = tygrid.interp(xi_rho = contx, eta_rho = conty)
    return conttx, contty

#function for creating depth averaged velocity-arrays about a single contour
def cont_vel(contx, conty, file):
    udarr = file["ubar"] 
    vdarr = file["vbar"] 
    contu = udarr.interp(xi_rho = contx, eta_rho = conty)
    contv = vdarr.interp(xi_rho = contx, eta_rho = conty)
    return contu, contv

#function for returning a dataArray of the depth values along contour
def cont_depth(contx, conty, grid):
    dgrid = grid["sm_h"] 
    contd = dgrid.interp(xi_rho = contx, eta_rho = conty)
    return contd
#interpoint distance function
def cont_dist(contx, conty, grid):
    dxgrid = 1 / grid["pm"] 
    dygrid = 1 / grid["pn"] 
    contdx = dxgrid.interp(xi_rho = contx, eta_rho = conty)
    contdy = dygrid.interp(xi_rho = contx, eta_rho = conty)
    return contdx, contdy
#fonh values on contour - hopefully quite similar throughout
def cont_fonh(contx, conty, grid):
    fhgrid = grid["fonh"]
    contfh = fhgrid.interp(xi_rho = contx, eta_rho = conty)
    return contfh
#function for preparing data for regridding using xesmf
def grid_prep(avg):
    #join dataArrays with similar dimensions
    x_set = xr.merge([avg["ubar"],avg["sustr"]])
    y_set = xr.merge([avg["vbar"],avg["svstr"]])
    #rename dims and coords to correspond with info-file names
    x_set = x_set.rename({"lat_u": "lat_r", "lon_u": "lon_r", "xi_u": "xi_rho","eta_u": "eta_rho"})
    y_set = y_set.rename({"lat_v": "lat_r", "lon_v": "lon_r", "xi_v": "xi_rho","eta_v": "eta_rho"})
    return x_set, y_set
"""
#function for extracting all necessary data for a contour to calculate circulation - returns two datasets
def cont_data(contour, grid, direc, t_steps, regrx, regry):
    #create xarray dataset for the contour
    #decompose contour in x and y coordinate dataArrays
    contx = xr.DataArray(contour[:,0], dims = "nz")
    conty = xr.DataArray(contour[:,1], dims = "nz")
    #create the constant dataArrays for the contour
    depth_set = cont_depth(contx,conty,grid)
    dx_set, dy_set = cont_dist(contx,conty,grid)
    #create dataset for the contour
    cont_data = xr.Dataset(data_vars = dict(
        contx = contx,
        conty = conty,
        depth = depth_set,
        distx = dx_set,
        disty = dy_set))
    #access the directory for U and V data
    files = glob.glob(direc + "/*.nc")[:t_steps]
    #o = cftime.datetime(2000,1,1,calendar="noleap") #temporal origin
    first = True
    for k in range(t_steps):
        #open NetCDFs 
        temp_file = xr.open_dataset(files[k])
        temp_file = temp_file.isel(ocean_time = 0)
        #prep files for regridding
        temp_x, temp_y = grid_prep(temp_file)
        #regrid the sets - must define the regridder objects in file where function is run
        x = regrx(temp_x)
        y = regry(temp_y)
        #merge regridded sets
        file = xr.merge([x,y])
        #create the dataArrays for single time steps
        print("lager tau")
        conttx, contty = cont_tau(contx,conty,file)
        print("lager vel")
        contu, contv = cont_vel(contx,conty,file)
        #include all dataArrays in temporary/current dataset
        #ufile.close()
        #vfile.close()
        print("lager dataset")
        curr_set = xr.Dataset(data_vars = dict(
            conttx = conttx,
            contty = contty,
            contu = contu,
            contv = contv
            ))       
        curr_set = curr_set.assign_coords(t_step = k) #zero indexed
        curr_set = curr_set.expand_dims("t_step")
        if first:
            vel_data = curr_set
            first = False
        else:
            #concatenate datasets along time dimension
            vel_data = xr.concat([vel_data,curr_set], dim = "t_step")
    #return the datasets with all necessary information to solve for circulation about the contour
    return cont_data, vel_data

def normalisers(cont_data):
    #calculate length of contour
    L_arr = (np.sqrt((cont_data["contx"].diff(dim = "nz")*cont_data["distx"][1:])**2 + 
                     (cont_data["conty"].diff(dim = "nz")*cont_data["disty"][1:])**2))
    L = L_arr.sum(dim = "nz")
    #calculate average depth along contour
    H_m = cont_data["depth"].mean(dim = "nz")
    return L, H_m

#function for solving for the average velocity about a contour, also returning the model velocity
def circ_solve(cont_data, vel, t_arr, delta_t, L, H_m):
    first = True
    for i in range(len(vel["t_step"])):
        t = t_arr[i]
        norm = np.exp(-1 * R * t / H_m.data) / L.data
        #norm = 1
        #print(norm)
        stress_da = (vel["conttx"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                              vel["contty"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        tau_da = (1 / (rho_o * cont_data["depth"][1:])) * np.exp(R * t / cont_data["depth"][1:]) * stress_da * delta_t
        #tau_da = (1 / (rho_o * cont_data["depth"][1:])) * np.exp(R * t / cont_data["depth"][1:]) * th_fac * (vel["conttx"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] + 
        #         vel["contty"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:]) * delta_t
        
        u_da = np.exp(R * t / cont_data["depth"][1:]) * (vel["contu"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                                                         vel["contv"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        #u_da = th_fac * (vel["meanu"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                         #vel["meanv"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        if first:
            stress_set = stress_da.sum(dim = "nz") / L.data
            tau_iter = tau_da.sum(dim = "nz") #must make two separate tau_arrays to be able to multiply in norm correctly
            tau_circ = tau_da.sum(dim = "nz") * norm
            u_circ = u_da.sum(dim = "nz") * norm
            first = False
        else:
            stress_set = xr.concat([stress_set,stress_da.sum(dim = "nz") / L.data], dim = "t_step")
            #tau_iter = xr.concat([tau_iter, tau_da.sum(dim = "nz")], dim = "t_step")
            tau_iter += tau_da.sum(dim = "nz")
            tau_circ = xr.concat([tau_circ, tau_iter * norm], dim = "t_step")
            #tau_circ = xr.concat([tau_circ, tau_iter.sum(dim = "t_step") * norm], dim = "t_step")
            u_circ = xr.concat([u_circ, u_da.sum(dim = "nz") * norm], dim = "t_step")
    #tau_circ = tau_circ.cumsum(dim = "t_step")
    return tau_circ, u_circ, stress_set

def file_merger(tau_circ, u_circ, stress, cont_data, vel_data, desc):
    #extract coords
    oc_t = stress.coords["ocean_time"]
    t_st = stress.coords["t_step"]
    #tau_circ does not have the same coordinates as the other two, create new tau with correct coords
    n_tau = xr.DataArray(data = tau_circ.data, dims = ["t_step"], coords = dict(ocean_time = (["t_step"], oc_t), t_step = t_st))
    n_tau = n_tau.to_dataset(name = "tau_circ")
    u_circ = u_circ.to_dataset(name = "u_circ")
    stress = stress.to_dataset(name = "stress")
    merged = xr.merge([n_tau,u_circ,stress,cont_data,vel_data])
    merged = merged.assign_attrs(dict(description = desc))
    return merged

def merged_to_disk(m_file, folder, name):
    m_file.to_netcdf(("{}{}{}".format(folder,name,".nc")))
    return
"""
#function for making datasets of the static contour info
def m_set(contours, grid, folder):
    for i in range(len(contours)):
        #decompose contour in x and y coordinate dataArrays
        contx = xr.DataArray(contours[i][:,0], dims = "nz")
        conty = xr.DataArray(contours[i][:,1], dims = "nz")
        #create the constant dataArrays for the contour
        depth_set = cont_depth(contx,conty,grid)
        #distance array
        dx_set, dy_set = cont_dist(contx,conty,grid)
        #fonh_array
        fh_set = cont_fonh(contx,conty,grid)
        #create dataset for the contour
        cont_data = xr.Dataset(data_vars = dict(
            contx = contx,
            conty = conty,
            depth = depth_set,
            fonh = fh_set,
            distx = dx_set,
            disty = dy_set))
        cont_data.to_netcdf(("{}{}{}{}".format(folder,"cont",str(i),".nc")))
    return
#function for obtaining the specific vorticity field in each time step
def vorticity(vgrid, grid):
    dvdx = vgrid["vbar"].diff(dim = "xi_rho") * grid["pn"].isel(xi_rho = slice(1,1602)) 
    dudy = vgrid["ubar"].diff(dim = "eta_rho") * grid["pm"].isel(eta_rho = slice(1,1202))
    vort = dvdx.isel(eta_rho = slice(1,1202)) - dudy.isel(xi_rho = slice(1,1602))
    #rvort = 
    return vort
#function for returning the specific vorticity for contour points
def cont_vort(contx,conty,vort):
    #vorticity field misses first row and column - reduce cont
    rcontx = contx - 1
    rconty = conty - 1
    contz = vort.interp(xi_rho = rcontx, eta_rho = rconty)
    return contz
#function for solving for the velocity and stress in time, for all contours
#goal: open each velocity dataset only once
#get static contour data from file directory
def conts_data(cdir, vdir, regrx, regry, folder):
    vfiles = sorted(glob.glob(vdir + "/*.nc")) #list of files w/ velocity data
    cfiles = glob.glob(cdir + "/*.nc") #list of files w/ contour data
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
        #interpolate to all contours
        for j in range(len(cfiles)):
            print("contour",j)
            #possibly necessary to open and close each static set in j-loop
            cont_d = xr.open_dataset(cfiles[j])
            #make tau
            conttx, contty = cont_tau(cont_d["contx"],cont_d["conty"],file)
            #make vel
            contu, contv = cont_vel(cont_d["contx"],cont_d["conty"],file)
            #make dataset
            curr_set = xr.Dataset(data_vars = dict(
                conttx = conttx,
                contty = contty,
                contu = contu,
                contv = contv
                )) 
            curr_set = curr_set.assign_coords({"t_step" : nr, "contour" : j}) #zero indexed
            curr_set = curr_set.expand_dims(["t_step","contour"])
            #write to (local) disk
            curr_set.to_netcdf(("{}{}{}{}{}{}{}{}".format(folder,"/cont",str(j),"/cont",str(j),"t_step",str(nr),".nc")))
    return
#rewrite the normalisers and circ_solve funcs to open data
def normalisers(cont_file):
    #open dataset
    cont_data = xr.open_dataset(cont_file)
    #calculate length of contour
    L_arr = (np.sqrt((cont_data["contx"].diff(dim = "nz")*cont_data["distx"][1:])**2 + 
                     (cont_data["conty"].diff(dim = "nz")*cont_data["disty"][1:])**2))
    L = L_arr.sum(dim = "nz")
    #calculate average depth along contour
    H_m = cont_data["depth"].mean(dim = "nz")
    return L, H_m

def circ_solve(cont_file, vel_dir, t_arr, delta_t, L, H_m, R):
    #open static dataset
    cont_data = xr.open_dataset(cont_file)
    #list of velocity files
    vel_arr = natsorted(glob.glob(vel_dir + "/*.nc"))
    #if arctic:
    vel_arr = vel_arr[700:]
    first = True
    for i in range(len(t_arr)):
        t = t_arr[i]
        vel = xr.open_dataset(vel_arr[i])
        vel = vel.isel(contour = 0, t_step = 0)
        norm = np.exp(-1 * R * t / H_m.data) / L.data
        stress_da = (vel["conttx"][1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                     vel["contty"][1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        tau_da = (1 / (rho_o * cont_data["depth"][1:])) * np.exp(R * t / cont_data["depth"][1:]) * stress_da * delta_t
        #tau_da = (1 / (rho_o * cont_data["depth"][1:])) * np.exp(R * t / cont_data["depth"][1:]) * th_fac * (vel["conttx"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] + 
        #         vel["contty"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:]) * delta_t
        
        u_da = np.exp(R * t / cont_data["depth"][1:]) * (vel["contu"][1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                                                         vel["contv"][1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        #u_da = th_fac * (vel["meanu"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                         #vel["meanv"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        if first:
            stress_set = stress_da.sum(dim = "nz") / L.data
            #include option to initialise tau-sum with u-data
            tau_iter = u_da.sum(dim = "nz")
            tau_circ = u_da.sum(dim = "nz") * norm
            #tau_iter = tau_da.sum(dim = "nz") #must make two separate tau_arrays to be able to multiply in norm correctly
            #tau_circ = tau_da.sum(dim = "nz") * norm
            u_circ = u_da.sum(dim = "nz") * norm
            first = False
        else:
            stress_set = xr.concat([stress_set,stress_da.sum(dim = "nz") / L.data], dim = "t_step")
            #tau_iter = xr.concat([tau_iter, tau_da.sum(dim = "nz")], dim = "t_step")
            tau_iter += tau_da.sum(dim = "nz")
            tau_circ = xr.concat([tau_circ, tau_iter * norm], dim = "t_step")
            #tau_circ = xr.concat([tau_circ, tau_iter.sum(dim = "t_step") * norm], dim = "t_step")
            u_circ = xr.concat([u_circ, u_da.sum(dim = "nz") * norm], dim = "t_step")
    #tau_circ = tau_circ.cumsum(dim = "t_step")
    return tau_circ, u_circ, stress_set

def file_merger(tau_circ, u_circ, stress, desc):
    #extract coords
    oc_t = stress.coords["ocean_time"]
    t_st = stress.coords["t_step"]
    #tau_circ does not have the same coordinates as the other two, create new tau with correct coords
    n_tau = xr.DataArray(data = tau_circ.data, dims = ["t_step"], coords = dict(ocean_time = (["t_step"], oc_t), t_step = t_st))
    n_tau = n_tau.to_dataset(name = "tau_circ")
    u_circ = u_circ.to_dataset(name = "u_circ")
    stress = stress.to_dataset(name = "stress")
    merged = xr.merge([n_tau,u_circ,stress])
    merged = merged.assign_attrs(dict(description = desc))
    return merged

def plotter(tau_c, u_c, s, days):
    #fig,axs = plt.subplots(2,1, sharex = True, figsize = (8,10))
    fig,axs = plt.subplots(2,1,sharex = True,gridspec_kw = {"height_ratios" : [3,1]},figsize=(8,6))
    #fig.suptitle("Circulation " + str(basin) + " fonh = " + str(fonh) + " R = "+ str(R))
    #fig.suptitle(inputname,fontsize = 25)
    axs[0].plot(days,-1*tau_c, label = "Linear model")
    #axs[0].set_ylim(-0.02,0.01)
    #axs[0].set_title("Mean contour velocity [m/s]", fontsize = 20)
    axs[0].grid(True)
    axs[0].set_ylabel(r"$u_{m}$ [m/s]", fontsize = 20)
    axs[0].plot(days,-1*u_c, label = "ROMS velocity data")
    axs[0].tick_params(labelsize = 15)
    #axs[0].legend()
    #axs[1].set_ylim(-0.01,0.02)
    axs[1].plot(days,-1*s)
    axs[1].tick_params(labelsize = 15)
    axs[1].set_ylabel(r"$\tau_{m}$ [Pa]", fontsize = 20)
    axs[1].grid(True)
    #axs[1].set_title("Contour average surface stress [Pa]", fontsize = 20)
    #if non arctic
    #axs[1].set_xlabel("days since 01.01.2012")
    #if arctic
    axs[1].set_xlabel("Days", fontsize = 20)
    fig.tight_layout(rect = (0,0,1,1))
    return fig

def plotter_v(tau_c, u_c, v_c, s, days, inputname): # basin, fonh, R):
    fig,axs = plt.subplots(2,1, sharex = True, gridspec_kw = {"height_ratios" : [3,1]}, figsize = (8,6))
    fig.suptitle(inputname,fontsize = "xx-large")
    #fig.suptitle("Circulation " + str(basin) + " fonh = " + str(fonh) + " R = "+ str(R))
    axs[0].plot(days,-1*tau_c, label = "Linear model")
    #axs[0].set_ylim(-0.2,0.2)
    axs[0].set_title("Mean contour velocity [m/s]",fontsize = "x-large")
    axs[0].grid(True)
    axs[0].plot(days,-1*u_c, label = "ROMS velocity data")
    axs[0].plot(days,-1*v_c, label = "Vorticity fluxes")
    #axs[0].plot(days,-1*tmv, label = "RHS total")
    axs[0].plot(days,-1*tau_c - v_c, label = "RHS total")
    #axs[0].legend()
    #axs[1].set_ylim(-0.01,0.02)
    axs[1].plot(days,-1*s)
    axs[1].grid(True)
    axs[1].set_title("Contour average surface stress [Pa]",fontsize = "x-large")
    axs[1].set_xlabel("days since 01.05.2016",fontsize = "large")
    fig.tight_layout(rect = (0,0,1,0.95))
    return fig
#function for a reduced plot of circulation
"""
def plotter_red(tau_c, u_c, s, days, inputname):
    fig,axs = plt.subplots(2,1,sharex = True,gridspec_kw = {"height_ratios" : [3,1]},figsize=(8,6))
    ax.set_title(inputname)
    twin = ax.twinx()
    pt, = twin.plot(days,-1*s,alpha = 0.4,color = "g")
    ax.plot(days,-1*tau_c, label = "Linear model")
    ax.plot(days,-1*u_c, label = "ROMS velocity data")
    ax.set_ylabel("Contour average velocity [m/s]")
    ax.set_xlabel("days since 01.12.2013")
    ax.grid(True) 
    ax.legend()
    #pt = twin.plot(days,-1*s, label = "Surface stress",alpha = 0.4)
    twin.grid(True)
    twin.yaxis.label.set_color(pt.get_color())
    twin.tick_params(axis='y', colors=pt.get_color())
    twin.set_ylabel("Surface stress [Pa]")
    fig.set_tight_layout(True)
    return fig
"""
#alternative variant of conts data that also includes the vorticity fluxes for a contour
def conts_data_v(cdir, vdir, grid, regrx, regry, folder):
    vfiles = sorted(glob.glob(vdir + "/*.nc")) #list of files w/ velocity data
    cfiles = glob.glob(cdir + "/*.nc") #list of files w/ contour data
    #grid_info = xr.open_dataset(grid)
    #iterate in time
    for i in range(len(vfiles)):
        nr = i 
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
        vort = vorticity(file, grid)
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
#same alternative variant for circ solve, that also solves for the vorticity flux
def circ_solve_v(cont_file, vel_dir, t_arr, delta_t, L, H_m, R):
    #open static dataset
    cont_data = xr.open_dataset(cont_file)
    #list of velocity files
    #vel_arr = glob.glob(vel_dir + "/*.nc")
    vel_arr = natsorted(glob.glob(vel_dir + "/*.nc"))
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
            #tau_iter = u_da.sum(dim = "nz")
            #tau_circ = u_da.sum(dim = "nz") * norm
            tau_iter = tau_da.sum(dim = "nz") #must make two separate tau_arrays to be able to multiply in norm correctly
            tau_circ = tau_da.sum(dim = "nz") * norm
            #could hard code the initial vorticity contribution to zero if tau_circ is initialised with u_0
            vort_iter = vort_da.sum(dim = "nz")
            vort_circ = vort_da.sum(dim = "nz") * norm
            u_circ = u_da.sum(dim = "nz") * norm
            tpv = u_da.sum(dim = "nz") * norm
            tmv = u_da.sum(dim = "nz") * norm
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
            tpv = xr.concat([tpv, (tau_iter + vort_iter) * norm], dim = "t_step")
            tmv = xr.concat([tmv, (tau_iter - vort_iter) * norm], dim = "t_step")
    #tau_circ = tau_circ.cumsum(dim = "t_step")
    return tau_circ, u_circ, vort_circ, tpv, tmv, stress_set
"""
def plotter_v(tau_c, u_c, v_c, s, days, inputname): # basin, fonh, R):
    fig,axs = plt.subplots(2,1, sharex = True, figsize = (8,10))
    fig.suptitle(inputname)
    #fig.suptitle("Circulation " + str(basin) + " fonh = " + str(fonh) + " R = "+ str(R))
    axs[0].plot(days,-1*tau_c, label = "Calculated tau contr")
    #axs[0].set_ylim(-0.2,0.2)
    axs[0].set_title("Mean contour velocity [m/s]")
    axs[0].grid(True)
    axs[0].plot(days,-1*u_c, label = "Model velocity data")
    axs[0].plot(days,v_c, label = "Vorticity fluxes")
    axs[0].plot(days,-1*tau_c + v_c, label = "RHS total")
    axs[0].legend()
    #axs[1].set_ylim(-0.01,0.02)
    axs[1].plot(days,-1*s)
    axs[1].grid(True)
    axs[1].set_title("Contour average surface stress [Pa]")
    axs[1].set_xlabel("days since 01.05.2016")
    return fig
"""
def file_merger_v(tau_circ, u_circ, v_circ, tpv, tmv, stress, desc):
    #extract coords
    oc_t = stress.coords["ocean_time"]
    t_st = stress.coords["t_step"]
    #tau_circ does not have the same coordinates as the other two, create new tau with correct coords
    n_tau = xr.DataArray(data = tau_circ.data, dims = ["t_step"], coords = dict(ocean_time = (["t_step"], oc_t), t_step = t_st))
    n_tau = n_tau.to_dataset(name = "tau_circ")
    n_vort = xr.DataArray(data = v_circ.data, dims = ["t_step"], coords = dict(ocean_time = (["t_step"], oc_t), t_step = t_st))
    n_vort = n_vort.to_dataset(name = "v_circ")
    n_tpv = xr.DataArray(data = tpv.data, dims = ["t_step"], coords = dict(ocean_time = (["t_step"], oc_t), t_step = t_st))
    n_tpv = n_tpv.to_dataset(name = "tpv")
    n_tmv = xr.DataArray(data = tmv.data, dims = ["t_step"], coords = dict(ocean_time = (["t_step"], oc_t), t_step = t_st))
    n_tmv = n_tmv.to_dataset(name = "tmv")
    u_circ = u_circ.to_dataset(name = "u_circ")
    stress = stress.to_dataset(name = "stress")
    merged = xr.merge([n_tau,u_circ,n_vort,stress,n_tpv,n_tmv])
    merged = merged.assign_attrs(dict(description = desc))
    return merged

#Function for calculating each term in equation Stokes (3.5)
def term_solver(cont_file,vel_dir,delta_t,R):
    #open static dataset
    cont_data = xr.open_dataset(cont_file)
    #list of velocity files
    #vel_arr = glob.glob(vel_dir + "/*.nc")
    vel_arr = natsorted(glob.glob(vel_dir + "/*.nc"))
    first = True
    for i in range(len(vel_arr)):
        vel = xr.open_dataset(vel_arr[i])
        vel = vel.isel(contour = 0, t_step = 0)
        if i == 0:
            unmo = vel["contu"]
            vnmo = vel["contv"]
            continue
        deru = (vel["contu"] - unmo) / delta_t
        derv = (vel["contv"] - vnmo) / delta_t 
        t1 = (deru[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] + 
              derv[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        t2 = vel["contz"][1:] * (vel["contu"][1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:] -
                                 vel["contv"][1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:])
        t3 = (1 / (rho_o * cont_data["depth"][1:])) * (vel["conttx"][1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                                                       vel["contty"][1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        t4 = (R / (cont_data["depth"][1:])) * (vel["contu"][1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                                               vel["contv"][1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        if first:
            der_arr = t1.sum(dim = "nz")
            vort_arr = t2.sum(dim = "nz")
            tau_arr = t3.sum(dim = "nz")
            R_arr = t4.sum(dim = "nz")
            first = False
        else:
            der_arr = xr.concat([der_arr, t1.sum(dim = "nz")], dim = "t_step")
            vort_arr = xr.concat([vort_arr, t2.sum(dim = "nz")], dim = "t_step")
            tau_arr = xr.concat([tau_arr, t3.sum(dim = "nz")], dim = "t_step")
            R_arr = xr.concat([R_arr, t4.sum(dim = "nz")], dim = "t_step")    
        unmo = vel["contu"]
        vnmo = vel["contv"]
    der_set = der_arr.to_dataset(name = "term1")
    vort_set = vort_arr.to_dataset(name = "term2")
    tau_set = tau_arr.to_dataset(name = "term3")
    R_set = R_arr.to_dataset(name = "term4")
    finished = xr.merge([der_set,vort_set,tau_set,R_set])
    return finished