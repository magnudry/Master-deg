#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 10:42:57 2021

@author: magnudry
"""

import numpy as np
import scipy.integrate as integrate
import scipy.special as special
from scipy.fft import fft, ifft
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr
import xarray as xr
import matplotlib.pyplot as plt
import glob
#import datetime as dt
import cftime


### GLOBAL CONSTANTS ###
#Drag coefficients
cdi = 0.0055 #ice-ocean
cda = 0.00125 #air-ocean
R = 5e-4 #ocean-bathymetry - need to find one that makes sense :)
#Densities
rho_a = 1.25 #air
rho_o = 1027.5 #ocean
#Frequencies
omega = (2 * np.pi) / 86400 #Earth's angular frequecy
#practicalities
th_fac = 5 #thinning factor (for slicing the grid)

### FUNCTIONS ###
#function for the surface stress
def tau(alpha, u_rel, u_a): #alpha = ice-concentration at cell
    t_i = rho_o * cdi * np.linalg.norm(u_rel) * u_rel
    t_a = rho_a * cda * np.linalg.norm(u_a) * u_a
    t = alpha * t_i + (1-alpha) * t_a
    return t

#function for the coriolis parameter
def f_cor(grid):
    grid = grid.assign(f = 2 * omega * np.sin(grid["nav_lat"] * np.pi/180))
    return grid

#function for creating a dataArray of depths to the dataset
def depth_grid(grid):
    grid = grid.assign(depth = grid["nav_lev"][grid["mbathy"]])
    #grid = grid.assign(smoothd = gaussian_filter(grid["depth"], sigma = 1))
    return grid
#function for smoothing the depth array with a gaussian filter
def sm_depth_grid(grid, sig):
    sm_depth = gaussian_filter(grid["depth"].data, sigma = sig) #arbitrary choice of sigma???
    da = xr.DataArray(sm_depth, dims = ["y","x"])
    grid["sm_depth"] = da
    return grid
#function for adding an f/H dataarray to the dataset
def fonh_grid(grid):
    grid = grid.assign(fonh = grid["f"] / grid["sm_depth"])
    #convert all zeros to nan, to avoid along-shore contour (takes a while)
    grid = grid.where(grid["fonh"] > 0.0)
    return grid

#assemble all addition functions in one - too heavy to run on local cpu apparently
def add_grids(grid):
    #coriolis parameter
    grid = grid.assign(f = 2 * omega * np.sin(grid["nav_lat"] * np.pi/180))
    #depth grid
    grid = grid.assign(depth = grid["nav_lev"][grid["mbathy"]])
    #smoothed depth grid
    sm_depth = gaussian_filter(grid["depth"].data, sigma = 5) #arbitrary choice of sigma???
    da = xr.DataArray(sm_depth, dims = ["y","x"])
    grid["sm_depth"] = da
    #fonh grid
    grid = grid.assign(fonh = grid["f"] / grid["sm_depth"])
    #convert all zeros to nan, to avoid along-shore contour (takes a while)
    grid = grid.where(grid["fonh"] != 0.)
    return grid

#function for making contours for a given fonh value (level)
def make_contours(grid, level):
    cs = plt.contour(grid["x"],grid["y"],grid["fonh"], levels = [level])
    dat = cs.allsegs
    #exclude small contours (less than 50 points)
    newdat = []
    for i in range(len(dat[0])):
        if len(dat[0][i]) < 50:
            continue
        newdat.append(dat[0][i]) #put approved arrays into newdat
    return newdat

#function for creating tau-arrays about a single contour
def cont_tau(contour, grid, ufile, vfile):
    #the files are already thinned    
    xg = grid["x"].data
    yg = grid["y"].data
    #fill nans w/ zero
    ufile["sozotaux"] = ufile["sozotaux"].fillna(0.0)
    vfile["sometauy"] = vfile["sometauy"].fillna(0.0)
    taux = ufile["sozotaux"].data
    fx = interpolate.interp2d(xg,yg,taux, kind = "linear")
    tauy = vfile["sometauy"].data
    fy = interpolate.interp2d(xg, yg, tauy, kind = "linear")
    conttaux = np.zeros(len(contour))
    conttauy = np.zeros(len(contour))
    for i in range(len(contour)):
        conttaux[i] = fx(contour[i][0],contour[i][1])
        conttauy[i] = fy(contour[i][0],contour[i][1])
    return conttaux, conttauy
"""
#function for creating velocity-arrays about a single contour - also depth averaged flow
def cont_vel(contour, grid, ufile, vfile):
    #the files are already thinned
    xg = grid["x"].data
    yg = grid["y"].data
    #fill all nans with 0.0 (even though bathymetry has value 0.0)
    ufile["vozocrtx"] = ufile["vozocrtx"].fillna(0.0)
    uvel = ufile["vozocrtx"].data
    vfile["vomecrty"] = vfile["vomecrty"].fillna(0.0)
    vvel = vfile["vomecrty"].data
    contu = np.zeros((len(contour),len(uvel[0])))
    contv = np.zeros((len(contour),len(vvel[0])))
    for i in range(len(uvel)): #each depth layer
        fxu = interpolate.interp2d(xg,yg,uvel[i], kind = "linear")
        fyv = interpolate.interp2d(xg,yg,vvel[i], kind = "linear")
        for j in range(len(contu)): #each contour point
            contu[j,i] = fxu(contour[j][0], contour[j][1])
            contv[j,i] = fyv(contour[j][0], contour[j][1])
    #find depth averaged velocities along contour
    meanu = np.true_divide(contu.sum(1),(contu!=0).sum(1))
    meanv = np.true_divide(contv.sum(1),(contv!=0).sum(1))
    #return contu, contv, meanu, meanv #return two 2D arrays and mean velocity arrays
    return meanu, meanv
#function for the depth at a contour point
def cont_depth(contour, grid):
    xg = grid["x"].data
    yg = grid["y"].data
    grid["sm_depth"] = grid["sm_depth"].fillna(0.0)
    fdepth = interpolate.interp2d(xg,yg,grid["sm_depth"].data, kind = "linear")
    depth_arr = np.zeros(len(contour))
    for i in range(len(contour)):
        depth_arr[i] = fdepth(contour[i,0],contour[i,1])
    return depth_arr    
#function for along-contour fonh-values
def cont_fonh(contour, grid):
    xg = grid["x"].data
    yg = grid["y"].data
    grid["fonh"] = grid["fonh"].fillna(0.0)
    ffonh = interpolate.interp2d(xg,yg,grid["fonh"].data, kind = "linear")
    fonh_arr = np.zeros(len(contour))
    for i in range(len(contour)):
        fonh_arr[i] = ffonh(contour[i,0],contour[i,1])
    return fonh_arr
#function for horizontal distances between points
def cont_dist(contour, hgrid):
    xg = hgrid["x"].data
    yg = hgrid["y"].data
    fe1 = interpolate.interp2d(xg,yg,hgrid["e1u"].data, kind = "linear")
    fe2 = interpolate.interp2d(xg,yg,hgrid["e2u"].data, kind = "linear")
    fe1_arr = np.zeros(len(contour))
    fe2_arr = np.zeros(len(contour))
    for i in range(len(contour)):
        fe1_arr[i] = fe1(contour[i,0],contour[i,1])
        fe2_arr[i] = fe2(contour[i,0],contour[i,1])
    return fe1_arr, fe2_arr
"""
### Variant functions of the contour functions above, they use xarray as much as possible ###
#function for creating tau-arrays about a single contour
def cont_tau_var(contx, conty, xfile, yfile):
    txgrid = xfile["sozotaux"] 
    tygrid = yfile["sometauy"] 
    conttx = txgrid.interp(x = contx, y = conty)
    contty = tygrid.interp(x = contx, y = conty)
    return conttx, contty
#function for creating velocity-arrays about a single contour - also depth averaged flow
#try to make a cont_vel with xarray functions without loading data unnecessarily - returns only xarray.dataArrays
def cont_vel_var(contx, conty, ufile, vfile):
    udarr = ufile["vozocrtx"] 
    vdarr = vfile["vomecrty"] 
    contu = udarr.interp(x = contx, y = conty)
    contv = vdarr.interp(x = contx, y = conty)
    #meanu, meanv = cont_mean_vel(contu, contv)
    meanu, meanv = cont_mean_upper(contu, contv)
    return contu, contv, meanu, meanv

def cont_mean_vel(contu,contv):
    #make weighted u and v arrays to accurately calulate depth averaged velocities
    uweight = xr.DataArray(np.diff(contu.coords["depthu"],prepend=0),dims = "depthu")
    vweight = xr.DataArray(np.diff(contv.coords["depthv"],prepend=0),dims = "depthv")
    #calculate the weighted sums in u and v
    wu = contu.weighted(uweight).sum(dim = "depthu")
    wv = contv.weighted(vweight).sum(dim = "depthv")
    #now locate the depth - count nonzero returns array of contour length with the indices corresponding to final nonzero value
    du = np.count_nonzero(contu, axis = 0)
    dv = np.count_nonzero(contv, axis = 0)
    #make depth arrays
    duarr = xr.DataArray(contu.coords["depthu"].isel(depthu = du).data,dims = "nz")
    dvarr = xr.DataArray(contv.coords["depthv"].isel(depthv = dv).data,dims = "nz")
    meanu = wu / duarr
    meanv = wv / dvarr
    return meanu, meanv
#function for calculating mean velocity of the upper 200 meters of column
def cont_mean_upper(contu,contv):
    uweight = xr.DataArray(np.diff(contu.coords["depthu"][:31],prepend=0),dims = "depthu")
    vweight = xr.DataArray(np.diff(contv.coords["depthv"][:31],prepend=0),dims = "depthv")
    wu = contu.isel(depthu = slice(0,31)).weighted(uweight).sum(dim = "depthu")
    wv = contv.isel(depthv = slice(0,31)).weighted(vweight).sum(dim = "depthv")
    meanu = wu / contu.coords["depthu"][30].data
    meanv = wv / contu.coords["depthu"][30].data
    return meanu, meanv

#function for returning a dataArray of the depth values along contour
def cont_depth_var(contx, conty, grid):
    dgrid = grid["sm_depth"] 
    contd = dgrid.interp(x = contx, y = conty)
    return contd

def cont_dist_var(contx, conty, hgrid):
    dxgrid = hgrid["e1u"] 
    dygrid = hgrid["e2u"] 
    contdx = dxgrid.interp(x = contx, y = conty)
    contdy = dygrid.interp(x = contx, y = conty)
    return contdx, contdy
#function for extracting all necessary data for a contour to calculate circulation - returns two datasets
def cont_data_var(contour, zgrid, hgrid, udir, vdir, tauxdir, tauydir, t_steps):
    #create xarray dataset for the contour
    #decompose contour in x and y coordinate dataArrays
    contx = xr.DataArray(contour[:,0], dims = "nz")
    conty = xr.DataArray(contour[:,1], dims = "nz")
    #create the constant dataArrays for the contour
    depth_set = cont_depth_var(contx,conty,zgrid)
    dx_set, dy_set = cont_dist_var(contx,conty,hgrid)
    #create dataset for the contour
    cont_data = xr.Dataset(data_vars = dict(
        contx = contx,
        conty = conty,
        depth = depth_set,
        distx = dx_set,
        disty = dy_set))
    #access the directiories for U and V data
    ufiles = glob.glob(udir + "/*.nc")[:t_steps] 
    vfiles = glob.glob(vdir + "/*.nc")[:t_steps] 
    txfiles = glob.glob(tauxdir + "/*.nc")[:t_steps]
    tyfiles = glob.glob(tauydir + "/*.nc")[:t_steps]
    #o = cftime.datetime(2000,1,1,calendar="noleap") #temporal origin
    first = True
    for k in range(t_steps):
        #open NetCDFs 
        txfile = xr.open_dataset(txfiles[k])
        tyfile = xr.open_dataset(tyfiles[k])
        ufile = xr.open_dataset(ufiles[k])
        vfile = xr.open_dataset(vfiles[k])
        #ufile.load()
        #vfile.load()
        #create the dataArrays for single time steps
        print("lager tau")
        conttx, contty = cont_tau_var(contx,conty,txfile,tyfile)
        print("lager vel")
        contu, contv, meanu, meanv = cont_vel_var(contx,conty,ufile,vfile)
        #include all dataArrays in temporary/current dataset
        #ufile.close()
        #vfile.close()
        print("lager dataset")
        curr_set = xr.Dataset(data_vars = dict(
            conttx = conttx,
            contty = contty,
            contu = contu,
            contv = contv,
            meanu = meanu,
            meanv = meanv
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
#function for interpolating tau_data in time to increase temporal resolution
def tau_interp(vel_data,times):
    ntaux = vel_data["conttx"].interp(t_step = times)
    ntauy = vel_data["contty"].interp(t_step = times)
    return ntaux, ntauy
def normalisers(cont_data):
    #calculate length of contour
    L_arr = (np.sqrt((cont_data["contx"].diff(dim = "nz")*cont_data["distx"][1:])**2 + 
                     (cont_data["conty"].diff(dim = "nz")*cont_data["disty"][1:])**2) * th_fac)
    L = L_arr.sum(dim = "nz")
    #calculate average depth along contour
    H_m = cont_data["depth"].mean(dim = "nz")
    return L, H_m
#function solving the circulation about contour for both the LHS and RHS of circulation equation
def circ_solve_tau(cont_data, ntaux, ntauy, t_arr, delta_t, L, H_m):
    first = True
    for i in range(len(ntaux["t_step"])):
        t = t_arr[i]
        norm = np.exp(-1 * R * t / H_m.data) / L.data
        stress_da = th_fac * (ntaux.isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                              ntauy.isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        tau_da = (1 / (rho_o * cont_data["depth"][1:])) * np.exp(R * t / cont_data["depth"][1:]) * stress_da * delta_t
        if first:
            stress_set = stress_da.sum(dim = "nz") / L.data
            tau_iter = tau_da.sum(dim = "nz") #must make two separate tau_arrays to be able to multiply in norm correctly
            tau_circ = tau_da.sum(dim = "nz") * norm
            first = False
        else:
            stress_set = xr.concat([stress_set,stress_da.sum(dim = "nz") / L.data], dim = "t_step")
            tau_iter += tau_da.sum(dim = "nz") #xr.concat([tau_iter, tau_da.sum(dim = "nz")], dim = "t_step")
            tau_circ = xr.concat([tau_circ, tau_iter * norm], dim = "t_step")
    return tau_circ, stress_set
def circ_solve_u(cont_data, vel, t_arr, delta_t, L, H_m):
    first = True
    for i in range(len(vel["t_step"])):
        t = t_arr[i]
        norm = np.exp(-1 * R * t / H_m.data) / L.data
        u_da = np.exp(R * t / cont_data["depth"][1:]) * th_fac * (vel["meanu"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                                                         vel["meanv"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        if first:
            u_circ = u_da.sum(dim = "nz") * norm
            first = False
        else:
            u_circ = xr.concat([u_circ, u_da.sum(dim = "nz") * norm], dim = "t_step")
    return u_circ
def circ_solve_var(cont_data, vel, t_arr, delta_t, L, H_m):
    first = True
    for i in range(len(vel["t_step"])):
        t = t_arr[i]
        norm = np.exp(-1 * R * t / H_m.data) / L.data
        #norm = 1
        #print(norm)
        stress_da = th_fac * (vel["conttx"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                              vel["contty"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
        tau_da = (1 / (rho_o * cont_data["depth"][1:])) * np.exp(R * t / cont_data["depth"][1:]) * stress_da * delta_t
        #tau_da = (1 / (rho_o * cont_data["depth"][1:])) * np.exp(R * t / cont_data["depth"][1:]) * th_fac * (vel["conttx"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] + 
        #         vel["contty"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:]) * delta_t
        
        u_da = np.exp(R * t / cont_data["depth"][1:]) * th_fac * (vel["meanu"].isel(t_step = i)[1:] * cont_data["contx"].diff(dim = "nz") * cont_data["distx"][1:] +
                                                         vel["meanv"].isel(t_step = i)[1:] * cont_data["conty"].diff(dim = "nz") * cont_data["disty"][1:])
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
#function for saving a contour dataset to disk
def circ_to_disk(tau_circ, u_circ, stress, folder, basin, fonh_val, sig):
    comb_set = xr.merge([tau_circ, u_circ, stress])
    fname = ("{}{}{}{}{}{}".format(folder,basin,str(fonh_val),"sig",str(sig),".nc"))
    comb_set.to_netcdf(fname)
    return

"""         
    tau_circ_arr = np.zeros(t_steps)
    u_circ_arr = np.zeros(t_steps)
    tau_circ = 0
    #u_circ = 0
    for m in range(t_steps):
        inc_tau = 0 #tau circulation increment
        u_circ = 0 #u circulation
        for l in range(len(contour)):
            #can calculate u_circulation here as well
            tscaler = np.exp(R*time_arr[m] / darr[l]) / (rho_o * darr[l])
            uscaler = np.exp(R*time_arr[m] / darr[l])
            if l == len(contour) -1:
                nxt = contour[0]
            else:
                nxt = contour[l+1]
            dl = np.subtract(nxt,contour[l])
            #tau x- and y-components
            x_comp = tscaler * taux_arr[m,l] * dl[0] * fe1arr[l]
            y_comp = tscaler * tauy_arr[m,l] * dl[1] * fe2arr[l]
            #u x- and y-components
            u_comp = uscaler * u_arr[m,l] * dl[0] * fe1arr[l]
            v_comp = uscaler * v_arr[m,l] * dl[1] * fe2arr[l]
            inc_tau += th_fac * (x_comp + y_comp)
            u_circ += th_fac * (u_comp + v_comp) 
        tau_circ += inc_tau * delta_t
        tau_circ_arr[m] = tau_circ
        u_circ_arr[m] = u_circ
    return tau_circ_arr, u_circ_arr
 
#function for solution of circulation for a given time and a given contour, direct integration used
def circ_solve(contour,zgrid,hgrid,udir,vdir,t_steps):
    darr = cont_depth(contour, zgrid) #depth array along contour
    fe1arr, fe2arr = cont_dist(contour, hgrid) #distance scales
    ufiles = glob.glob(udir + "/*.nc")[:t_steps] #start with 10 files
    vfiles = glob.glob(vdir + "/*.nc")[:t_steps] 
    o = cftime.datetime(2000,1,1,calendar="noleap") #temporal origin
    delta_t = 5*24*3600
    time_arr = np.zeros(t_steps)
    taux_arr = np.zeros((t_steps,len(contour)))
    tauy_arr = np.zeros((t_steps,len(contour)))
    #open each dataset only once
    for k in range(t_steps):
        ufile = xr.open_dataset(ufiles[k],use_cftime=True)
        ufile = ufile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
        vfile = xr.open_dataset(vfiles[k],use_cftime=True)
        vfile = vfile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
        time_arr[k] = (ufile["time_centered"].data - o).total_seconds()
        ctaux, ctauy = cont_tau(contour, zgrid, ufile, vfile)
        taux_arr[k] = ctaux
        tauy_arr[k] = ctauy
        #get velocities from final dataset
        if k == t_steps - 1:
            print(time_arr[k])
            contu, contv, meanu, meanv = cont_vel(contour, zgrid, ufile, vfile)
        #close files, to hopefully save memory
        ufile.close()
        vfile.close()
    #all the necessary data has been extracted and created in above loop
    #now calculate
    #xcomp_arr = np.zeros((t_steps,len(contour)))
    #ycomp_arr = np.zeros((t_steps,len(contour)))
    tau_circ = 0
    u_circ = 0
    for l in range(len(contour)):
        x_sum = 0
        y_sum = 0
        for m in range(t_steps):
            scaler = (np.exp(R*time_arr[m] / darr[l]) * delta_t) / (rho_o * darr[l])
            x_sum += scaler * taux_arr[m,l]
            y_sum += scaler * tauy_arr[m,l]
        if l == len(contour) -1:
            nxt = contour[0]
        else:
            nxt = contour[l+1]
        dl = np.subtract(nxt,contour[l])
        #calculated circulation
        tau_circ += th_fac * ((x_sum * dl[0] * fe1arr[l]) + (y_sum * dl[1] * fe2arr[l]))
        #circulation with model u_o-vals - for comparison
        u_circ += np.exp(R*time_arr[-1] / darr[l]) * th_fac * ((meanu[l] * dl[0] * fe1arr[l]) + (meanv[l] * dl[1] * fe2arr[l]))
    return tau_circ, u_circ
        
#try a variant of circ_solve where the order of summation is switched - potentially much faster to run.
def circ_solve_var(contour,zgrid,hgrid,udir,vdir,t_steps):
    #contx = xr.dataArray(contour[:][0], dims = "nz")
    #conty = xr.dataArray(contour[:][1], dims = "nz")
    darr = cont_depth(contour, zgrid) #depth array along contour
    fe1arr, fe2arr = cont_dist(contour, hgrid) #distance scales
    ufiles = glob.glob(udir + "/*.nc")[:t_steps] #start with 10 files
    vfiles = glob.glob(vdir + "/*.nc")[:t_steps] 
    o = cftime.datetime(2000,1,1,calendar="noleap") #temporal origin
    delta_t = 5*24*3600
    time_arr = np.zeros(t_steps)
    taux_arr = np.zeros((t_steps,len(contour)))
    tauy_arr = np.zeros((t_steps,len(contour)))
    u_arr = np.zeros((t_steps,len(contour)))
    v_arr = np.zeros((t_steps,len(contour)))
    #open each dataset only once
    for k in range(t_steps):
        print(k)
        ufile = xr.open_dataset(ufiles[k],use_cftime=True)
        ufile = ufile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
        vfile = xr.open_dataset(vfiles[k],use_cftime=True)
        vfile = vfile.isel(time_counter = 0, x = slice(0,1580,th_fac), y = slice(0,1801,th_fac))
        time_arr[k] = (ufile["time_centered"].data - o).total_seconds()
        #get tau_values
        print("lager tau")
        ctaux, ctauy = cont_tau(contour, zgrid, ufile, vfile)
        taux_arr[k] = ctaux
        tauy_arr[k] = ctauy
        #get velocities    
        print("lager u")
        meanu, meanv = cont_vel(contour, zgrid, ufile, vfile)
        u_arr[k] = meanu
        v_arr[k] = meanv
        #close files, to hopefully save memory
        print("lukker filer")
        ufile.close()
        vfile.close()
    #all the necessary data has been extracted and created in above loop
    #now calculate
    tau_circ_arr = np.zeros(t_steps)
    u_circ_arr = np.zeros(t_steps)
    tau_circ = 0
    #u_circ = 0
    for m in range(t_steps):
        inc_tau = 0 #tau circulation increment
        u_circ = 0 #u circulation
        for l in range(len(contour)):
            #can calculate u_circulation here as well
            tscaler = np.exp(R*time_arr[m] / darr[l]) / (rho_o * darr[l])
            uscaler = np.exp(R*time_arr[m] / darr[l])
            if l == len(contour) -1:
                nxt = contour[0]
            else:
                nxt = contour[l+1]
            dl = np.subtract(nxt,contour[l])
            #tau x- and y-components
            x_comp = tscaler * taux_arr[m,l] * dl[0] * fe1arr[l]
            y_comp = tscaler * tauy_arr[m,l] * dl[1] * fe2arr[l]
            #u x- and y-components
            u_comp = uscaler * u_arr[m,l] * dl[0] * fe1arr[l]
            v_comp = uscaler * v_arr[m,l] * dl[1] * fe2arr[l]
            inc_tau += th_fac * (x_comp + y_comp)
            u_circ += th_fac * (u_comp + v_comp) 
        tau_circ += inc_tau * delta_t
        tau_circ_arr[m] = tau_circ
        u_circ_arr[m] = u_circ
    return tau_circ_arr, u_circ_arr
            
"""
    