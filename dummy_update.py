#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 13:51:53 2022

@author: magnudry
"""
#update the dummy data test to suit the var-functions in solution_lin and the new routine for circulation calculations
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from solution_lin import *
#%%
#constants and arrays
#Circular bathymetric contours about North Pole, and circular wind stress on top.
#quadratic dataArrays for bathymetry and wind
H0 = 5000 #max depth
zgrid = np.zeros((99,99)) #the bathymetry grid
mid = 49 #midpoint of grid (north pole)
dis = 3000 #prescribed distance between points in grid
test_r = 1e5
t_steps = 150
ix = np.arange(99)
iy = np.arange(99)
#Make the meshgrid and the corresponding radius and angle arrays
x = np.linspace(-49,49,99)
y = np.linspace(-49,49,99)
xs, ys = np.meshgrid(x, y, sparse=True)
r_rad = np.sqrt(xs**2 + ys**2) * dis #raw radius
rad = np.sqrt(xs**2 + ys**2) / (np.sqrt(2)*49) #normalised radius
alpha = np.arctan2(ys,xs)
#make a circular wind stress grid on top
#the wind stress should oscillate sinusoidally with time 
xgrid = np.zeros(((99,99,t_steps))) #xtau-grid
ygrid = np.zeros(((99,99,t_steps))) #ytau-grid
omega_s = 6 * np.pi / (365*24*3600) #frequency for the surface stress
tmax = 0.1 #max wind-stress
tmin = 0.25 #min wind-stress
dt = 24*3600 #time step
times = np.arange(t_steps) * dt
days = np.arange(t_steps)
#dt2 = 2.5*24*3600
#%%
#print(rad)
print(times[-1])
#%%
plt.pcolormesh(alpha, cmap = "viridis")
#%%
print(alpha[:,48])
#%%
#Grid functions
#function for depth at given radius
def bell_h(r): #r is a 2D array or scalar of the radius from the NP
    H = H0 * np.exp(-r**2/(7.2e9))
    return H 
#function for the stress at given time
def stress_mg(r,angle):
    #tau = tmax * r * np.sin(omega_s * t)
    tau = tmax * r
    tx = tau * np.sin(angle)
    ty = -1 * tau * np.cos(angle)
    return tx, ty
def stress(gx,gy,r,angle): #using func_stress on every point in the tau-arrs, and return filled arrays
    for c in range(t_steps):
        time = times[c]
        gx[:,:,c],gy[:,:,c] = stress_mg(r,angle,time)
    return gx,gy
#function for creating the constant depth contours of the dummy data set
def cont_H(grid,x,y,r):
    level = bell_h(r)
    contdat = plt.contour(x,y,grid,levels=[level])
    contour = contdat.allsegs[0][0]
    return contour
#simplified version of cont_data_var, no need for cont dist and cont depth
def cont_data_d(contour, txfile, tyfile):
    contx = xr.DataArray(contour[:,0], dims = "nz")
    conty = xr.DataArray(contour[:,1], dims = "nz")
    conttx, contty = cont_tau_var(contx,conty,txfile,tyfile)
    tau_data = xr.Dataset(data_vars = dict(
                            conttx = conttx,
                            contty = contty))
    """
    first = True
    for k in range(t_steps):
        conttx, contty = cont_tau_var(contx,conty,txfile.isel(time = k),tyfile.isel(time = k))
        curr_set = xr.Dataset(data_vars = dict(
            conttx = conttx,
            contty = contty))
        curr_set = curr_set.assign_coords(t_step = k) #zero indexed
        curr_set = curr_set.expand_dims("t_step")
        if first:
            tau_data = curr_set
            first = False
        else:
            tau_data = xr.concat([tau_data,curr_set],dim = "t_step")
    """
    return tau_data, contx, conty
#dummy circ_solve
def circ_solve_d(tau_d, t_arr, H, contx, conty):
    first = True
    for i in range(len(t_arr)):
        t = t_arr[i]
        norm = np.exp(-R * t / H) / (2 * np.pi * test_r)
        #constant stress
        stress_da = (tau_d["conttx"][1:] * contx.diff(dim = "nz") * dis + tau_d["contty"][1:] * conty.diff(dim = "nz") * dis)
        #sinusoidal stress
        #stress_da = (tau_d["conttx"].isel(time = i)[1:] * contx.diff(dim = "nz") * dis + tau_d["contty"].isel(time = i)[1:] * conty.diff(dim = "nz") * dis)
        temp_tau = (1/(rho_o * H)) * np.exp(R * t / H) * dt * stress_da
        if first:
            stress_set = stress_da.sum(dim = "nz") / (2 * np.pi * test_r)
            #tau_iter = temp_tau.sum(dim = "nz")
            #tau_circ = temp_tau.sum(dim = "nz") * norm
            #test initialisation
            tau_iter = 0
            tau_circ = temp_tau.sum(dim = "nz") * 0
            first = False
        else:
            stress_set = xr.concat([stress_set,stress_da.sum(dim = "nz") / (2 * np.pi * test_r)], dim = "t_step")
            tau_iter += temp_tau.sum(dim = "nz")
            tau_circ = xr.concat([tau_circ,tau_iter * norm], dim = "t_step")
    return tau_circ, stress_set    
#Analytical function
def an_circ(r,t,H):
    #H = bell_h(r)
    r_n = r / (mid * dis)
    #print(H)
    #K1 = np.exp(-R*t/H) * 1/(rho_o * H)
    K1 = (2 * np.pi * r) / (rho_o * H)
    K2 = ((tmax - tmin)/2) * r_n /((R/H)**2 + omega_s**2)
    K3 = (tmax/2) * (H/R)
    #print(K1,K2,K3)
    lhs = K1 * (K2 * (np.exp(R*t/H) * ((R/H)*np.sin(omega_s * t) - omega_s * np.cos(omega_s * t)) + omega_s) + K3 * (np.exp(R*t/H) - 1))
    return lhs

#new analytical sol
def an_circ_n(r,t,H):
    #for tau = Ar
    lhs = (np.exp(R * t / H) - 1) * 2 * np.pi * tmax * (r**2) / (np.sqrt(2) * mid * dis * rho_o * R) 
    #for tau = Ar sin(wt)
    #lhs = (np.exp(R * t / H) * ((R/H)*np.sin(omega_s*t) - omega_s*np.cos(omega_s*t)) + omega_s) * 2 * np.pi * tmax * (r**2) / (np.sqrt(2) * mid * dis * rho_o * H * ((R/H)**2 + omega_s**2)) 
    #lhs = 2 * np.pi * (tmax / (mid * dis)) * r**2 / (rho_o * R) * (np.exp(R * t / H)-1)
    return lhs

#%%
#make depth array
depths = bell_h(r_rad)
xtau,ytau = stress_mg(rad, alpha)
#xtau,ytau = stress(xgrid,ygrid,rad,alpha)
#%%
print(depths[49,49])
#print(xtau[45,16])
#%%
#make datasets in xarray of the relevant variables - to test the functions of solution_lin

txfile = xr.Dataset(data_vars = dict(
        sozotaux = (["y","x"], xtau)), 
    attrs = dict(
        description = "dummy tau_x variable"))
tyfile = xr.Dataset(data_vars = dict(
        sometauy = (["y","x"], ytau)), 
    attrs = dict(
        description = "dummy tau_y variable"))
zfile = xr.Dataset(data_vars = dict(
        fonh = (["y","x"], depths)),
    attrs = dict(
        description = "Phony fonh-values grid - pure topography"))
#%%
print(txfile)
#%%
cont = cont_H(zfile["fonh"],ix,iy,test_r)
H = bell_h(test_r)
#%%
print(H)
#%%
cont_tau, contx, conty = cont_data_d(cont, txfile, tyfile)
#%%
print(cont_tau["conttx"])
print(contx.diff(dim = "nz"))
#%%
circ_tau, stress_arr = circ_solve_d(cont_tau, times, H, contx, conty)
#%%
print(np.abs(circ_tau))
#%%
an_set = np.zeros(t_steps)
for i in range(t_steps):
    norm = np.exp(-R * times[i] / H) / (2 * np.pi * test_r)
    an_set[i] = an_circ_n(test_r, times[i], H) * norm
#%%
#norm = np.exp(-R * times[i] / H) / (2 * np.pi * test_r)
#%%
fig, axs = plt.subplots(2,1, sharex = True, gridspec_kw = {"height_ratios" : [3,1]}, figsize=(8,6))
#fig.suptitle("Circulation at radius 100 km, depth 1246.8 m",fontsize = "xx-large")
#axs[0].set_title("Mean velocity about contour [m/s]",fontsize = "x-large")
axs[0].plot(days,-1 * circ_tau, label = "Numerical")
axs[0].plot(days,an_set, label = "Analytical")
axs[0].set_ylabel(r"$u_{m}$ [m/s]", fontsize = 20)
axs[0].tick_params(labelsize = 15)
axs[0].grid(True)
#axs[0].legend()
#axs[1].set_title("Surface stress [Pa]",fontsize = "x-large")
axs[1].plot(days,-1*stress_arr)
axs[1].set_xlabel("Days",fontsize = 20)
axs[1].set_ylabel(r"$\tau_{m}$ [Pa]", fontsize = 20)
axs[1].tick_params(labelsize = 15)
axs[1].grid(True)
fig.tight_layout(rect = (0,0,1,1))
#%%
folder = "/uio/lagringshotell/geofag/students/metos/magnudry/Master/Figures/"
name = "dummy_lin_sol_tinit.pdf"
fig.savefig(("{}{}".format(folder,name)))
#%%
print(an_set)
#%%
print((-1 * circ_tau)/an_set)
