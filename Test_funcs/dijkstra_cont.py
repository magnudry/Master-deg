#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:54:07 2021

@author: magnudry
"""
import numpy as np
import xarray as xr
"""
def dijkstra_q(grid, orig):
    #mark all nodes unvisited
    grid = grid.assign(visited = grid["fonh"]*False) #visited False is 0.0
"""    
#function for creating f/H-contours - index (i,j), not coordinates
def contour_q(grid, orig):
    #grid: the topography dataset
    #rang: the error margin range for a contour - positive number
    #orig: numpy array of spatial origin of the contour - index (i,j), not coordinates
    # want to make contours of equal values for fonh - must use a range for each contour.
    # don't know the length of the contour arrays - gonna be slow (use of append)
    x0 = orig[0]
    y0 = orig[1]
    fonharr = np.array([[x0,y0]])
    #check a 5x5 matrix about point in question
    #next point is the one with fonh-value closest to point in question
    while True:
        mindiff = 1e-4
        for i in range(x0-3,x0+4):
            for j in range(y0-3,y0+4):
                #print("hei")
                #don't check point itself or previous points
                pot = np.array([i,j])
                
                if pot in fonharr:
                    continue
                #print(pot)
                diff = np.abs(grid["fonh"].isel(x = orig[0], y = orig[1]) - grid["fonh"].isel(x = i, y = j))
                if diff < mindiff: #range requirement can break loop, use with care
                    mindiff = diff
                    nxt = np.array([i,j])
                    print(nxt)
                    #print("hei")
        #print(nxt)
        if nxt[0] == orig[0] and nxt[1] == orig[1]:
            break
        fonharr = np.append(fonharr,[nxt], axis = 0)
        #print(fonharr)
        x0 = nxt[0]
        y0 = nxt[1]
    return fonharr

def contour_new(grid, orig):
    #problem in contour_q is if whole environment about point is in fonharr, while-loop iterates forever
    #try to exclude the environment about former point in fonharr for each new point
    
    x0 = orig[0]
    y0 = orig[1]
    fonharr = np.array([[x0,y0]])
    #counter = 0
    envir = np.zeros((9,2),dtype = int) #array of the points surrounding previous contour point
    while True: 
        mindiff = 1 #arbitrary high number
        nenvir = np.zeros((9,2),dtype = int) #array of the points surrounding point in question (create in loop)
        counter = 0
        for i in range(x0-1,x0+2):
            for j in range(y0-1,y0+2):
                pot = np.array([i,j])
                nenvir[counter] = pot
                counter += 1
                if pot in envir:
                    continue
                #now check the difference in value between candidate points and value at origin to determine next
                #point in contour
                diff = abs(grid["fonh"].isel(x = orig[0], y = orig[1]) - grid["fonh"].isel(x = i, y = j))
                if diff < mindiff:
                    nxt = np.array([i,j]) #potential next point in fonharr
                    mindiff = diff #adjust minimum difference
        if nxt[0] == orig[0] and nxt[1] == orig[1]:
            break
        envir = nenvir
        fonharr = np.append(fonharr,[nxt], axis = 0)
        x0 = nxt[0]
        y0 = nxt[1]
    return fonharr
                
                
                
                
                
    