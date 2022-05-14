# -*- coding: utf-8 -*-
"""
Created on Mon May  2 19:34:55 2022

@author: jmcti
"""
#This is an old two component model from section 1 of the pdf, only using numpy arrays. When compared to the
#newer 2 component code the speed is not comparable.

import numpy as np
import matplotlib.pyplot as plt
import time


tic = time.perf_counter()

#4th order forward discrete laplacian
def laplacian(Z):
    L = ((-1/12)*np.roll(Z,-2, axis = 0) + (4/3)*np.roll(Z,-1, axis = 0) + (4/3)*np.roll(Z,1, axis=0)-(1/12)*np.roll(Z,2, axis = 0)
           -(1/12)*np.roll(Z,-2, axis = 1) +(4/3)*np.roll(Z,-1, axis = 1) + (4/3)*np.roll(Z,1, axis=1)-(1/12)*np.roll(Z,2, axis = 1)
            -5*Z)
    return L
#Time evolution
def tderiv(arr):
    var_deriv = laplacian(arr)+1/delta**2 * (arr-arr**3)
    diff = (-k/tau * laplacian(var_deriv))
    return diff

#Runge Kutta 4th order time difference calculation
def tstep(u_0):
    k1 = tderiv(u_0)
    u_1 = u_0+k1*dt/2
    k2 = tderiv(u_1)
    u_2 = u_0+k2*dt/2
    k3 = tderiv(u_2)
    u_3 = u_0+k3*dt
    k4 = tderiv(u_3)
    arr = u_0 + 1/6 * dt * (k1+2*k2+2*k3+k4)
    return arr

#Initializing array such that it is between -1 and 1
def initialize(length,mu,sigma):
    while True:
        uarr = np.random.normal(mu, sigma, (length,length))
        if np.amax(uarr)<=1 and np.amin(uarr)>=-1:
            break
    return uarr

# Parameters matching that of the newer model
k = .1
delta = 1
tau=1
ntime = 100000
lat_length = 100
dt = .2
ini_dens = -.1
var = .1

u_arr=initialize(lat_length,ini_dens,var)
fig, ax = plt.subplots()

# Evolving the system through time
for i in range(ntime):
    u_arrn=tstep(u_arr)        
    u_arr=u_arrn
toc = time.perf_counter()
print(toc-tic)

#uncomment savefig command to save figure
def main(name):
    v = ax.imshow(u_arr, interpolation='gaussian', cmap='viridis')
    fig.colorbar(v, ax=ax)
    plt.title('%s, %s, %s' %(k/tau,delta,ini_dens))
    #fig.savefig(name, dpi = 1200)

main('ktau%s_delta%s_final.png')