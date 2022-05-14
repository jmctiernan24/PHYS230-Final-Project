# -*- coding: utf-8 -*-
"""
Created on Fri May 13 12:37:46 2022

@author: jmcti
"""
#This code is modeling liquid-liquid phase separation for a two component system with the Flory-Huggins free energy
#from section 6.2 of the Physics 230 Project pdf. In this case the density ranges from 0-1.

#This code is also identical to rk_fe_numba_combined.py, with slight modifications involving the parameters used, form of the time
#derivative and initialization of the array. So I will only comment regions that have changed.


import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import time


tic = time.perf_counter()

# See section 6.2 of the Physics 230 Project pdf.
D = .1
lmbda = 1
chi12 = 4

ntime = 4000000
lat_length = 100
dt = .01

ini_dens = .4
var = .05

graph='' # 'yes' if you want graphs shown in interactive window
g_iter= 1000 #time steps between each graph
tcls = '' #Specifies if you want to generate an array corresponding to simulation time
time_vals = [] # Used to determine simulation time at various time steps
tcalc_iter=1000 # time steps between each simulation time check

@njit(fastmath=True)
def laplacian(Z):
    L=np.zeros((lat_length,lat_length))
    for i in range(lat_length):
       iup1 = (i+1)%lat_length
       idown1 = (i-1)%lat_length
       iup2 = (i+2)%lat_length
       idown2 = (i-2)%lat_length
       for j in range(lat_length):
           jup1 = (j+1)%lat_length
           jdown1 = (j-1)%lat_length
           jup2 = (j+2)%lat_length
           jdown2 = (j-2)%lat_length
           L[i,j]+= (-Z[iup2,j]+16*Z[iup1,j]-30*Z[i,j]+16*Z[idown1,j]-Z[idown2,j] -Z[i,jup2]+16*Z[i,jup1]-30*Z[i,j]+16*Z[i,jdown1]-Z[i,jdown2])/12
           
    return L

#Slightly different variational derivative when compared to original two component model
@njit(fastmath=True)
def tderiv(arr):
    var_deriv = np.log(arr/(1-arr))-2*chi12*arr-2*lmbda**2*chi12*laplacian(arr)
    diff = (D*laplacian(var_deriv))
    return diff

@njit(fastmath=True)
def tstep(u_0):
    arr = u_0 + tderiv(u_0)*dt
    return arr

@njit(fastmath=True)
def tsteprk(u_0):
    k1 = tderiv(u_0)
    u_1 = u_0+k1*dt/2
    k2 = tderiv(u_1)
    u_2 = u_0+k2*dt/2
    k3 = tderiv(u_2)
    u_3 = u_0+k3*dt
    k4 = tderiv(u_3)
    arr = u_0 + 1/6 * dt * (k1+2*k2+2*k3+k4)
    return arr

#Only modification here is that the density only ranges from 0 to 1.
@njit(fastmath=True)
def initialize(length,mu,sigma):
    while True:
        uarr = np.random.normal(mu, sigma, (length,length))
        if np.amax(uarr)<=1 and np.amin(uarr)>=0:
            break
    return uarr

#Only modification is with the plot titles
def start(u_arr,gcls,tcalc):
    ngraph=0
    for i in range(ntime):
        u_arrn=tstep(u_arr)
        #u_arrn=tsteprk(u_arr)
        if np.isnan(u_arrn[0,0])==True:
            print('Broken: Time step needs to be smaller')
            raise SystemExit
        if i%g_iter==0 and gcls=='yes':
            ngraph+=1
            fig, ax = plt.subplots()
            v=ax.imshow(u_arrn, interpolation='gaussian', cmap='viridis')
            plt.colorbar(v,ax=ax)
            plt.title('u with: D=%s, lambda=%s, chi12=%s, initial density=%s, time=%s' %(D,lmbda,chi12,ini_dens,dt*i))
            #fig.savefig('D%s_lambda%s_chi12%s_inidens%s_time%s.png' %(D,lmbda,chi12,ini_dens,dt*i),dpi=1200)
        if tcalc=='yes' and i%tcalc_iter==0:
            toc = time.perf_counter()
            time_vals.append(toc-tic)
        u_arr=u_arrn
    return u_arr


def main():
    u_arr=initialize(lat_length,ini_dens,var)
    u_arrf=start(u_arr,graph,tcls)
    return u_arrf

final = main()
fig, ax = plt.subplots()
v=ax.imshow(final, interpolation='gaussian', cmap='viridis')
plt.colorbar(v,ax=ax)
plt.title('u with: D=%s, lambda=%s, chi12=%s, initial density=%s, time=%s' %(D,lmbda,chi12,ini_dens,dt*ntime))
#fig.savefig('D%s_lambda%s_chi12%s_inidens%s_final.png' %(D,lmbda,chi12,ini_dens),dpi=1200)
toc = time.perf_counter()
print(time_vals)
print('The total time of simulation was:',toc-tic)