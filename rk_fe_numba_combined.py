# -*- coding: utf-8 -*-
"""
Created on Fri May 13 11:51:10 2022

@author: jmcti
"""
#With this code we are simulating the differential equation shown in Section 1 of the Phys 230 Project pdf. 
#Specifically, this is liquid-liquid phase separation for two components on a flat surface with a density array ranging from 
#from -1 to 1. We also want to compare the simulation time between 4th order Runge Kutta and Forward Euler.

#Importing relevant packages
import numpy as np
import matplotlib.pyplot as plt
from numba import njit #Used to speed up loops
import time


tic = time.perf_counter() #To find initial time of simulation

#See Section 1 of Phys 230 Project pdf
k = .1
delta = 1
tau=1

ntime = 4000000 #Number of time steps
lat_length = 100 #Lattice ranges from 0-lat_length in both x and y directions 
dt = .1 #Time step

ini_dens = -.5 #Initial average density of array
var = .05 #Standard deviation used to initialize array

graph='' # 'yes' specifies if you want graphs shown in interactive window throughout simulation
g_iter= 1000 #Specifies the number of time steps between showing graphs
tcls = '' #Specifies if you want to generate an array corresponding to simulation time at various steps
time_vals = [] # Array to hold simulation time at various time steps
tcalc_iter=10000 #Number of steps between each simulation time check 

#Use 4th order forward difference for laplacian with periodic boundary conditions
@njit(fastmath=True)
def laplacian(Z):
    L=np.zeros((lat_length,lat_length)) #Initializing array for laplacian
    for i in range(lat_length):
       iup1 = (i+1)%lat_length #These relations are used to enforce PBC, since going over lattice length goes back to 1
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

#Finding time derivative for density array
@njit(fastmath=True)
def tderiv(arr):
    var_deriv = laplacian(arr)+1/delta**2 * (arr-arr**3)
    diff = (-k/tau * laplacian(var_deriv))
    return diff

#Using Forward Euler method to evolve the system forward a time step
@njit(fastmath=True)
def tstep(u_0):
    arr = u_0 + tderiv(u_0)*dt
    return arr

#Using 4th order Runge Kutta to evolve system forward a time step
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

#Initializing a single array 
@njit(fastmath=True)
def initialize(length,mu,sigma):
    while True:
        uarr = np.random.normal(mu, sigma, (length,length))
        if np.amax(uarr)<=1 and np.amin(uarr)>=-1: #Making sure array fits within -1 and 1
            break
    return uarr

#Creating start function that evolves the system through time. Can uncomment/comment tstep(u_arr)
#to switch between the Runge Kutta and Forward Euler methods
def start(u_arr,gcls,tcalc):
    ngraph=0 #Numbering for printed out graphs if needed
    for i in range(ntime):
        u_arrn=tstep(u_arr)
        #u_arrn=tsteprk(u_arr)
        if np.isnan(u_arrn[0,0])==True: #Making sure simulation doesn't break
            print('Broken: Time step needs to be smaller')
            raise SystemExit
        if i%g_iter==0 and gcls=='yes': #If wanted this generates a graph for each number of time steps. Can uncomment fig.savefig command to save image
            ngraph+=1
            fig, ax = plt.subplots()
            v=ax.imshow(u_arrn, interpolation='gaussian', cmap='viridis')
            plt.colorbar(v,ax=ax)
            plt.title('u with: k/tau=%s, delta=%s, initial density=%s, time=%s' %(k/tau,delta,ini_dens,dt*i))
            #fig.savefig('ktau%s_delta%s_inidens%s_gnum%s.png' %(k/tau,delta,ini_dens,ngraph),dpi=1200)
        if tcalc=='yes' and i%tcalc_iter==0: #If wanted, this tells the code to determine simulation time throughout simulation
            toc = time.perf_counter()
            time_vals.append(toc-tic)
        u_arr=u_arrn
    return u_arr

#Creating main function to run simulation
def main():
    u_arr=initialize(lat_length,ini_dens,var)
    u_arrf=start(u_arr,graph,tcls)
    return u_arrf

final = main() 

#Plotting final array, can save by uncommmenting savefig command.
fig, ax = plt.subplots()
v=ax.imshow(final, interpolation='gaussian', cmap='viridis')
plt.colorbar(v,ax=ax)
plt.title('u with: k/tau=%s, delta=%s, initial density=%s, time=%s' %(k/tau,delta,ini_dens,dt*ntime))
#fig.savefig('ktau%s_delta%s_inidens%s_final.png' %(k/tau,delta,ini_dens),dpi=1200)
toc = time.perf_counter()
print('Simulation time array:',time_vals) #Print out array of simulation time values if needed
print('The total time of simulation was:',toc-tic)