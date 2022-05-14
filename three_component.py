# -*- coding: utf-8 -*-
"""
Created on Thu May 12 23:39:32 2022

@author: jmcti
"""
# With this code we are modeling liquid-liquid phase separation for a three component system based off the 
# Flory-Huggins free energy density shown in Section 5.1/6.3 (ui goes from 0-1). This is a simpler code just meant for finding the
# final configuration of our system for the given parameters. 

#Start by importing the relevant packages
import numpy as np
import matplotlib.pyplot as plt
from numba import njit #Used to make our loops faster
import time


tic = time.perf_counter()  #This represents the simulation start time

#See equation in Phys 230 Project pdf Section 5.1/6.3
D = .1  
lmbda = 1
chi12 = 4
chi13 = 4
chi23 = 4

#Below are used to initialize densities
ini_dens1 = .25 #Average initial density for u1
ini_dens2 = .25 #Average initial density for u2
var1 = .05  #Variation corresponding to normal distributions
var2 = .05

ntime = 1000000  #Number of time steps
lat_length = 100  #Lattice ranges from 0-lat_length in both x and y directions
dt = .005  #Size of time step


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

#Finding the time derivative for both density arrays
@njit(fastmath=True)
def tderiv(arr):
    term1 = -2*chi13*arr[0]+(chi12-chi13-chi23)*arr[1]
    term2 = (chi12-chi23-chi13)*arr[0]-2*chi23*arr[1]
    var_deriv1 = np.log(arr[0]/(1-arr[0]-arr[1]))+term1+lmbda**2*laplacian(term1)
    var_deriv2 = np.log(arr[1]/(1-arr[0]-arr[1]))+term2+lmbda**2*laplacian(term2)
    diff1 = (D*laplacian(var_deriv1))
    diff2 = (D*laplacian(var_deriv2))
    return diff1,diff2

#Using Forward Euler method
@njit(fastmath=True)
def tstep(u_0):
    time_deriv = tderiv(u_0)
    return u_0[0] + time_deriv[0]*dt,u_0[1] + time_deriv[1]*dt

#Initializing our two density arrays based off a normal distribution such that u1+u2+u3=1 (u3=1-u1-u2)
@njit(fastmath=True)
def initialize(length,mu1,sigma1,mu2,sigma2):
    while True:
        u1arr = np.random.normal(mu1, sigma1, (length,length))
        u2arr = np.random.normal(mu2,sigma2,(length,length))
        if np.amin(u1arr)>=0 and np.amin(u2arr)>=0 and np.amax(u1arr+u2arr)<=1: #Making sure arrays fit with ui definitions
            break
    return u1arr,u2arr

#Looping over our time steps. The njit command does not impact simulation time in this location.
def start(u_arr):
    for i in range(ntime):
        u_arrn=tstep(u_arr)
        if np.isnan(u_arrn[0][0][0])==True: #Making sure code doesn't fail mid simulation
            print('Broken: Time step needs to be smaller')
            raise SystemExit
        u_arr=u_arrn
    return u_arr

#Main function to start simulation
def main():
    u_arr=initialize(lat_length,ini_dens1,var1,ini_dens2,var2)
    u_arrf=start(u_arr)
    return u_arrf

#Determine evolved density arrays
final = main()

#Plotting arrays. Can uncomment fig.savefig command to generate high resolution final images
fig, ax = plt.subplots()
ax.imshow(final[0], interpolation='gaussian', cmap='viridis')
plt.title('u1 with: D=%s, lambda=%s, chi12=%s, chi13=%s, chi23=%s, u1 initial=%s, u2 initial=%s, time=%s' %(D,lmbda,chi12,chi13,chi23,ini_dens1,ini_dens2,ntime*dt))

#fig.savefig('u1_D%s_lambda%s_chi12%s_chi13%s_chi23%s_u1dens%s_u2dens%s_final.png' %(D,lmbda,chi12,chi13,chi23,ini_dens1,ini_dens2),dpi=1200)

fig, ax = plt.subplots()
ax.imshow(final[1], interpolation='gaussian', cmap='viridis')
plt.title('u2 with: D=%s, lambda=%s, chi12=%s, chi13=%s, chi23=%s, u1 initial=%s, u2 initial=%s, time=%s' %(D,lmbda,chi12,chi13,chi23,ini_dens1,ini_dens2,ntime*dt))

#fig.savefig('u2_D%s_lambda%s_chi12%s_chi13%s_chi23%s_u1dens%s_u2dens%s_final.png' %(D,lmbda,chi12,chi13,chi23,ini_dens1,ini_dens2),dpi=1200)

toc = time.perf_counter()
print('The total time of simulation was:',toc-tic) #Determining final time of simulation