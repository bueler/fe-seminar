#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#importing packages
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import brentq
from scipy.fft import rfft, rfftfreq

#reading data from CSV
Data = pd.read_csv('ReceiverData')
Time = pd.read_csv('ReceiverTime')

#turning them into numpy float arrays
Data=(Data.columns.values).astype(float)
Time=(Time.columns.values).astype(float)

#calculating a fft
N = len(Data)
yf = rfft(Data)
xf = rfftfreq(N, 1 / (1/0.001)) #timestep is 0.001 

#resonate frequency for a spherical shell cavity
def func_resfreq(w):
    a=10.               #inner radius
    b=10.+25.*0.1       #outer radius
    
    #trancendental equation from Jackson EM problem 8.7
    return (w**2+1/(b*a))/(w**2+(a*b*(w**2-1/a**2)*(w**2-1/b**2)))-np.tan(w*(b-a))/(w*(b-a))

#function for finding roots of given function
#credit goes to GertVdE on stack exchange
#https://scicomp.stackexchange.com/questions/27314/how-to-solve-the-transcendental-equation-tanx-frac2xx2-1
def roots(N):
    # Find N roots of the equation

    #Allocate space
    roots = np.zeros(N)

    #Margin to stay away from poles
    margin = 1e-8

    #First root
    roots[0] = brentq(func_resfreq, 1.0 + margin, np.pi/2 - margin, rtol=1e-14)

    #Subsequent N-1 roots
    for i in range(1,N):
        left = (2*i - 1)*np.pi/2
        right = (2*i + 1)*np.pi/2
        roots[i] = brentq(func_resfreq, left + margin, right - margin, rtol = 1e-14)

    return roots 

#get the first 10 roots
result = roots(10)

#plotting
fig, (ax1, ax2) = plt.subplots(2,dpi=300)

#plot the time series
ax1.plot(Time, Data)

#plot the fft
ax2.plot(xf,np.abs(yf))

#plot suspected resonate frequencies
ax2.axvline(np.pi*np.sqrt(1/(25.*0.1)**2),10**-3,10**3,c="grey",ls="-.")
ax2.axvline(np.pi*np.sqrt(4/(25.*0.1)**2),10**-3,10**3,c="grey",ls="-.")
ax2.axvline(np.pi*np.sqrt(9/(25.*0.1)**2),10**-3,10**3,c="grey",ls="-.")
ax2.axvline(np.pi*np.sqrt(16/(25.*0.1)**2),10**-3,10**3,c="grey",ls="-.")
ax2.axvline(np.pi*np.sqrt(25/(25.*0.1)**2),10**-3,10**3,c="grey",ls="-.")
ax2.axvline(np.pi*np.sqrt(36/(25.*0.1)**2),10**-3,10**3,c="grey",ls="-.")
ax2.axvline(result[0],10**-3,10**3,c="k",ls=":")
ax2.axvline(result[1],10**-3,10**3,c="k",ls=":")
ax2.axvline(result[2],10**-3,10**3,c="k",ls=":")

#labeling and axis parameters
ax2.set_ylim(10**-2,10**2)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax1.set_ylabel(r"$\phi$")
ax1.set_xlabel("Time")
ax2.set_ylabel(r"$\phi^2$/Hz")
ax2.set_xlabel("Hz")
plt.show()