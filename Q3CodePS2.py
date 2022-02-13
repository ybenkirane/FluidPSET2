# -*- coding: utf-8 -*-
"""
Insert description of the code

@autho: Yacine Benkirane
@collab: Jasmine Zhang and Alexandre Stuart
Created on Sat Feb 8 15:23:02 2022
"""

import numpy as np
import matplotlib.pyplot as plt
import imageio as imo

dT = 0.1 #the timestep of our simulation
Nsteps = 60 #total number of timesteps

## Setting up initial conditions (vortex centres and circulation)
# Vortex rings

y_v = [-1,1,-1,1]
x_v = [-1,-1,1,1]
k_v = [-1,1,-1,1]


#Setting up the plot 

plt.ion()
fig,ax = plt.subplots(1,1)


#let's mark vortex coordinates

p, = ax.plot(x_v, y_v, '.m', markersize=10)


# draw the initial velocity streamline
ngrid = 4 #dimension of your simulation grid
Y, X = np.mgrid[-ngrid:ngrid:1080j, -ngrid:6:1080j]
#360j sets the resolution of the cartesian grid; play around with it as you see fit

velocity_x = np.zeros(np.shape(X))#this holds x−velocity
velocity_y = np.zeros(np.shape(Y))#this holds y−velocity


#for loop instancing through each vortex element
for n in range(len(x_v)): 
    
    #computing the velocity field
    for m in range(len(velocity_x[0])):         #columns
        for l in range(len(velocity_x)):        #rows
            diff_y = y_v[n] - Y[l][m]
            diff_x = x_v[n] - X[l][m]
            
            #integrating advection velocity into vortex
            if (y_v[n] - Y[l][m])**2 + (x_v[n] - X[l][m])**2 !=0:
                
                #ensuring vortex doesn't advecate itself
                velocity_x[l][m]+=k_v[n]*diff_y/(diff_x**2+diff_y**2)
                velocity_y[l][m]+=-k_v[n]*diff_x/(diff_x**2+diff_y**2)       

#Boundary and initial conditions
ax.set_xlim([-ngrid, 6])
ax.set_ylim([-ngrid, ngrid])    

plt.streamplot(X,Y,velocity_x,velocity_y,density=[1,1],color='g')
plt.savefig('HDframe0.png')

#Time Evolution of Vortex Interaction

tally = 0

while tally < Nsteps:
    
    #initial velocity conditions
    velocity_x_vortex=[0,0,0,0]
    velocity_y_vortex=[0,0,0,0]
    
    #sum velocity vectors
    for n in range(len(x_v)):
        for m in range(len(x_v)):
            diff_y = y_v[m] - y_v[n] 
            diff_x = x_v[m] - x_v[n]
            if diff_y**2+diff_x**2 !=0:
                velocity_x_vortex[n]+=k_v[m]*diff_y/(diff_x**2+diff_y**2)
                velocity_y_vortex[n]+=-k_v[m]*diff_x/(diff_x**2+diff_y**2) 
                
    #Position update post-frame
    for n in range(len(x_v)):    
        x_v[n]=x_v[n]+velocity_x_vortex[n]*dT
        y_v[n]=y_v[n]+velocity_y_vortex[n]*dT

    #post-frame velocity update
    velocity_x = np.zeros(np.shape(X))
    velocity_y = np.zeros(np.shape(Y))
    
    #Stream-line Velocity Initializing
    for n in range(len(x_v)):
        for m in range(len(velocity_x[0])):
            for l in range(len(velocity_x)):
                diff_y = y_v[n] - Y[l][m]
                diff_x = x_v[n] - X[l][m]
                if (y_v[n]-Y[l][m])**2 + (x_v[n]-X[l][m])**2 !=0:
                    velocity_x[l][m]+=k_v[n]*diff_y/(diff_x**2+diff_y**2)
                    velocity_y[l][m]+=-k_v[n]*diff_x/(diff_x**2+diff_y**2)
    
    ##Plot Update
    #End of Streamline Frame
    ax.collections = []
    ax.patches = []
    
    p.set_xdata(x_v)
    p.set_ydata(y_v)
    
    #Plotting Position and velocity streamlines 
    plot=ax.streamplot(X,Y,velocity_x,velocity_y,density=[1,1],color='g')
    plt.savefig('HDframe'+ str(tally+1) +'.png')
    
    #tallying to next frame before loop reset
    tally+=1

# Building GIF of all frames combined 
        
with imo.get_writer('VortexInteractionAnimation.gif', mode='I') as writer:
    for filename in ['HDframe'+ str(tally) +'.png' for tally in range(100)]:
        image = imo.imread(filename)
        writer.append_data(image)