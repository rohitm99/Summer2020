#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 12:46:22 2020

@author: rohit
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# If plot doesn't show, please try Tools>IPython Console>Graphics>Graphics backend>Qt5'

N=10
size=10
pause = 1 # frame buffer time
dt=.004
distance = .25 # distance for interaction
quarantine = 0

xx,yy,xm,ym = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
# randomly places N points within plot
def pos():
    for i in range(N):
        xx[i]+= np.random.random()*size
        yy[i]+= np.random.random()*size
    return xx,yy


def velocity(): # gives random initial velocity
    for i in range(N-quarantine):
        xm[i]=((np.random.random())-np.random.random())*1
        ym[i]=((np.random.random())-np.random.random())*1
    return xm,ym

x, y = pos() #randomized initial positions
x_v, y_v = velocity() #x_v and y_v are respective initial velocities

def particles(x,y,x_v,y_v):
    for i in range(N-quarantine):
        j=N-1-quarantine
        while j>i: 
            if (np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)<distance):
                  a,b = x_v[i],y_v[i]
                  x_v[i] = x_v[j]
                  x_v[j] = a
                  y_v[i] = y_v[j]
                  y_v[j] = b        
            else:
                x,y,x_v,y_v = stepper_tangible(x, y, x_v, y_v) 
            j-=1
        i+=1
    return x,y,x_v,y_v
    
def stepper_tangible(x,y,x_v,y_v): # this is a box with tangible walls
    for i in range(N-quarantine):  
        x[i] += dt*x_v[i] 
        y[i] += dt*y_v[i]
        # for a finite box with tangible edges(bouncing)
        if (x[i]) >= size or x[i]<=0:
            x_v[i] = -x_v[i]
            x[i] += x_v[i]*dt
            
        if (y[i]) >= size or y[i]<0:
            y_v[i] = -y_v[i]
            y[i] += y_v[i]*dt
            
    return x,y,x_v,y_v
   
def animate(i):
    for j in range(pause):
        xx,yy,xv,yv=particles(x,y,x_v,y_v) # You can change this to stepper_tangible for bouncing!
    plt.cla()
    plt.axis([0,size,0,size])
    plt.title(str(N) + ' Particle simulation')
    plt.scatter(xx,yy)
  
        
anim = ani.FuncAnimation(plt.gcf(),animate,interval=20)
plt.show()