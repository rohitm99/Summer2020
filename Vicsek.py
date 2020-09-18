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

N=30
size=5
pause = 1 # frame buffer time
dt=.0005
distance = .35 # distance for interaction


xx,yy,xm,ym,orientation = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
# randomly places N points within plot
def pos():
    for i in range(N):
        xx[i]+= np.random.random()*size
        yy[i]+= np.random.random()*size
    return xx,yy


def velocity(): # gives random initial velocity
    for i in range(N):
        xm[i]=((np.random.random())-np.random.random())*.5
        ym[i]=((np.random.random())-np.random.random())*.5
        orientation[i] = np.sin(np.random.random()-np.random.random())
    return xm,ym,orientation

x, y = pos() #randomized initial positions
x_v, y_v,orientation = velocity() #x_v and y_v are respective initial velocities

def particles(x,y,x_v,y_v):
    for i in range(N):
        j=N-1
        while j>i: 
            if (np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)<=distance):
                   x_v[i] = (.7*x_v[i]+.3*x_v[j])
                   y_v[i] = (.7*y_v[i]+.3*y_v[j])  
            else:
                x,y,x_v,y_v = stepper(x, y, x_v, y_v) 
            j-=1
        i+=1
    return x,y,x_v,y_v
    
def stepper(x,y,x_v,y_v): # intangible continous edges
    for i in range(N):
        if x[i] > size:
            x[i] -= size
        if y[i] > size:
            y[i] -=size 
        if x[i]<0:
            x[i]+=size
        if y[i]<0:
            y[i]+=size
        if np.random.random()>.99996:
            orientation[i]+=.5
            x_v[i]+=.09*np.cos(orientation[i])
            y_v[i]+=.09*np.sin(orientation[i])
        x[i] += dt*x_v[i] # new pos is velocity times small dt added to initial pos
        y[i] += dt*y_v[i]
       
    return x,y,x_v,y_v

def animate(i):
    for j in range(pause):
        xx,yy,xv,yv=particles(x,y,x_v,y_v) # You can change this to stepper_tangible for bouncing!
    plt.cla()
    plt.axis([0,size,0,size])
    plt.title(str(N) + ' Particle simulation')
    #plt.scatter(xx,yy)
    plt.quiver(xx,yy,xv,yv)
  
anim = ani.FuncAnimation(plt.gcf(),animate,interval=20)
plt.show()
