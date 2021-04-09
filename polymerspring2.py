#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 15:43:35 2021

@author: rohit
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from numpy import random

# If plot doesn't show, please try Tools>IPython Console>Graphics>Graphics backend>Qt5'
N=4
size=10
pause = 4 # frame buffer time
dt=.01
distance = .2 # distance for elastic collision
Rr = 1 # radius of particle
I = np.identity(2) # identity tensor
Fg=np.array([0,-10]) # force of gravity

constant = 3/4 # scalar component of force should be multiplied here
def stokes1(r,Fc):
    RR= np.outer(r*Rr,r*Rr)
    rmag = np.sqrt(np.dot(r,r))*Rr
    if np.dot(r,r)>.5: # boundary condition
        return  constant*(np.dot(Fc,I)/rmag+np.dot(Fc,RR)/rmag**3) # returns non diemsnional velocity
    else:
        return  constant*(.5*np.dot(Fc,I)+.5*np.dot(Fc,RR)) # boundary condition

x_v,y_v = np.zeros(N),np.zeros(N)-1

# initial positions for stokes function
x = np.array([3.0,4.0,6.0,7.0])
y = np.array([9.5,8.0,8.0,9.5])

# initial positions used for spring dx
X = np.array(x)
Y = np.array(y)

def honey(x,y,x_v,y_v):
    R =[]
    Rv = [] 
    vorticity = np.zeros(N)
    ks = 10.0*Rr/Fg[1] # non dimensional linear spring constant
    kt = 2 # torsional spring constant
    for i in range(len(x)):
        R.append([x[i],y[i]])
        Rv.append([x_v[i],y_v[i]])
    for j in range(len(x)):
        for k in range(N):
            if k!=j: # so self interaction due to stokeslet is ignored
               # non dimensional spring force = ks*dd
                r = np.array(R[j]) - np.array(R[k])
                if np.abs(j-k)==1: # since spring interactions only with adjacent particles
                    d = np.sqrt((x[j]-x[k])**2+(y[j]-y[k])**2) # current distance between 2 particles
                    dd =d-np.sqrt((X[j]-X[k])**2+(Y[j]-Y[k])**2)  # difference from equilibrium
                    Fs = [ks*dd*(x[j]-x[k])/(Rr*d),ks*dd*(y[j]-y[k])/(Rr*d)]
                    theta = (np.arcsin((y[j]-y[k])/(Rr*d))) # angle between interacting particles in rad
                    d_theta = theta - np.arctan((Y[j]-Y[k])/(X[j]-X[k]))
                    w=np.array([.5*d_theta*np.sin(theta)*kt,-.5*d_theta*np.cos(theta)*kt])
                    
                else:
                    Fs=[0,0]
                Rv[j] += stokes1(r,[Fs[0],Fg[1]/-Fg[1]+Fs[1]])
                vorticity[j]+=Vorticity(r)
               
    
    Rv[3]+=w
    Rv[0]+= [-w[0],w[1]]
    return image(R,np.array(Rv),np.array(vorticity))

def Vorticity(r):
    vt = stokes1(r+[0,.05],Fg/-Fg[1])
    vb = stokes1(r+[0,-.05],Fg/-Fg[1])
    vr = stokes1(r+[.05,0],Fg/-Fg[1])
    vl = stokes1(r+[-.05,0],Fg/-Fg[1])
    
    duxdy = (vr[1]-vl[1])/.1
    duydx = (vt[0]-vb[0])/.1
    
    omega = .5*(duydx-duxdy)
    return omega
    
def vortfield(vorticity,r):
    
    vfield = []
    rmag = np.sqrt(np.dot(r,r))
    vfield.append(np.cross([0,0,8*np.pi*vorticity*Rr**3],r)/(8*np.pi*rmag**3))
    return vfield[0]
    

def image(RR,RV,vorticity): # adds the velocity resulting from the image particles "beneath" bottom
    R_image =[]
    for i in range(N):
         R_image.append([x[i],-y[i]])
    for j in range(N): # real particles
        for k in range(N): # image particles
                r = np.array(RR[j])-np.array(R_image[k])
                RV[j]+=stokes1(r,Fg/Fg[1])
                RV[j]+=random.normal(0,1)*np.sqrt(2*1.380649e-23*298/(.95*6*np.pi*dt*1))*(.95*6*np.pi*1)/-Fg[1] # brownian motion term
                #vorticity[j]+=Vorticity(r)
    return np.array(RV),np.array(vorticity)



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
        x[i] += dt*x_v[i] # new pos is velocity times small dt added to initial pos
        y[i] += dt*y_v[i]
       
    return x,y,x_v,y_v


def animate(i):
    for i in range(pause):
        rv,vort = honey(x,y,x_v,y_v)
        xv,yv = np.hsplit(rv,2)
        xx,yy,vx,vy=stepper(x,y,xv,yv) 
    plt.cla()
    plt.axis([0,size,0,size])
    XX = np.round((X[1]-xx[1]),5)
    plt.title('Horizontal displacement: ' + str(XX))
    plt.scatter(xx,yy)
    plt.quiver(xx,yy,0,vort)
    plt.plot(xx,yy)
    plt.grid()
    plt.show()
  
        
anim = ani.FuncAnimation(plt.gcf(),animate,interval=20)
plt.show()