#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 22:15:17 2020

@author: rohit
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 00:09:39 2020

@author: rohit
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# If plot doesn't show, please try Tools>IPython Console>Graphics>Graphics backend>Qt5'

N=4
size=12
pause = 4 # frame buffer time
dt=.001
distance = .2
colormap = np.array(['k','b','g','y','r'])
status = np.zeros((N,), dtype=int) # colormap to keep track of vorticity direction

I = np.identity(2) # identity tensor

Fc = [0,-1] # basic force unit vector

constant = 25/(8*np.pi) # scalar component of force should be multiplied here
def stokes(r):
    RR= np.outer(r,r)
    rmag = np.sqrt(np.dot(r,r))
    #Fc = F/np.abs(F)
    if np.dot(r,r)>.5: # boundary condition
        return  constant*(np.dot(Fc,I)/rmag+np.dot(Fc,RR)/rmag**3)
    else:
        return  constant*(.5*np.dot(Fc,I)+.5*np.dot(Fc,RR)) # boundary condition
    
xx,yy,xm,ym = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
x_v,y_v = np.zeros(N),np.zeros(N)

x = np.array([5.5,6.5,5.5,6.5])
y = np.array([8.5,8.5,9.5,9.5])
vorticity = np.zeros(N)
def honey(x,y):
    R =[]
    Rv = []
    for i in range(len(x)):
        R.append([x[i],y[i]])
        Rv.append([x_v[i],y_v[i]])
    for j in range(len(x)):
        for k in range(len(y)):
            r = np.array(R[j]) - np.array(R[k])
            Rv[j]+= stokes(r)
    return np.array(Rv)



def particles(x,y,x_v,y_v):
    for i in range(N):
        j=N-1
        while j>i: 
            if (np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)<distance):
                  a,b = x_v[i],y_v[i]
                  x_v[i] = x_v[j]
                  x_v[j] = a
                  y_v[i] = y_v[j]
                  y_v[j] = b        
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
        x[i] += dt*x_v[i] # new pos is velocity times small dt added to initial pos
        y[i] += dt*y_v[i]
       
    return x,y,x_v,y_v


# R = np.array([xx[0],yy[0]])
# plt.scatter(x,y)
# RR = honey(x,y)
# xv,yv=np.hsplit(RR,2)
# plt.quiver(x,y,xv,yv)

def animate(i):
    for i in range(pause):
        RR = honey(x,y)
        xv,yv = np.hsplit(RR,2)
        xx,yy,vx,vy=particles(x,y,xv,yv) # You can change this to stepper_tangible for bouncing!

    plt.cla()
    plt.axis([0,size,0,size])
    plt.title(str(N) + ' Particle simMulAtion')
    plt.scatter(xx,yy,c=colormap[status])
  
        
anim = ani.FuncAnimation(plt.gcf(),animate,interval=20)
plt.show()