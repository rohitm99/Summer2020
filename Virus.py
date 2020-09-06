#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 12:29:33 2020

@author: rohit@ArchLinux

"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# If plot doesn't show, please try Tools>IPython Console>Graphics>Graphics backend>Qt5'

N=20
size=16
pause = 1# BEWARE changing this can alter the effective rec_time
dt=.1
distance = 2 # distance for infection
rec_time = 50 # time until someone recovers, 2.25 units per second
initial_infected = 2
masks = int(np.round(N*0.9))
masked = np.zeros(N)

quarantine =int(np.round(N*0.2)) # adjust the fraction within to change percentage sheltered
colormap = np.array(['g','r','b'])

for i in range(masks-quarantine):
    masked[i+initial_infected]+=np.random.random()
    i+=1

status = np.zeros((N,), dtype=int)
for i in range(initial_infected):
    status[i]+=1
    i+=1

xx,yy,xm,ym,days = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
# randomly places N points within plot
def pos():
    for i in range(N):
        xx[i]+= np.random.random()*size
        yy[i]+= np.random.random()*size
    return xx,yy


def velocity(): # gives random initial velocity
    for i in range(N-quarantine):
        xm[i]=((np.random.random())-np.random.random())
        ym[i]=((np.random.random())-np.random.random())
    return xm,ym

x, y = pos() #randomized initial positions
x_v, y_v = velocity() #x_v and y_v are respective initial velocities
def stepper(x,y,x_v,y_v): # intangible continous edges
    for i in range(N-quarantine):
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

def stepper_tangible(x,y,x_v,y_v): # this is a box with tangible walls
    for i in range(N-quarantine):  
        x[i] += dt*x_v[i] 
        y[i] += dt*y_v[i]
        # for a finite box with tangible edges(bouncing)
        if (x[i]) >= size-.2 or x[i]<=.2:
            x_v[i] = -x_v[i]
            x[i] += x_v[i]*dt
            
        if (y[i]) >= size-.2 or y[i]<=.2:
            y_v[i] = -y_v[i]
            y[i] += y_v[i]*dt
            
    return x,y,x_v,y_v

infected_total = []
infected_total.append(np.sum(status))
uninfected = []

def infect(status,days,infected_total,uninfected):
    current = np.count_nonzero(status==1)
    for i in range(N-quarantine):
        if status[i]==1:
            days[i]+=1/(pause*10)
        j=N-1-quarantine
        while j>i: # this finds all the distances from all the points
            if (status[i]+status[j]==1 and np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)<=distance ):
                if status[i]==1 and masked[i]+masked[j]<np.random.random()*.5 and status[j]<1:
                    status[j]+=1
                if status[j]==1 and masked[i]+masked[j]<np.random.random()*.5 and status[i]<1:
                    status[i]+=1
            j-=1
        if days[i]>=rec_time and status[i]==1: # immunity 
            status[i]+=1
        if np.random.random()>.99998 and status[i]>1: # small chance of remission from 'recovered' 
            status[i]-=1 # number^ is small since this specific if statement runs thousands of times
            days[i]-=rec_time
        if current!=np.count_nonzero(status==1):
            infected_total.append(np.count_nonzero(status==1))
            uninfected.append(N-infected_total[-1])
        i+=1
    
    return status,days,infected_total,uninfected
         
def animate(i):
    for j in range(pause):
        xx,yy,xv,yv=stepper(x,y,x_v,y_v) # You can change this to stepper_tangible for bouncing!
        state,day,infected_num,uninfected_num = infect(status,days,infected_total,uninfected)
    plt.subplot(1,2,1)
    plt.cla()
    plt.axis([0.1,size-.1,0.1,size-.1])
    plt.scatter(xx,yy,c=colormap[state])
    # live graph
    plt.subplot(1,2,2)
    plt.cla()
    plt.title('Number of infected')
    plt.plot(infected_num,'r')
    plt.plot(uninfected_num,'g')
    j+=1
        
anim = ani.FuncAnimation(plt.gcf(),animate,interval=20)
plt.show()
#%%
# This is to find number of distances that are calculated per frame
# slightly less than n^2
def distances(particles):
    dists = [0,1]
    for i in range(particles):
        dists.append(dists[-1]+(dists[-1]-dists[-2])+1)
        i+=1
    print('Number of distances to calculate for ',particles,' particles: ',dists[-3])
    return dists
efficient = distances(N)
q = np.arange(0,N+2,1)
f = lambda a: a**2
inefficient = f(q)
plt.plot(q,inefficient,'r',label='N^2')
plt.plot(q,efficient,'b',label='My function')
plt.legend()
plt.title('Number of distance calculations required')

#%%
