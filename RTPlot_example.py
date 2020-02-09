# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 17:07:46 2020

@author: Benedict
"""
import matplotlib.pyplot as plt
from core import molecule
import matplotlib.animation as animation
from matplotlib import style
import random as rd
import potential_funcs as pf

style.use('fivethirtyeight')
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

ScaleFactor = 10e12
DELTAX  =   10e-13*10e12 
m0      =   1.6e-27*10e12
q0      =   1.6e-19*10e11
N       =   100
molecules=[None]*N
for i in range(N):
    if i%2==0:
        molecules[i]=molecule([0,0],[DELTAX*4*rd.random(),DELTAX*4*rd.random()],m0,-q0,i,[pf.grad_potential_coulomb])
    else:
        molecules[i]=molecule([0,0],[DELTAX*4*rd.random(),DELTAX*4*rd.random()],m0,-q0,i,[pf.grad_potential_coulomb])
    molecules[i].calc_grid_pos()
posx=[None]*N
posy=[None]*N

def animate(j):
    for i in molecules:
        i.calc_grid_pos()
    for i in molecules:
        i.calc_new_pos(molecules)
    for index,i in enumerate(molecules):
        posx[index]=i.pos[0]
        posy[index]=i.pos[1]
    ax1.clear()
    ax1.scatter(posx,posy)
    print(j)

ani = animation.FuncAnimation(fig, animate)
plt.show()