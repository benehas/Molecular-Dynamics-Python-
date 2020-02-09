# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 17:22:01 2020

@author: Bened
"""
import matplotlib.pyplot as plt
from core import molecule
import random as rd
import potential_funcs as pf

ScaleFactor = 10e12
DELTAX  =   10e-13*10e12 
m0      =   1.6e-27*10e12
q0      =   1.6e-19*10e12

molecules=[None]*100
for i in range(100):
    if i%2==0:
        molecules[i]=molecule([0,0],[DELTAX*4*rd.random(),DELTAX*4*rd.random()],m0,-q0,i,[pf.grad_potential_coulomb])
    else:
        molecules[i]=molecule([0,0],[DELTAX*4*rd.random(),DELTAX*4*rd.random()],m0,+q0,i,[pf.grad_potential_coulomb])
    molecules[i].calc_grid_pos()

posxl=[list()]*40
posyl=[list()]*40


for j in range(300):
    for i in molecules:
        i.calc_grid_pos()
    for i in molecules:
        i.calc_new_pos(molecules)
    if j%4==0:
        for k in range(40):
            posxl[k].append(molecules[k].pos[0])
            posyl[k].append(molecules[k].pos[1])
for k in range(40):
    plt.scatter(posxl[k],posyl[k])
