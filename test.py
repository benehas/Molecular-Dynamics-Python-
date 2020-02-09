# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 17:29:05 2020

@author: Bened
"""

import matplotlib.pyplot as plt
from core import molecule
import matplotlib.animation as animation
from matplotlib import style
import random as rd
import potential_funcs as pf


ScaleFactor = 10e12
DELTAX  =   10e-13*10e12 
m0      =   1.6e-27*10e12
q0      =   1.6e-19*10e12

molecules=[None]*3
for i in range(3):
    molecules[i]=molecule([0,0],[0,0],m0,-q0,i,[pf.grad_potential_coulomb])
    molecules[i].calc_grid_pos()

for i in molecules:
    print(i.pos)
    print(i.get_molecules_in_range(molecules))

