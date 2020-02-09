# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 17:14:33 2020

@author: Bened
"""

ScaleFactor = 10e12
e0      =   8.8e-12*10e12
G       =   6.6e-11*10e12
epsilon =   10e-38
rm      =   10e-10*10e12


def grad_potential_grav(molecule1,molecule2,r):
    x=molecule1.pos[0]-molecule2.pos[0]
    y=molecule1.pos[1]-molecule2.pos[1]
    return (-G*molecule1.mass*molecule1.mass*2*x*pow(x*x+y*y,-1.5),-G*molecule1.mass*molecule1.mass*2*y*pow(x*x+y*y,-1.5))


def grad_potential_coulomb(molecule1,molecule2):
    x=(molecule1.pos[0]-molecule2.pos[0])
    y=(molecule1.pos[1]-molecule2.pos[1])
    return (1/(4*3.14*e0)*molecule1.charge*molecule1.charge*2*x*pow(x*x+y*y,-1.5),1/(4*3.14*e0)*molecule1.charge*molecule1.charge*2*y*pow(x*x+y*y,-1.5))
