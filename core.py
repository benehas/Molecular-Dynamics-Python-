# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 17:41:23 2020

@author: Bened
"""

import math

ScaleFactor = 10e12
DELTAX  =   10e-13*10e12 
DELTAT  =   10e-14*10e13 # timespacing
e0      =   8.8e-12*10e12
m0      =   1.6e-27*10e12
q0      =   1.6e-19*10e12
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

class molecule:
    __slots__=['vel','pos','mass','charge','id','grid_pos','grad_potential_funcs']
    
    def __init__(self,vel,pos,mass,charge,id,grad_potential_funcs):
        self.vel=vel
        self.pos=pos
        self.mass=mass
        self.charge=charge
        self.id=id
        self.grad_potential_funcs=grad_potential_funcs

    def calc_grid_pos(self):
        self.grid_pos=[int(self.pos[0]/DELTAX),int(self.pos[1]/DELTAX)]

    def calc_force(self,molecules_in_range):
        Fx=0
        Fy=0
        for mol in molecules_in_range:
            for funcs in self.grad_potential_funcs:                
                A,B=funcs(self,mol)
                Fx=Fx+A
                Fy=Fy+B
            return (Fx,Fy)
    def get_molecules_in_range(self,molecules):
        return [x for x in molecules if abs(x.grid_pos[0]-self.grid_pos[0])<=1 and abs(x.grid_pos[1]-self.grid_pos[1])<=1 and self.id!=x.id]
            
    def calc_new_pos(self,molecules):
        F=self.calc_force(self.get_molecules_in_range(molecules))
        try:
            self.pos[0]=self.pos[0]+DELTAT*self.vel[0]+0.5*DELTAT*DELTAT*F[0]/self.mass
            self.pos[1]=self.pos[1]+DELTAT*self.vel[1]+0.5*DELTAT*DELTAT*F[1]/self.mass
            self.vel[0]=self.vel[0]+DELTAT*F[0]/self.mass
            self.vel[1]=self.vel[1]+DELTAT*F[1]/self.mass
        except:
            self.pos[0]=self.pos[0]+DELTAT*self.vel[0]
            self.pos[1]=self.pos[1]+DELTAT*self.vel[1]
