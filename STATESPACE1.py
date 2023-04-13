# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 22:23:34 2023

DRONE STATE SPACE SIMULATION AND CONTROL TESTING ENVIRONMENT

"""

import numpy as np
import sys
import scipy as sp
import control.matlab as con
from control.matlab import TransferFunction as tf
import sympy
import filterpy
from matplotlib import pyplot as plt
plt.close('all')

#TODO: check paper specs on drone orientation (x or +)
#TODO: find a set of good first estimates for the constants initalized below
#TODO: determine how to simulate a disturbance.
#      in theory, the equation becomes dx = Ax + Bu + Fd
#TODO: work towards implementation of place() or acker()

'INITIALIZE CONSTANTS'

g = 9.81
Ixx =0.1
Iyy = 0.1
Izz = 0.1
m = 0.25
Jxx = 0.1
Jyy = 0.1
Jzz = 0.1
b = 1.0
l = 0.2
d = 1.0

'''
DEFINE STATE SPACE

In the form 
dx = Ax + Bu; Y = Cx

where the state matrix
x = [x,y,z,phi,theta,psi,dx,dy,dz,d(phi),d(theta),d(psi)]

phi:   pitch
theta: roll
psi:   yaw

the input matrix
u = [thrust (N), roll (Nm), pitch (Nm), yaw (Nm)]

'''

A = np.array([
    [0,0,0,0,0,0,Ixx,0,0,0,0,0],
    [0,0,0,0,0,0,0,Iyy,0,0,0,0],
    [0,0,0,0,0,0,0,0,Izz,0,0,0],
    [0,0,0,0,0,0,0,0,0,Ixx,0,0],
    [0,0,0,0,0,0,0,0,0,0,Iyy,0],
    [0,0,0,0,0,0,0,0,0,0,0,Izz],
    [0,0,0,0,0,g,  0,0,0,0,0,0],
    [0,0,0,0,-g, 0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0  ],
    [0,0,0,0,0,0,0,0,0,0,0,0  ],
    [0,0,0,0,0,0,0,0,0,0,0,0  ],
    [0,0,0,0,0,0,0,0,0,0,0,0  ]
    ])

B = np.array([
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],    
    [1/m,  0,0,0],
    [0,1/Jxx,0,0],
    [0,0,1/Jyy,0],
    [0,0,0,1/Jzz]    
    ])

C = np.array([
    [Ixx,0,0,0,0,0,0,0,0,0,0,0],
    [0,Iyy,0,0,0,0,0,0,0,0,0,0],
    [0,0,Izz,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,0,0,0,0,0  ],
    ])

D = np.zeros((4,4))

sys = con.ss(A,B,C,D) # OLTF

# Eigenvalues of A = poles of the system
eig_A = np.linalg.eigvals(A) # Consists of 12 zeros. 
# In ECS we would now manually place 12 poles according to some criteria, such
# as response time, setting time, and steady state error.

'SIMULATE'

t_sim = np.linspace(0,5,1000)
u_sim = np.tile(np.array([[0.1,0,0,0]]),(1000,1)) # constant input for a first test

yout, T, xout = con.lsim(sys,U=u_sim,T=t_sim)

'PLOT RESULTS'

fig = plt.figure(1)

plt.subplot(3,4,1)
plt.plot(T,xout[:,0])
plt.grid()
plt.title('x')

plt.subplot(3,4,2)
plt.plot(T,xout[:,1])
plt.grid()
plt.title('y')

plt.subplot(3,4,3)
plt.plot(T,xout[:,2])
plt.grid()
plt.title('z')

plt.subplot(3,4,4)
plt.plot(T,xout[:,3])
plt.grid()
plt.title('phi')

plt.subplot(3,4,5)
plt.plot(T,xout[:,4])
plt.grid()
plt.title('theta')

plt.subplot(3,4,6)
plt.plot(T,xout[:,5])
plt.grid()
plt.title('psi')

plt.subplot(3,4,7)
plt.plot(T,xout[:,6])
plt.grid()
plt.title('dx')

plt.subplot(3,4,8)
plt.plot(T,xout[:,7])
plt.grid()
plt.title('dy')

plt.subplot(3,4,9)
plt.plot(T,xout[:,8])
plt.grid()
plt.title('dz')

plt.subplot(3,4,10)
plt.plot(T,xout[:,9])
plt.grid()
plt.title('d(phi)')

plt.subplot(3,4,11)
plt.plot(T,xout[:,10])
plt.grid()
plt.title('d(theta)')

plt.subplot(3,4,12)
plt.plot(T,xout[:,11])
plt.grid()
plt.title('d(psi)')

plt.tight_layout()
