# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 22:23:34 2023

DRONE STATE SPACE SIMULATION AND CONTROL TESTING ENVIRONMENT

LINEARIZED

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

A: contains the physics (state-propagation)
B: contains the input relation
Y: the outputs (what is measured)
C: defines (and scales) which state is measured
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
eig_A = np.linalg.eigvals(A) 
print('Poles:',eig_A)# Consists of 12 zeros. 
# In ECS we would now manually place 12 poles according to some criteria, such
# as response time, setting time, and steady state error.

#########################

# CALCULATE ACTIVE POLES VIA DESIGN CRITERIA

zeta = 0.9
#NATURAL FREQUENCY (wn) GIVEN DAMPING (zeta) AND RISE TIME (tr)
tr = 0.1
wn = float((np.pi - np.arctan(np.sqrt(1-zeta**2)/zeta))/np.sqrt(1-zeta**2)/tr) 
# DAMPED FREQUENY (wd)
wd = wn*np.sqrt(1-np.sqrt(zeta**2))
pc1 = complex(-zeta*wn,wd)
pc2 = complex(-zeta*wn,-wd)

# # Choose 10 additional poles (this is fucking stupid)
q=4#10.1
r=5#10.2
s=6#10.3
P03 = q *-zeta*wn
P04 = q *-zeta*wn
P05 = q *-zeta*wn
P06 = q *-zeta*wn
P07 = r *-zeta*wn
P08 = r *-zeta*wn
P09 = r *-zeta*wn
P10 = r *-zeta*wn
P11 = s *-zeta*wn
P12 = s *-zeta*wn
J = [ pc1,pc2,P03,P04,P05,P06,P07,P08,P09,P10,P11,P12 ]

print()
print("Design specifications")
print("---------------------------------------------")
print("Damping = ",round(zeta,2))
print("Rise time = ", round(tr,2),"s")
print("Natural frequency = ", round(wn/(2*np.pi),2),"Hz")
print("Damped frequency = ", round(wd/(2*np.pi),2),"Hz")
print("-zeta*wm = ",-zeta*wn)
print("Design poles = ",J)

# 3. Determine the state feedback controller K 

#K = con.acker(A,B,J)
K = con.place(A,B,J)

print()
print("Gain matrix")
print("---------------------------------------------")
#print("Acker result K=",con.acker(A,B,J))
print("Place result K=",con.place(A,B,J))
print()


# 4. Calculate the feedforward matrix

invN = C*np.linalg.inv(A-B*K)*B*-1
N= np.linalg.pinv(invN)
print()
print("Feedforward matrix")
print("---------------------------------------------")
print('N=', N )
print()


# 5. Close the loop

A_cl = A - B*K
print() # Check the closed loop matrix to verify pole placement
print("Closed loop poles")
print("---------------------------------------------")
print(np.linalg.eigvals(A_cl))
print()

sys_CL = con.ss(A_cl,B*N,C,D*N)

#########################

'SIMULATE'

t_sim = np.linspace(0,1,100000)
u_sim = np.tile(np.array([[0.1,0,0,0]]),(100000,1)) 
# constant input for a first test

yout, T, xout = con.lsim(sys_CL,U=u_sim,T=t_sim)

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
