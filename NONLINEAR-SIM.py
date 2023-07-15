# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 14:52:38 2023

@author: leona

NON-LINEAR QUADCOPTER SIMULATION

ref: https://www.hindawi.com/journals/je/2022/2449901/

TODO:
    implement missing resolution of all states:
        - inertial vs body frame transforms
        - ode compute on linear and rotationl equations?
        - figure out state feedback

"""
import numpy as np
from numpy import sin as s
from numpy import cos as c
from numpy import tan
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from time import time_ns

'''mechanical constants'''
g = 9.81
m = 0.250
l = 0.065                  # distance between rotor and center of mass
Ixx = 0.0001295
Iyy = 0.0001295
Izz = 0.0002494 
Jr = 0.5*0.0014*0.03175**2 # rotor inertia approximated as a thin disk
kt = 0.1863*10**(-6)       # thrust coefficient
kdx = 0.2                  # drag coefficients chosen arbitrarily similar to
kdy = 0.2                  #    the values in the paper with the drag model 
kdz = 0.2                  #    simplified to Fd = -kd*v

'''SET-POINTS'''
z_des   = 10
the_des = 0
phi_des = 0
psi_des = 0

# Inputs
w1 = 2100
w2 = 2100
w3 = 2100
w4 = 2100

# Propeller speed to thrust and torque
T    =   kt*(w1**2 + w2**2 + w3**2 + w4**2)
Tphi = l*kt*(-w2^2 + w4^2)
Tthe = l*kt*(-w1^2 + w3^2)
Tpsi = l*kt*(-w1^2 + w2^2 - w3^2 + w4^2)
wr   = -w1 + w2 - w3 + w4

# INITIAL CONDITIONS FOR ALL STATES
dp  = dq  = dr  = 0  # angular accelerations roll, pitch, yaw 
p   = q   = r   = 0  # rates of rotation     roll, pitch, yaw 
phi = the = psi = 0  # angle of pose         roll, pitch, yaw 
ddx = ddy = ddz = 0  # linear accelerations  forward, leftward, upward
dx  = dy  = dz  = 0  # linear speeds         forward, leftward, upward
x   = y   = z   = 0  # global position       forward, leftward, upward

def thrust_to_speed(Tphi,Tthe,Tpsi,T):
    #TODO: calculate speeds from desired thrusts and torques
    # ...
    return([w1,w2,w3,w4])

def controller(x):
    #TODO: given state, and PID params calculate desired thrust and torque
    # w1 =
    # w2 = 
    # w3 = 
    # w4 = 
    print('TODO')
    return(thrust_to_speed(Tphi,Tthe,Tpsi,T))

def ode(X,t):    
    # STATES
    p = X[0]
    q = X[1]
    r = X[2]
    
    dx = X[3]
    dy = X[4]
    dz = X[5]
    
    phi = X[6]
    the = X[7]
    psi = X[8]
    
    
    
    # ROTATIONAL DYNAMICS
    dpdt = q*r*(Iyy-Izz)/Ixx - Jr*( q/Ixx)*wr + Tphi/Ixx
    dqdt = p*r*(Izz-Ixx)/Iyy - Jr*(-p/Iyy)*wr + Tthe/Iyy
    drdt = p*q*(Ixx-Iyy)/Izz - Jr*( 0    )*wr + Tpsi/Izz
    # LINEAR DYNAMICS    
    ddxddt = -g*0 + T/m * (c(psi)*s(the)*c(phi) + s(psi)*s(phi)) - (1/m)*kdx*dx
    ddyddt = -g*0 + T/m * (s(psi)*s(the)*c(phi) + c(psi)*s(phi)) - (1/m)*kdy*dy
    ddzddt = -g*1 + T/m * (c(the)*c(phi))                        - (1/m)*kdz*dz
    # ANGLE TRANSFORMS
    dphidt = p + s(phi)*tan(the)*q + c(phi)*tan(the)*r
    dthedt = 0 + c(phi)*q          - s(phi)*r
    dpsidt = 0 + (s(phi)/c(the))*q + (c(phi)/c(the))*r
    # RETURN VALUES
    return([dpdt,dqdt,drdt,ddxddt,ddyddt,ddzddt,dphidt,dthedt,dpsidt])


X0 = [p,q,r,dx,dy,dz,phi,the,psi,dp,dq,dr,ddx,ddy,ddz,x,y,z]
# t_sim = np.linspace(0,5,1000)
# x = odeint(ode,X0,t_sim)
# print(x)

'''SIMULATION PARAMETERS'''
t_end = 5               # (s) total simulation time 
t_dly = 10 *(10**-3)    # (s) time-delay: how stale is the sensed state? 
t_cpu = 100*(10**-3)    # (s) compute-time: speed of the control loop

t_sim = 0               # (s) tracks simulation time. always starts at zero.
while(t_sim<=t_end):    # TODO
    # TO BE TRACKED: X, T_calc, t_clk
    # 1. get state snapshot, X_snp at time t = t_clk
    if(t_sim == 0):
        X = X0 # if first loop: import initial conditions
    X_sens = X
    # 2. using the real state X, simulate ode for time t = t_dly
    X = odeint(ode,X,np.linspace(t_clk, t_sim+t_dly, 100)
    # 3. calculate controller outputs T_des using X_sens
    W = controller(X_sens)
    
    # 4. using X and T_calc simulate ode for time t = t_cpu-t_dly


    t_sim += t_cpu - t_dly









