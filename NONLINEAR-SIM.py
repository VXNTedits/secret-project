# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 14:52:38 2023

@author: leona

NON-LINEAR QUADCOPTER SIMULATION

ref: https://www.hindawi.com/journals/je/2022/2449901/

TODO:
    * implement missing resolution of all states:
        - for some states you might need to compute the remaining
          inertial vs body frame transforms
          
        - ode figure out how to solve the 2nd ord ode to obtain positions!!!
        
    * figure out control scheme

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

# INITIAL CONDITIONS FOR ALL STATES AND INPUTS
# States
dp  = dq  = dr  = 0  # angular accelerations roll, pitch, yaw 
p   = q   = r   = 0  # rates of rotation     roll, pitch, yaw 
phi = the = psi = 0  # angle of pose         roll, pitch, yaw 
ddx = ddy = ddz = 0  # linear accelerations  forward, leftward, upward
dx  = dy  = dz  = 0  # linear speeds         forward, leftward, upward
x   = y   = z   = 0  # global position       forward, leftward, upward
X0 = [p,q,r,dx,dy,dz,phi,the,psi,dp,dq,dr,ddx,ddy,ddz,x,y,z]
# Inputs
w1 = 0
w2 = 0
w3 = 0
w4 = 0
W0 = [w1,w2,w3,w4]
# Propeller speed to thrust and torques
Tthr =   kt*(w1**2 + w2**2 + w3**2 + w4**2)
Tphi = l*kt*(-w2^2 + w4^2)
Tthe = l*kt*(-w1^2 + w3^2)
Tpsi = l*kt*(-w1^2 + w2^2 - w3^2 + w4^2)
T0 = [Tphi,Tthe,Tpsi,Tthr]

def controller(X):
    #TODO: given state and PID gains, calculate desired thrust and torque
    
    # CONTROLLER GAINS
    # ...
    # ...
    
    # ACTUATOR LIMITS
    w_max = 24200*(2*np.pi/60) # (rad/s)    if Imax = 3.5A
    dw_max = 4200              # (rad/s^2)  this is a severe underestimation
                               #            based on 3Nmm of torque, but safe 
                               #            for a first design attempt

    ############### PLACE HOLDER - REPLACE WITH CONTROL LAWS ###############    
    w1 = 1905
    w2 = 1905
    w3 = 1905
    w4 = 1905
    # speed to thrust and torque transformations
    Tthr =   kt*(w1**2 + w2**2 + w3**2 + w4**2)
    Tphi = l*kt*(-w2^2 + w4^2)
    Tthe = l*kt*(-w1^2 + w3^2)
    Tpsi = l*kt*(-w1^2 + w2^2 - w3^2 + w4^2)
    ############### PLACE HOLDER - REPLACE WITH CONTROL LAWS ############### 
    
    return([w1,w2,w3,w4],[Tphi,Tthe,Tpsi,Tthr])

def ode(X,t):
    # TODO
    # possible option: include W as a state with simple constraints, take
    # dwdt = cst.
    # in this case, the controller outputs no-longer match the actuator effort
    # so you'd have to rewrite to make the distinction between 
    # Tdes, Wdes and T,W
    #
    
    # UNPACK RELEVANT STATES
    p = X[0]
    q = X[1]
    r = X[2]
    dx = X[3]
    dy = X[4]
    dz = X[5]    
    phi = X[6]
    the = X[7]
    psi = X[8]  
    
    # UNPACK INPUTS
    Tphi = T[0]
    Tthe = T[1]
    Tpsi = T[2]
    Tthr = T[3]
    w1 = W[0]
    w2 = W[1]
    w3 = W[2]
    w4 = W[3]
    
    # ROTATIONAL DYNAMICS
    dp = q*r*(Iyy-Izz)/Ixx - Jr*( q/Ixx)*(-w1+w2-w3+w4) + Tphi/Ixx
    dq = p*r*(Izz-Ixx)/Iyy - Jr*(-p/Iyy)*(-w1+w2-w3+w4) + Tthe/Iyy
    dr = p*q*(Ixx-Iyy)/Izz - Jr*( 0    )*(-w1+w2-w3+w4) + Tpsi/Izz
    
    # LINEAR DYNAMICS    
    ddx = -g*0 + Tthr/m * (c(psi)*s(the)*c(phi) + s(psi)*s(phi)) - (1/m)*kdx*dx
    ddy = -g*0 + Tthr/m * (s(psi)*s(the)*c(phi) + c(psi)*s(phi)) - (1/m)*kdy*dy
    ddz = -g*1 + Tthr/m * (c(the)*c(phi))                        - (1/m)*kdz*dz
    
    # ANGLE TRANSFORMS
    dphi = p + s(phi)*tan(the)*q + c(phi)*tan(the)*r
    dthe = 0 + c(phi)*q          - s(phi)*r
    dpsi = 0 + (s(phi)/c(the))*q + (c(phi)/c(the))*r    
    
    return([dp,dq,dr,ddx,ddy,ddz,dphi,dthe,dpsi])


'''SIMULATION PARAMETERS'''
t_end = 5               # (s) total simulation time 
t_dly = 10 *(10**-3)    # (s) time-delay: how stale is the sensed state? 
t_cpu = 100*(10**-3)    # (s) compute-time: speed of the control loop
resol = 1  *(10**-3)    # (s) simulation resolution

# initialize arrays to record simulation data
X_log = np.zeros([1,18])
W_log = np.zeros([1,4])
T_log = np.zeros([1,4])
t_log = np.zeros([1])

t_sim = 0
while(t_sim<=t_end):    # TODO
    # The state is stored in this format:    
    # X = [p,q,r,dx,dy,dz,phi,the,psi,dp,dq,dr,ddx,ddy,ddz,x,y,z].
    
    # 0. If first loop: import initial conditions.
    if(t_sim == 0):
        X = X0 
        W = W0
        T = T0
        
    # 1. Get state snapshot at time t = t_sim.
    X_sens = X 
    
    # 2. Using the real state X, simulate ode for time t = t_dly.
    #    X_sim should take the form: [p,q,r,dx,dy,dz,phi,the,psi]
    dt1 = np.linspace(t_sim, t_sim+t_dly, round(t_dly/resol))
    Xt1 = odeint(ode,X,dt1)
    t_log = np.append(dt1)
    
    # 3. Using the sim results, compute the remaining derived states (?),  
    #    and append the whole thing to X_log, where every row Rx is a complete
    #    set of states X. Every Rx must correspond to the timestamp present in 
    #    the same index in t_log.
    X_log = np.append(Xt1) #TODO
    T_log = np.append(T)   #TODO
    W_log = np.append(W)   #TODO
    
    # 4. Update controller outputs using X_sens.  
    W, T = controller(X_sens)
    
    # 5. Using X and W,T simulate ode for time t = t_cpu-t_dly.
    dt2 = np.linspace(t_sim, t_sim+t_cpu-t_dly, round((t_cpu-t_dly)/resol))
    Xt2 = odeint(ode,X,dt2)
    t_log = np.append(dt2)
    
    # 6. Repeat step 3
    X_log = np.append(Xt2) #TODO
    T_log = np.append(T)   #TODO
    W_log = np.append(W)   #TODO
    
    t_sim += (t_cpu - t_dly)


# finally, plot results






