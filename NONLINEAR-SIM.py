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
from matplotlib import pyplot as plt
plt.close('all')

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

'''set-points'''
z_des   = -10              # z axis points downward
the_des = 0
phi_des = 0
psi_des = 0

'''initial conditions for all states and inputs'''
# States
dp  = dq  = dr  = 0  # angular accelerations roll, pitch, yaw 
p   = q   = r   = 0  # rates of rotation     roll, pitch, yaw 
phi = the = psi = 0  # angle of pose         roll, pitch, yaw 
ddx = ddy = ddz = 0  # linear accelerations  forward, leftward, upward
dx  = dy  = dz  = 0  # linear speeds         forward, leftward, upward
x   = y   = z   = 0  # global position       forward, leftward, upward
X0 = [p,q,r, x,y,z, dx,dy,dz, phi,the,psi, dp,dq,dr, ddx,ddy,ddz]
# Inputs
w1 = 0
w2 = 0
w3 = 0
w4 = 0
W0 = [w1,w2,w3,w4]
# Equivalent thrust
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
    
    x = X[3]
    y = X[4]
    z = X[5]
    
    dx = X[6]
    dy = X[7]
    dz = X[8]    
    
    phi = X[9]
    the = X[10]
    psi = X[11]  
    
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
    ddx = -g*0 + Tthr/m * (c(psi)*s(the)*c(phi) + s(psi)*s(phi)) - 1/m * kdx*dx
    ddy = -g*0 + Tthr/m * (s(psi)*s(the)*c(phi) + c(psi)*s(phi)) - 1/m * kdy*dy
    ddz = -g*1 + Tthr/m * (c(the)*c(phi))                        - 1/m * kdz*dz
    
    # s1 = dx # Reduced order method generates a set of coupled 1st-order ODE's. 
    # s2 = dy # The new set of ODE's looks like just assigning variables 
    # s3 = dz # unneccessarily, but explicitly adding ODE's is required.    
    # s1 = -g*0 + Tthr/m * (c(psi)*s(the)*c(phi) + s(psi)*s(phi)) - (1/m)*kdx*x
    # s2 = -g*0 + Tthr/m * (s(psi)*s(the)*c(phi) + c(psi)*s(phi)) - (1/m)*kdy*y
    # s3 = -g*1 + Tthr/m * (c(the)*c(phi))                        - (1/m)*kdz*z
    
    # ANGLE TRANSFORMS
    dphi = p + s(phi)*tan(the)*q + c(phi)*tan(the)*r
    dthe = 0 + c(phi)*q          - s(phi)*r
    dpsi = 0 + (s(phi)/c(the))*q + (c(phi)/c(the))*r    
    
    #return([dp,dq,dr, dx,dy,dz, ddx,ddy,ddz, dphi,dthe,dpsi])
    return([dp,dq,dr, dx,dy,dz, ddx,ddy,ddz, dphi,dthe,dpsi])


def deltas(X_log,Xt,dt):
    n_rows, n_cols = Xt.shape
    deltas = np.zeros((n_rows, 6))
    for row_index in range(n_rows):
        pt = Xt[row_index,0]
        qt = Xt[row_index,1]
        rt = Xt[row_index,2]
        dxt = Xt[row_index,6]
        dyt = Xt[row_index,7]
        dzt = Xt[row_index,8]
        try:
            # This should always be ok, since we have initial conditions
            p0 = X_log[-1,0]
            q0 = X_log[-1,1]
            r0 = X_log[-1,2]
            dx0 = X_log[-1,6]
            dy0 = X_log[-1,8]
            dz0 = X_log[-1,8]
        except IndexError:
            # This should never happen
            print("Error in deltas operation: Index out of bounds.")
            p0 = q0 = r0 = dx0 = dy0 = dz0 = 0
        finally:
            dpdt = (pt-p0)/dt
            dqdt = (qt-q0)/dt
            drdt = (rt-r0)/dt
            ddxddt = (dxt-dx0)/dt
            ddyddt = (dyt-dy0)/dt
            ddzddt = (dzt-dz0)/dt
        new_row = np.array([dpdt,dqdt,drdt, ddxddt,ddyddt,ddzddt])
        deltas[row_index, :] = new_row     
    return(deltas)
    

'''SIMULATION PARAMETERS'''
t_end = 5               # (s) total simulation time 
t_dly = 10 *(10**-3)    # (s) time-delay: how stale is the sensed state? 
t_cpu = 100*(10**-3)    # (s) compute-time: speed of the control loop
resol = 1  *(10**-3)    # (s) simulation resolution

# initialize arrays with initial conditions to record simulation data
X_log = np.zeros([1,18])
W_log = np.zeros([1,4])
T_log = np.zeros([1,4])
t_log = np.zeros([1])
W = W0
T = T0
t_sim = 0
while(t_sim<=t_end):
    # The state is stored in this format:   
    #        [0,1,2  3,4,5   6, 7, 8,   9, 10, 11, 12,13,14,  15, 16, 17]
    # X    = [p,q,r, x,y,z, dx,dy,dz, phi,the,psi, dp,dq,dr, ddx,ddy,ddz].
    # Xsim = [p,q,r, x,y,z, dx,dy,dz, phi,the,psi]
    #
    # 0. If first loop: import initial conditions.
    if(t_sim == 0):
        X_log[0,:] = X0
        W_log[0,:] = W0
        T_log[0,:] = T0
        
    # 1. SENSOR READINGS. Get state snapshot at time t = t_sim.
    X_sens = X_log[-1] 
    
    # 2. PROCESSING LATENCY. Using the real state X, simulate ode for time 
    #    t = t_dly. 
    #    X_sim should take the form: [p,q,r, x,y,z, dx,dy,dz, phi,the,psi]
    dt1 = np.linspace(t_sim, t_sim+t_dly, round(t_dly/resol))
    Xt1 = odeint(ode, X_log[-1][:12], dt1)

    # 3. Using the sim results, derive the remaining states and append the 
    #    whole thing to X_log, where every row R is a complete set of states X. 
    #    Every R must correspond to the timestamp present in the same index in 
    #    t_log.
    dels = deltas(X_log, Xt1, resol)
    Xt1 = np.concatenate((Xt1, dels), axis=1)
    X_log = np.concatenate((X_log, Xt1), axis=0)
    t_log = np.concatenate((t_log, dt1), axis=0)
    
    # create array Tt1 with number of rows = len(dt1)
    # fill each row with T_log[-1]
    Tt1 = np.tile(T_log[-1], (len(dt1), 1))
    # append Tt1 onto the end of T_log
    T_log = np.concatenate((T_log, Tt1), axis=0)
    # repeat for W
    Wt1 = np.tile(W_log[-1], (len(dt1), 1))
    W_log = np.concatenate((W_log, Wt1), axis=0)
    
    # 4. CONTROL ACTION. Update controller outputs using X_sens.  
    W, T = controller(X_sens)
    
    # 5. COMPUTATION FREQUENCY. Using X and W,T simulate ode for time 
    #    t = t_cpu-t_dly.
    dt2 = np.linspace(t_sim, t_sim+t_cpu-t_dly, round((t_cpu-t_dly)/resol))
    Xt2 = odeint(ode, X_log[-1][:12], dt2)
    
    # 6. Repeat step 3, this time for the second simulated delay, Xt2.
    dels = deltas(X_log, Xt2, resol)
    Xt2 = np.concatenate((Xt2, dels), axis=1)
    X_log = np.concatenate((X_log, Xt2), axis=0)
    t_log = np.concatenate((t_log, dt2), axis=0)
    # This time append the newly calculated T and W values!
    Tt1 = np.tile(T, (len(dt1), 1))  
    Wt1 = np.tile(W, (len(dt1), 1))
    W_log = np.concatenate((W_log, Wt1), axis=0)
    
    # 7. Track simulation time
    t_sim += (t_cpu - t_dly)


# finally, plot results 
fig, ax = plt.subplots()
ax.plot(t_log, X_log[:,5])
plt.show()



