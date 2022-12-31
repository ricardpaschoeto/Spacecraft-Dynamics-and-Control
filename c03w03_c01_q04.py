from sympy import *
import numpy as np
from Controllers import *
from Attitude_Control import *

time = symbols('time')

f = 0.05
s0 = 0.2*sin(f*time)
s1 = 0.3*cos(f*time)
s2 = -0.3*sin(f*time)

sigma_rn = np.array([0, 0, 0])
sigma_bn = np.array([0.1, 0.2, -0.1])
omega_bn = np.deg2rad([30, 10,-20]) # rad/sec
L = np.array([0, 0, 0])
#L1 = np.array([0.5, -0.3, 0.2])

controller = Full_Controller()
sim_time = 120
ctrl = Attitude_Control(controller, 100, 75, 80, sigma_rn, sigma_bn,omega_bn, 5.0,10.0,L)
ctrl.simulator(sim_time, 0.01, 30)

# Norm sigma_br at  30 s:  0.1947365934979193