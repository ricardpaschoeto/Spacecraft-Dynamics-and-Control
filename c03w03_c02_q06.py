from sympy import *
import numpy as np
from Controllers import *
from Attitude_Control import *

time = symbols('time')

f = 0.05
s0 = 0.2*sin(f*time)
s1 = 0.3*cos(f*time)
s2 = -0.3*sin(f*time)

sigma_rn = np.array([s0, s1, s2])
sigma_bn = np.array([0.1, 0.2, -0.1])
omega_bn = np.deg2rad([30, 10,-20]) # rad/sec
L = np.array([0.5, -0.3, 0.2])
delta_L = np.array([0, 0, 0])

controller = Torqueless_Controller()
sim_time = 120
ctrl = Attitude_Control(controller, 100, 75, 80, sigma_rn, sigma_bn,omega_bn, 5.0,10.0,L, delta_L)
ctrl.simulator(sim_time, 0.01, 80)

# Norm sigma_br at  80 s:  0.01701046199647912