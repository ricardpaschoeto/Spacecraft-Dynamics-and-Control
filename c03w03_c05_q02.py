from sympy import *
import numpy as np
from Controllers import *
from Attitude_Control import *

sigma_rn = np.array([0, 0, 0])
sigma_bn = np.array([0.1, 0.2, -0.1])
omega_bn = np.deg2rad([30, 10,-20]) # rad/sec

L = np.array([0, 0, 0])
delta_L = np.array([0, 0, 0])

K = 5
P = np.array([22.3607,19.3649,20.0])

controller = Gain_Controller()
sim_time = 120
ctrl = Attitude_Control(controller, 100, 75, 80, sigma_rn, sigma_bn,omega_bn, K,P,L, delta_L)
ctrl.simulator(sim_time, 0.01, 30)

# Norm sigma_br at  30 s:  0.1328