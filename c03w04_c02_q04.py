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
L = np.array([0, 0, 0])
delta_L = np.array([0, 0, 0])

K = 5
P = 10

controller = Saturated_Controller(1.0)
sim_time = 180
target = 60
dt = 0.01

ctrl = Attitude_Control(controller, 100, 75, 80, sigma_rn, sigma_bn,omega_bn, K,P,L, delta_L)
ctrl.simulator(sim_time, dt, target)

# Norm sigma_br at  60 s:  0.5368617551817896