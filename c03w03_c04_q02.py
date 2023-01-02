from sympy import *
import numpy as np
from Controllers import *
from Attitude_Control import *
import matplotlib.pyplot as plt

time = symbols('time')

f = 0.05
s0 = 0.2*sin(f*time)
s1 = 0.3*cos(f*time)
s2 = -0.3*sin(f*time)

sigma_rn = np.array([s0, s1, s2])
sigma_bn = np.array([0.1, 0.2, -0.1])
omega_bn = np.deg2rad([3, 1,-2]) # rad/sec

L = np.array([0, 0, 0])
delta_L = np.array([0.5, -0.3, 0.2])

K = 5
P = 10
Ki = 0.005

sim_time = 240
target = 45
dt = 0.01
tvec = np.linspace(0, sim_time, int(sim_time/dt))

controller = Integral_Controller(omega_bn, Ki)
ctrl = Attitude_Control(controller, 100, 75, 80, sigma_rn, sigma_bn,omega_bn, K,P,L, delta_L)
ctrl.simulator(sim_time, dt, target)
z_hist = np.array(controller.z_history)
print("predicted steady state error", z_hist[-1])
plt.plot(tvec, z_hist[:, 0], 'g')
plt.plot(tvec, z_hist[:, 1], 'b')
plt.plot(tvec, z_hist[:, 2], 'r')
plt.legend(["z0", "z1", "z2"], loc='upper right')
plt.show()

# Norm sigma_br at  45 s:  0.01443
# predicted steady state error [ 9.99987896 -6.00029445  3.99819536]