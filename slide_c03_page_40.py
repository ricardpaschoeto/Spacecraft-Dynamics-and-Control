from sympy import *
import numpy as np
from Controllers import *
from Attitude_Control import *
import matplotlib.pyplot as plt

sigma_rn = np.array([0, 0, 0])
sigma_bn = np.array([-0.3, -0.4, 0.2])
omega_bn = np.array([0.2, 0.2,0.2]) # rad/sec

L = np.array([0, 0, 0])
delta_L = np.array([0.05, 0.10, -0.10])

K = 1
P = 3
Ki = 0.01

sim_time = 150
target = 35
dt = 0.01
tvec = np.linspace(0, sim_time, int(sim_time/dt))

controller = Integral_Controller(omega_bn, Ki)
ctrl = Attitude_Control(controller, 10, 10, 10, sigma_rn, sigma_bn,omega_bn, K,P,L, delta_L)
ctrl.simulator(sim_time, dt, target)

z_hist = np.array(controller.z_history)
print("predicted steady state error", np.round(z_hist[-1], 3))
plt.plot(tvec, z_hist[:, 0], 'g')
plt.plot(tvec, z_hist[:, 1], 'b')
plt.plot(tvec, z_hist[:, 2], 'r')
plt.legend(["z0", "z1", "z2"], loc='upper right')
plt.show()