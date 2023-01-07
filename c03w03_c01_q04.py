from sympy import *
import numpy as np
from Controllers import *
from Attitude_Control import *

sigma_rn = np.array([0, 0, 0])
sigma_bn = np.array([0.1, 0.2, -0.1])
omega_bn = np.deg2rad([30, 10,-20]) # rad/sec
L = np.array([0, 0, 0])
delta_L = np.array([0, 0, 0])

controller = Full_Controller()
sim_time = 120
target  = 30
dt = 0.01

tvec = np.linspace(0, sim_time, int(sim_time/dt))

ctrl = Attitude_Control(controller, 100, 75, 80, sigma_rn, sigma_bn,omega_bn, 5.0,10.0,L, delta_L)
ctrl.simulator(sim_time, dt, target)

u_hist = np.array(controller.u_history)
print("Control State", u_hist[-1])
plt.plot(tvec, u_hist[:, 0], 'g')
plt.plot(tvec, u_hist[:, 1], 'b')
plt.plot(tvec, u_hist[:, 2], 'r')
plt.legend(["u0", "u1", "u2"], loc='upper right')
plt.grid()
plt.show()

# Norm sigma_br at  30 s:  0.1947365934979193