from sympy import *
import numpy as np
from Controllers import *
from Attitude_Control import *

sigma_rn = np.array([0, 0, 0])
sigma_bn = np.array([-0.3, -0.4, 0.2])
omega_bn = np.array([0.2, 0.2,0.2]) # rad/sec

K = 1.0
P = 3.0

L = np.array([0, 0, 0])
delta_L = np.array([0., 0., 0.])

controller = CLD_Controller()
sim_time = 15
target  = 10
dt = 0.01

tvec = np.linspace(0, sim_time, int(sim_time/dt))

ctrl = Attitude_Control(controller, 30, 20, 10, sigma_rn, sigma_bn,omega_bn, K,P,L, delta_L)
ctrl.simulator(sim_time, dt, target)

u_hist = np.array(controller.u_history)
print("Control State", u_hist[-1])
plt.plot(tvec, u_hist[:, 0], 'g')
plt.plot(tvec, u_hist[:, 1], 'b')
plt.plot(tvec, u_hist[:, 2], 'r')
plt.legend(["u0", "u1", "u2"], loc='upper right')
plt.grid()
plt.show()