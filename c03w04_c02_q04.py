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
umax = np.array([1., 1., 1.])

K = 5
P = 10

controller = Saturated_Controller(umax)
sim_time = 180
target = 60
dt = 0.01

tvec = np.linspace(0, sim_time, int(sim_time/dt))

ctrl = Attitude_Control(controller, 100, 75, 80, sigma_rn, sigma_bn,omega_bn, K,P,L, delta_L)
ctrl.simulator(sim_time, dt, target)

u_hist = np.array(controller.u_history)
print("Control State", u_hist[-1])
plt.plot(tvec, u_hist[:, 0], 'g')
plt.plot(tvec, u_hist[:, 1], 'b')
plt.plot(tvec, u_hist[:, 2], 'r')
plt.legend(["u0", "u1", "u2"], loc='upper right')
plt.show()

# Norm sigma_br at  60 s:  0.520774528734927