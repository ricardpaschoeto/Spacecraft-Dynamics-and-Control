import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from scipy.spatial.transform import Rotation as R

time = symbols('time')
class Attitude_Control:
    def __init__(self, controller, I1, I2, I3, sigma_rn, sigma_bn, omega_bn, K, P, L, delta_L):
        
        # Controller type
        self.controller = controller

        # Principal Inertias
        self.I = np.diag([I1, I2, I3])

        # Initial states
        self.sigma_rn = sigma_rn
        self.omega_rn = np.array([0,0,0])
        self.sigma_bn = sigma_bn
        self.omega_bn = omega_bn
        
        self.omega_br = self.omega_bn - self.omega_rn
        self.sigma_br = np.array([0,0,0])

        if self.sigma_rn.dtype == 'O':
            self.s_rn_dot = self.sigma_rn_dot()
        else:
            self.s_rn_dot = np.array([0,0,0])
        
        # Gain Matrices and Torque
        self.K = K
        if isinstance(P, np.ndarray):
            self.P = np.diag(P)
        else:
            self.P = P * np.eye(3)
            
        self.L = L
        self.delta_L = delta_L

        # Errors List
        self.omega_br_list = []
        self.sigma_br_list = []

        # Reference List
        self.omega_rn_list = []
        self.sigma_rn_list = []

        # Body trajectory List
        self.omega_bn_list = []
        self.sigma_bn_list = []        

    def skew_matrix(self, value):
        value = np.squeeze(value)
        m = np.array([[        0, -value[2],  value[1]],
                      [ value[2],         0, -value[0]],
                      [-value[1],  value[0],         0]])

        return m

    def sigma_rn_dot(self):
        
        if self.sigma_rn.dtype != 'O':
            return self.sigma_rn

        time = symbols('time')
        s0_dot = diff(self.sigma_rn[0],time)
        s1_dot = diff(self.sigma_rn[1],time)
        s2_dot = diff(self.sigma_rn[2],time)

        sbn = np.array([s0_dot, s1_dot, s2_dot])
        
        return sbn

    def sigma_bn_dot(self):
        
        sdot = 0.25 * np.dot(((1 - np.dot(self.sigma_bn, self.sigma_bn)) * np.eye(3) + 2 * self.skew_matrix(self.sigma_bn) + 2 * np.outer(self.sigma_bn,  self.sigma_bn)), self.omega_bn)

        return sdot

    def sigma_shadow(self, sigma):
        norm = np.linalg.norm(sigma) ** 2
        return np.array([-i / norm for i in sigma])

    def B(self, sigma_rn):

        ss = np.dot(sigma_rn, sigma_rn)
        matrix_B = np.zeros((3,3))
        # First line
        matrix_B[0][0] = 1 - ss + 2 * sigma_rn[0]**2
        matrix_B[0][1] = 2 * (sigma_rn[0]*sigma_rn[1] - sigma_rn[2])
        matrix_B[0][2] = 2 * (sigma_rn[0]*sigma_rn[2] + sigma_rn[1])
        # Second line
        matrix_B[1][0] = 2 * (sigma_rn[1]*sigma_rn[0] + sigma_rn[2])
        matrix_B[1][1] = 1 - ss + 2 * sigma_rn[1]**2
        matrix_B[1][2] = 2 * (sigma_rn[1]*sigma_rn[2] - sigma_rn[0])
        # Third line
        matrix_B[2][0] = 2 * (sigma_rn[2]*sigma_rn[0] - sigma_rn[1])
        matrix_B[2][1] = 2 * (sigma_rn[2]*sigma_rn[1] + sigma_rn[0])
        matrix_B[2][2] = 1 - ss + 2 * sigma_rn[2]**2

        return matrix_B

    def omg_rn(self, t):

        if self.sigma_rn.dtype == 'O':
            sigma_rn = np.array(lambdify(time, list(self.sigma_rn))(t), dtype=float)            
        else:
            sigma_rn = self.sigma_rn

        if self.s_rn_dot.dtype != 'O':
            sigma_rn_dot = self.s_rn_dot
        else:
            sigma_rn_dot = np.array(lambdify(time, list(self.s_rn_dot))(t), dtype=float)

        B = self.B(sigma_rn)
        
        return 4 * np.dot(np.linalg.inv(B), sigma_rn_dot)

    def omg_rn_dot(self, t, dt):

        omega_rn_adv = self.omg_rn(t+dt)
        omega_rn_actual = self.omg_rn(t)

        return (omega_rn_adv - omega_rn_actual) / dt

    def omega_bn_dot(self, t, dt):

        if self.sigma_rn.dtype == 'O':
            sigma_rn = np.array(lambdify(time, list(self.sigma_rn))(t), dtype=float)            
        else:
            sigma_rn = self.sigma_rn
    
        self.sigma_br = -R.from_matrix(np.dot(R.from_mrp(self.sigma_bn).as_matrix().T, R.from_mrp(sigma_rn).as_matrix())).as_mrp()

        self.omega_rn = self.omg_rn(t)
        omega_rn_dot = self.omg_rn_dot(t,dt)

        BR = R.from_mrp(self.sigma_br).as_matrix().T

        omega_rn_bf = np.dot(BR, self.omega_rn)
        omega_rn_bf_0 = np.dot(BR, self.omg_rn(0))
        omega_rn_dot_bf = np.dot(BR, omega_rn_dot)
        omega_br = self.omega_bn - omega_rn_bf

        u = self.controller.control(self.K, self.P, self.L, self.I, self.sigma_br, self.sigma_bn, omega_br, self.omega_bn, omega_rn_bf_0, self.omega_rn, omega_rn_bf, omega_rn_dot_bf, t, dt)

        omega_bn_dot =  np.dot(np.linalg.inv(self.I) ,(-np.cross(self.omega_bn, np.dot(self.I , self.omega_bn)) + u + self.L + self.delta_L))

        return omega_bn_dot

    def simulator(self, sim_time, dt, t_target):

        tvec = np.linspace(0, sim_time, int(sim_time/dt + 1))
        t_target_vec = int(t_target/dt + 1)

        if self.sigma_rn.dtype == 'O':
            sigma_rn = np.array(lambdify(time, list(self.sigma_rn))(0), dtype=float)
        else:
            sigma_rn = self.sigma_rn

        self.omega_rn = self.omg_rn(0)
        self.sigma_br =  -R.from_matrix(np.dot(R.from_mrp(self.sigma_bn).as_matrix().T, R.from_mrp(sigma_rn).as_matrix())).as_mrp()
        
        self.sigma_rn_list.append(sigma_rn)
        self.sigma_br_list.append(self.sigma_br)
        self.sigma_bn_list.append(self.sigma_bn)
        self.omega_bn_list.append(self.omega_bn)
        self.omega_rn_list.append(self.omega_rn)     
        
        for t in tvec[1:]:

            if self.sigma_rn.dtype == 'O':
                sigma_rn = np.array(lambdify(time, list(self.sigma_rn))(t), dtype=float)            
            else:
                sigma_rn = self.sigma_rn
            
            self.sigma_rn_list.append(sigma_rn)
                
            self.sigma_br = -R.from_matrix(np.dot(R.from_mrp(self.sigma_bn).as_matrix().T, R.from_mrp(sigma_rn).as_matrix())).as_mrp()
            self.sigma_br_list.append(self.sigma_br)

            self.omega_rn = self.omg_rn(t)
            self.omega_rn_list.append(self.omega_rn)

            self.sigma_bn = self.sigma_bn + self.sigma_bn_dot() * dt
            self.omega_bn = self.omega_bn + self.omega_bn_dot(t,dt) * dt

            if np.dot(self.sigma_bn, self.sigma_bn) > 1:
                self.sigma_bn = self.sigma_shadow(self.sigma_bn)

            self.sigma_bn_list.append(self.sigma_bn)
            self.omega_bn_list.append(self.omega_bn)

            if t % 25 == 0:
                print("Simulated {} seconds".format(t))
        
        # Plot Results
        self.sigma_bn_list = np.array(self.sigma_bn_list)
        self.sigma_rn_list = np.array(self.sigma_rn_list)
        self.sigma_br_list = np.array(self.sigma_br_list)
        self.omega_bn_list = np.array(self.omega_bn_list)
        self.omega_rn_list = np.array(self.omega_rn_list)
        
        error_norm = [np.sqrt(i**2 + j**2 + k**2) for i, j, k in self.sigma_br_list]
        print("Norm sigma_br at ", t_target, "s: ", error_norm[t_target_vec])
        print("Steady State error sigma_ss at ", tvec[-1], "s: ", np.round(self.sigma_br_list[-1], 2))

        _, axs = plt.subplots(2,2)
        axs[0,0].plot(tvec, self.sigma_bn_list[:, 0], 'g')
        axs[0,0].plot(tvec, self.sigma_rn_list[:, 0], 'g--')
        axs[0,0].plot(tvec, self.sigma_bn_list[:, 1], 'b')
        axs[0,0].plot(tvec, self.sigma_rn_list[:, 1], 'b--')
        axs[0,0].plot(tvec, self.sigma_bn_list[:, 2], 'r')
        axs[0,0].plot(tvec, self.sigma_rn_list[:, 2], 'r--')
        axs[0,0].set_title('Attitude (sigma) history')
        axs[0,0].legend(["sigma_1", "sigma_1_target", "sigma_2", "sigma_2_target", "sigma_3", "sigma_3_target"], loc='upper right')
        axs[0,0].grid()

        axs[0,1].plot(tvec, self.omega_bn_list[:, 0], 'b')
        axs[0,1].plot(tvec, self.omega_rn_list[:, 0], 'b--')
        axs[0,1].plot(tvec, self.omega_bn_list[:, 1], 'r')
        axs[0,1].plot(tvec, self.omega_rn_list[:, 1], 'r--')
        axs[0,1].plot(tvec, self.omega_bn_list[:, 2], 'g')
        axs[0,1].plot(tvec, self.omega_rn_list[:, 2], 'g--')
        axs[0,1].set_title("Rate (w) history")
        axs[0,1].legend(["w_1", "w_1_target", "w_2", "w_2_target", "w_3", "w_3_target"], loc='upper right')
        axs[0,1].grid()

        axs[1,0].plot(tvec, error_norm, 'b')
        axs[1,0].set_title("Error Norm")
        axs[1,0].grid()

        axs[1,1].plot(tvec, self.sigma_br_list[:,0], 'b')
        axs[1,1].plot(tvec, self.sigma_br_list[:,1], 'r')
        axs[1,1].plot(tvec, self.sigma_br_list[:,2], 'g')
        axs[1,1].legend(["sigma_1", "sigma_2", "sigma_3"])
        axs[1,1].set_title("Sigma BR")
        axs[1,1].grid()

        plt.show()

