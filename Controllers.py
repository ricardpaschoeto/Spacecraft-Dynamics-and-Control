import numpy as np

class Controllers:

    def control(self, K, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf):
        pass


class Full_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf):

        return (-K * sigma_br - np.dot(P, omega_br) + np.dot(I, omega_rn_dot_bf - np.cross(omega_bn, omega_rn_bf)) + np.cross(omega_bn, np.dot(I, omega_bn)) - L)

class Proportional_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf):

        return (-K * sigma_br - np.dot(P, omega_br))