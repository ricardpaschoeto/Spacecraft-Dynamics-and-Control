import numpy as np

class Controllers:

    def control(self):
        pass

class Full_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):

        return (-K * sigma_br - np.dot(P, omega_br) + np.dot(I, omega_rn_dot_bf - np.cross(omega_bn, omega_rn_bf)) + np.cross(omega_bn, np.dot(I, omega_bn)) - L)

class Proportional_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):

        return (-K * sigma_br - np.dot(P, omega_br))

class Integral_Controller(Controllers):
    def __init__(self, omega_bn_0, Ki):

        self.omega_bn_0 = omega_bn_0
        self.omega_rn_0 = np.array([0,0,0])
        self.Ki = Ki
        self.s = np.array([0,0,0])
        self.z_history = []

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t, dt):

        if t == 0:
            self.omega_rn_0 = omega_rn_bf

        self.s = self.s + sigma_br * dt
        z = K * self.s + np.dot(I, (omega_bn - omega_rn_bf) - (self.omega_bn_0 - self.omega_rn_0))
        self.z_history.append(z)

        return (-K * sigma_br - np.dot(P, omega_br) + np.dot(I, omega_rn_dot_bf - np.cross(omega_bn, omega_rn_bf)) + np.cross(omega_bn, np.dot(I, omega_bn)) - L - np.dot(P, self.Ki * z))

class Gain_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):

        return (-K * sigma_br - np.dot(P, omega_br) + np.cross(omega_bn, np.dot(I, omega_bn)) - L)

class Saturated_Controller(Controllers):
    def __init__(self, umax):
        self.umax = umax

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):

        u_uns = (-K * sigma_br - np.dot(P, omega_br) + np.dot(I, omega_rn_dot_bf - np.cross(omega_bn, omega_rn_bf)) + np.cross(omega_bn, np.dot(I, omega_bn)) - L)
    
        if np.abs(u_uns).all() < self.umax :
            return u_uns
        else:
            return self.umax * np.sign(u_uns)
 
