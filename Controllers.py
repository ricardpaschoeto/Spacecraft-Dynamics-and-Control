import numpy as np

class Controllers:

    def __init__(self):
        self.u_history = []

    def control(self):
        pass

class Full_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):

        u = (-K * sigma_br - np.dot(P, omega_br) + np.dot(I, omega_rn_dot_bf - np.cross(omega_bn, omega_rn_bf)) + np.cross(omega_bn, np.dot(I, omega_bn)) - L)
        self.u_history.append(u)

        return u

class Proportional_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):
        
        u = (-K * sigma_br - np.dot(P, omega_br))
        self.u_history.append(u)

        return u 

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

        u = (-K * sigma_br - np.dot(P, omega_br) + np.dot(I, omega_rn_dot_bf - np.cross(omega_bn, omega_rn_bf)) + np.cross(omega_bn, np.dot(I, omega_bn)) - L - np.dot(P, self.Ki * z))
        self.u_history.append(u)

        return 

class Gain_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):

        u = (-K * sigma_br - np.dot(P, omega_br) + np.cross(omega_bn, np.dot(I, omega_bn)) - L)
        self.u_history.append(u)

        return u

class Saturated_Controller(Controllers):

    def __init__(self, umax):
        self.umax = umax
        super().__init__()

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):

        u_temp = np.zeros((3,))
        u_uns = (-K * sigma_br - np.dot(P, omega_br) + np.dot(I, omega_rn_dot_bf - np.cross(omega_bn, omega_rn_bf)) + np.cross(omega_bn, np.dot(I, omega_bn)) - L)
    
        for ii in range(len(u_uns)):
            if np.abs(u_uns[ii]) <= self.umax[ii]:
                u_temp[ii] = u_uns[ii]
            else:
                u_temp[ii] = self.umax[ii] * np.sign(u_uns[ii]) 
        self.u_history.append(u_temp)

        return u_temp

class CLD_Controller(Controllers):

    def control(self, K, P, L, I, sigma_br, omega_br, omega_bn, omega_rn_bf, omega_rn_dot_bf, t = 0, dt = 0.):

        p1 = - np.dot(I, np.dot(P,omega_bn))
        p2 = - np.dot(np.outer(omega_bn, omega_bn) + np.dot(4*K/(1 + np.dot(sigma_br, sigma_br)) - np.dot(omega_bn,omega_bn)/2, np.eye(3,3)), sigma_br)
        p3 = np.dot(I, p2)
        u = p1 + p3 + np.cross(omega_bn, np.dot(I, omega_bn))
        self.u_history.append(u)
 
        return u