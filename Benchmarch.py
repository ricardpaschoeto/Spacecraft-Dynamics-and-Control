from numpy import linalg
from scipy.spatial.transform import Rotation
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos


### Given constants and other configuration

mrp = np.array([0.1, 0.2, -0.1]).T
mrp_history = [mrp]
w = np.array([np.deg2rad(i) for i in [30, 10, -20]]).T
w_history = [w]
w0 = np.copy(w)
mrp0 = np.copy(mrp)

K = 5 # Nm
P = 10 * np.eye(3) #NMs
#P = np.diag([22.3607,19.3649,20.000])
#Ki = 0.005
#Ki = 0

I1 = 100
I2 = 75
I3 = 80 # kg m^2

I = np.array([[I1, 0, 0], [0, I2, 0],[0, 0, I3]])
#DL = np.array([0.5,-0.3,0.2]).T
DL = np.array([0,0,0]).T
L = np.array([0.5, -0.3, 0.2]).T

return_0 = False

###

### Utility functions

def tilde(x):
    x = np.squeeze(x)
    return np.array([[0, -x[2], x[1]],
                     [x[2], 0, -x[0]],
                     [-x[1], x[0], 0]
                     ])


def mrp_to_rotation_matrix(sigma):
    sigma_squared = np.inner(sigma, sigma)
    q0 = (1 - sigma_squared) / (1 + sigma_squared)
    q = [2 * sigma_i / ( 1 + sigma_squared) for sigma_i in sigma]
    q.extend([q0])
    return Rotation.from_quat(q).as_matrix().T

def rotmat_to_mrp(matrix):
    zeta = np.sqrt(np.trace(matrix) + 1)
    constant = 1 / (zeta**2 + 2 * zeta)
    s1 = constant * (matrix[1, 2] - matrix[2, 1])
    s2 = constant * (matrix[2, 0] - matrix[0, 2])
    s3 = constant * (matrix[0, 1] - matrix[1, 0])
    return np.array([s1, s2, s3])

def mrp_shadow(mrp):
    norm = np.linalg.norm(mrp) ** 2
    return np.array([-i / norm for i in mrp])


def target(tval):
    if return_0:
        return np.array([0, 0, 0])
    f = 0.05
    s1 = 0.2 * np.sin(f * tval)
    s2 = 0.3 * np.cos(f * tval)
    s3 = -0.3 * np.sin(f * tval)
    return np.array([s1, s2, s3])

def target_rate(tval):
    if return_0:
        return np.array([0, 0, 0])
    f = 0.05
    s1_dot = 0.2 * f * cos(f * tval)
    s2_dot = -0.3 * f * sin(f * tval)
    s3_dot = -0.3 * f * cos(f * tval)
    sigma_dot = np.array([s1_dot, s2_dot, s3_dot])
    sigma = target(tval)
    A = mrp_dot_matrix(sigma)
    w = 4 * np.dot(np.linalg.inv(A), sigma_dot)

    return w

def target_rate_rate(tval, dt):
    if return_0:
        return np.array([0, 0, 0])
    w1 = target_rate(tval)
    w2 = target_rate(tval - dt)
    return (w1 - w2)/dt

def mrp_dot(mrp, w):
    return 0.25 * np.dot(((1 - np.dot(mrp, mrp)) * np.eye(3) + 2 * tilde(mrp) + 2 * np.outer(mrp,  mrp)), w)

def mrp_dot_matrix(mrp):
    ss = np.dot(mrp, mrp)
    A = np.zeros((3,3))
    A[0, 0] = 1 - ss + 2 * mrp[0] **2
    A[1, 0] = 2*(mrp[1] * mrp[0] + mrp[2])
    A[2, 0] = 2*(mrp[2] * mrp[0] - mrp[1])
    A[0, 1] = 2*(mrp[0] * mrp[1] - mrp[2])
    A[1, 1] = 1 - ss + 2 * mrp[1] ** 2
    A[2, 1] = 2*(mrp[2] * mrp[1] + mrp[0])
    A[0, 2] = 2*(mrp[0] * mrp[2] + mrp[1])
    A[1, 2] = 2*(mrp[1] * mrp[2] - mrp[0])
    A[2, 2] = 1 - ss + 2 * mrp[2] ** 2
    return A


def control(t, dt, mrp, w):
    sigma_r_n = target(t)
    sigma_b_r = rotmat_to_mrp(np.dot(mrp_to_rotation_matrix(mrp), mrp_to_rotation_matrix(sigma_r_n).T))
    w_r_n = target_rate(t)
    w_r_n_dot = target_rate_rate(t, dt)
    DCM_b_r = mrp_to_rotation_matrix(sigma_b_r)
    w_r_n_body_frame = np.dot(DCM_b_r, w_r_n)
    w_r_n_dot_body_frame = np.dot(DCM_b_r, w_r_n_dot)
    w_b_r = w - w_r_n_body_frame

    # delta_w = w - w_r_n
    # delta_w0 = w0 - target_rate(0)
    # v0 = sigma_integral(t, sigma_b_r[0])
    # v1 = sigma_integral(t, sigma_b_r[1])
    # v2 = sigma_integral(t, sigma_b_r[2])
    # v_int = np.array([v0,v1,v2])
    # z = v_int + I @ (delta_w - delta_w0)
    # control is
    # -K*mrp_br - P * w_br + I * (w_rn_dot - w_bn X w_rn) + w_bn X Iw_bn
    u = -K * sigma_b_r - np.dot(P, w_b_r) + np.dot(I, w_r_n_dot_body_frame - np.cross(w, w_r_n_body_frame)) + np.cross(w, np.dot(I, w)) - L
    #u = -K * sigma_b_r - np.dot(P, w_b_r) + np.dot(I, w_r_n_dot_body_frame - np.cross(w, w_r_n_body_frame)) + np.cross(w, np.dot(I, w)) - np.dot(np.dot(P, Ki),z) - L
    #u = -K * sigma_b_r - np.dot(P, w_b_r) 
    #if linalg.norm(u) >= 1:
        #u = u*np.sign(w)

    return u

def wdot(t, dt, mrp, w):
    u = control(t, dt, mrp, w)
    # Iw_dot = -w X Iw + Q
    w_dot = np.dot(np.linalg.inv(I), (-np.cross(w, np.dot(I, w)) + u + DL + L))
    return w_dot

def integrand(sigma_b_r):
    return sigma_b_r

def sigma_integral(t, sigma_b_r):
    sigma_int = quad(integrand, 0, t)[0]

    return sigma_int


h = 0.01
time = 120
t_target = 70

tvec = np.linspace(0, time, int(time/h + 1))
t_target_vec = int(t_target/h + 1)
prev_t = 0
target_histories = [target(0)]
target_rate_history = [target_rate(0)]
error_history = [rotmat_to_mrp(np.dot(mrp_to_rotation_matrix(mrp), mrp_to_rotation_matrix(target(0)).T))]
#print(error_history)
for ti in tvec[1:]:
    dt = ti - prev_t
    prev_t = ti
    sigma_r_n = target(ti)
    target_histories.append(sigma_r_n)
    sigma_b_r = rotmat_to_mrp(np.dot(mrp_to_rotation_matrix(mrp), mrp_to_rotation_matrix(sigma_r_n).T))
    #print(sigma_b_r)
    error_history.append(sigma_b_r)
    w_r_n = target_rate(ti)
    target_rate_history.append(w_r_n)
    #w_r_n_dot = target_rate_rate(ti, dt)
    #w_b_r = w - w_r_n

    # calculate and apply dots
    mrp = mrp + mrp_dot(mrp, w) * dt
    w = w + wdot(ti, dt, mrp, w) * dt

    if np.dot(mrp, mrp) > 1:
        mrp = mrp_shadow(mrp)

    mrp_history.append(mrp)
    w_history.append(w)
    if ti % 25 == 0:
        print("Simulated {} seconds".format(ti))


mrp_history = np.array(mrp_history)
w_history = np.array(w_history)
target_histories = np.array(target_histories)
target_rate_history = np.array(target_rate_history)
error_history = np.array(error_history)
error_norm = [np.sqrt(i**2 + j**2 + k**2) for i, j, k in error_history]
print("Norm sigma_br at ", t_target, "s: ", error_norm[t_target_vec])

_, axs = plt.subplots(2,2)
axs[0,0].plot(tvec, mrp_history[:, 0], 'g')
axs[0,0].plot(tvec, target_histories[:, 0], 'g--')
axs[0,0].plot(tvec, mrp_history[:, 1], 'b')
axs[0,0].plot(tvec, target_histories[:, 1], 'b--')
axs[0,0].plot(tvec, mrp_history[:, 2], 'r')
axs[0,0].plot(tvec, target_histories[:, 2], 'r--')
axs[0,0].set_title('Attitude (sigma) history')
axs[0,0].legend(["sigma_1", "sigma_1_target", "sigma_2", "sigma_2_target", "sigma_3", "sigma_3_target"])
axs[0,0].grid()

axs[0,1].plot(tvec, w_history[:, 0], 'b')
axs[0,1].plot(tvec, target_rate_history[:, 0], 'b--')
axs[0,1].plot(tvec, w_history[:, 1], 'r')
axs[0,1].plot(tvec, target_rate_history[:, 1], 'r--')
axs[0,1].plot(tvec, w_history[:, 2], 'g')
axs[0,1].plot(tvec, target_rate_history[:, 2], 'g--')
axs[0,1].set_title("Rate (w) history")
axs[0,1].legend(["w_1", "w_1_target", "w_2", "w_2_target", "w_3", "w_3_target"])
axs[0,1].grid()

axs[1,0].plot(tvec, error_norm)
axs[1,0].set_title("Error Norm")
axs[1,0].grid()

axs[1,1].plot(tvec, error_history[:, 0], 'b')
axs[1,1].plot(tvec, error_history[:, 1], 'r')
axs[1,1].plot(tvec, error_history[:, 2], 'g')
axs[1,1].set_title("Sigma BR history")
axs[1,1].legend(["sigma_1", "sigma_2", "sigma_3"])
axs[1,1].grid()

plt.show()