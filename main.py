import numpy as np
from math import sin

def pendulum_rk4(theta_0, omega_0, dt, t_max):
    theta, omega = theta_0, omega_0
    t = 0
    theta_res, omega_res, t_res = [], [], []
    while t < t_max:
        theta_res.append(theta)
        omega_res.append(omega)
        t_res.append(t)
        k1_theta = omega
        k1_omega = -(g / L) * sin(theta)
        k2_theta = omega + (dt/2) * k1_omega
        k2_omega = -(g / L) * sin(theta + (dt/2) * k1_theta)
        k3_theta = omega + (dt/2) * k2_omega
        k3_omega = -(g / L) * sin(theta + (dt/2) * k2_theta)
        k4_theta = omega + dt * k3_omega
        k4_omega = -(g / L) * sin(theta + dt * k3_theta)
        theta = theta + (dt/6) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta)
        omega = omega + (dt/6) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega)
        t += dt
    return np.array(theta_res), np.array(omega_res), np.array(t_res)

g = 9.81 
L = 1.0 
theta_0 = np.pi/2 
omega_0 = 0 
dt = 0.01 
t_max = 10 
theta, omega, t = pendulum_rk4(theta_0, omega_0, dt, t_max)
print("Angle:")
print(theta)
print("Angular Speed:")
print(omega)
print("Time:")
print(t)
