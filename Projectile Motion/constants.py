# modules
import numpy as np

color = ['r', 'b', 'k']
m = float(2.9)  # mass [kg]
A = float(7.8e-3)  # cross-section [m²]
v_0 = float(400)  # The initial velocity  [m/s]
g = float(9.8)  # gravitational acceleration [m/s²]
C = [0, float(0.05)]  # The dimensionless constant of the drag force
# It's an array, which enables us to include and exclude the drag force

h = float(1 / 8)
x_0 = float(0)
y_0 = float(0)
theta_0 = float(45 * np.pi / 180)
# somewhere in the code, it's useful to have θ in degree
theta_deg = np.linspace(0, 90, 50, endpoint=True)
# and somewhere it's not
theta = theta_deg * np.pi / 180
vx_0 = float(v_0 * np.cos(theta_0))
vy_0 = float(v_0 * np.sin(theta_0))

