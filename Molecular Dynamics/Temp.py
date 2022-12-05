# modules
import constants
from initialization import init
from Update import update
import numpy as np
import matplotlib.pyplot as plt


t = 0.0
h = constants.h
t_max = constants.t_max
temper = [0.0]
time = [0.0]

x = np.zeros(constants.N, dtype=float)  # The position,
y = np.zeros(constants.N, dtype=float)
vx = np.zeros(constants.N, dtype=float)  # velocity and
vy = np.zeros(constants.N, dtype=float)
ax = np.zeros(constants.N, dtype=float)  # acceleration of particles.
ay = np.zeros(constants.N, dtype=float)
bx = np.zeros(constants.N, dtype=float)  # The temporary acceleration of particles in (n+1) time step.
by = np.zeros(constants.N, dtype=float)

x, y, vx, vy, ax, ay = init(x, y, vx, vy, ax, ay)

while t <= t_max:
    s = 0
    for i in range(constants.N):
        s += vx[i] ** 2 + vy[i] ** 2
    temper.append(2 * s / (constants.N - 1))
    t, x, y, vx, vy, ax, ay = update(t, x, y, vx, vy, ax, ay, bx, by)
    time.append(t)

plt.title(" Temperature Fluctuation ")
plt.ylabel(" T [K] ")
plt.xlabel(" t [ps] ")
plt.plot(time, temper)
plt.savefig("TemperatureVsTime")
