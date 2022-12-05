# modules
import numpy as np
import random
import constants
from Acc import acc


def init(x, y, vx, vy, ax, ay):
    n = constants.N
    length = constants.L
    t0 = constants.T0
    sqn = int(round(np.sqrt(n)))
    delta = float(0.5 * length / np.sqrt(n))

    for i in range(n):
        j = int(i / sqn)
        k = int(i % sqn)
        x[i] = k * delta    # Initialize xᵢ
        y[i] = j * delta    # and yᵢ

    for i in range(n):
        vx[i] = random.random()     # By default the random number generator ...
        vy[i] = random.random()     # ... uses the current system time.

    sx = float(0)
    sy = float(0)

    for i in range(n):
        sx += vx[i]
        sy += vy[i]

    sx /= n
    sy /= n

    for i in range(n):
        vx[i] -= sx
        vy[i] -= sy

    s = 0
    for i in range(n):
        s += vx[i] ** 2 + vy[i] ** 2

    r = np.sqrt(2 * t0 * (n - 1) / s)

    for i in range(n):
        vx[i] *= r
        vy[i] *= r

    ax, ay = acc(x, y, ax, ay)

    return [x, y, vx, vy, ax, ay]
