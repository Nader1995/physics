# modules
import constants
from RD import rd
from LJforce import force


def acc(x, y, ax, ay):
    n = constants.N
    for i in range(n):
        ax[i] = 0
        ay[i] = 0

    for i in range(n - 1):
        for j in range(i + 1, n):
            dx = rd(x[i], x[j])
            dy = rd(y[i], y[j])

            fx, fy = force(dx, dy)
            ax[i] += fx
            ay[i] += fy
            ax[j] -= fx
            ay[j] -= fy

    return [ax, ay]
