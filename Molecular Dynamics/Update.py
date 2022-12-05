# modules
from Acc import acc
import constants
from PB import pb


def update(t, x, y, vx, vy, ax, ay, bx, by):
    n = constants.N
    h = constants.h
    for i in range(n):
        x[i] = pb(x[i] + h * vx[i] + 0.5 * ax[i] * (h ** 2))
        y[i] = pb(y[i] + h * vy[i] + 0.5 * ay[i] * (h ** 2))

    bx, by = acc(x, y, bx, by)

    for i in range(n):
        vx[i] += 0.5 * h * (ax[i] + bx[i])
        vy[i] += 0.5 * h * (ay[i] + by[i])

        ax[i] = bx[i]
        ay[i] = by[i]

    t += h
    return [t, x, y, vx, vy, ax, ay]
