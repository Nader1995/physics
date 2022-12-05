# modules
import numpy as np


def force(dx, dy):
    dr = np.sqrt(dx ** 2 + dy ** 2)
    g = 24 * (2 * (dr ** -12) - (dr ** -6)) * (dr ** -2)
    fx = -dx * g
    fy = -dy * g
    return [fx, fy]
