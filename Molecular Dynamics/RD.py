# modules
import constants


def rd(x1, x2):
    length = constants.L
    d = x2 - x1
    if length/2 >= d >= -length/2:
        return d
    elif d > length/2:
        return d - length
    elif d < -length/2:
        return d + length
