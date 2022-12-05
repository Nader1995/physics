# modules
import constants


def pb(x):
    length = constants.L
    if length >= x >= 0:
        return float(x)
    elif x > length:
        return float(x - length)
    elif x < 0:
        return float(x + length)
