# modules
from period import period
import constants

# constants
J = constants.J
h = constants.h


def e_flip(s, i, j):

    result = s[i][period(j+1)] + s[i][period(j-1)] + s[period(i+1)][j] + s[period(i-1)][j]
    return 2 * s[i][j] * (h + J * result)
