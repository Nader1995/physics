# Modules
import numpy as np

# Constants
a = 6.5e-3  # K/m
alpha = 2.5
T_0 = 273.15 + 0.85  # According to 2020 records by NASA
Y_0 = 10 ** 4


def linear(y):
    return 1.2 - 8e-5 * y


def isothermal(y):
    return 1.2 * np.exp(-y/Y_0)


def adiabatic(y):
    return 1.2 * (1 - a * y / T_0) ** alpha
