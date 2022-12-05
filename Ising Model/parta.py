# ********************************************************************
# *                                                                  *
# * Ising, v 1.02                                                    *
# * Date: 13990605                                                   *
# * Nader Safarinia, Student ID: 981104                              *
# * Copyleft 🄯 2020 N. Safari. All left reserved!                    *
# *                                                                  *
# ********************************************************************

# modules
from matplotlib import cm
from celluloid import Camera
import numpy as np
import math
import random
from datetime import datetime
from eflip import e_flip
import constants
import matplotlib.pyplot as plt

# constants
J = constants.J
h = constants.h
L = constants.L
N = constants.N
runTime = constants.runTime
T = float(4)      # Temperature; a) T = 1.5, 2.0, 2.25, 4.0

# variables
N_1 = np.arange(0, runTime, dtype=int)           # Monte Carlo step per spin for part (a),
M = np.zeros(runTime)                            # To store magnetization for part (a)
s = np.ones((L, L))

# I use (a) and (b) to build the lattice
a = np.arange(0, N)
for i in range(N):
    a[i] = a[i] / L

b = np.arange(0, N)
for i in range(N):
    b[i] = b[i] % L

points = [a, b]
camera = Camera(plt.figure())

counter = 0
while counter < runTime:

    # To show the lattice
    num = 0
    c = np.zeros(N)
    for n in range(L):
        for m in range(L):
            # (c) contains arrays of -1 and 1,
            # each of which indicates a specific color
            c[num] = s[n][m]
            num += 1

    # c = 0 and c = 1 are red and blue colors
    c += 1
    colors = cm.rainbow(c)
    plt.scatter(*points, c=colors, s=10)
    camera.snap()

    # To calculate the magnetization
    for n in range(L):
        for m in range(L):
            M[counter] += s[n][m]
    M[counter] = M[counter] / N

    for k in range(N):
        random.seed(datetime.now())
        i = random.randint(0, L - 1)
        random.seed(datetime.now())
        j = random.randint(0, L - 1)
        e = e_flip(s, i, j)

        if e <= 0:
            s[i][j] = - s[i][j]
        else:
            random.seed(datetime.now())
            r = random.random()
            if r <= math.exp(-e/T):
                s[i][j] = - s[i][j]
            else:
                pass

    counter += 1
# To depict the result in an animation mode
anim = camera.animate(blit=True)
anim.save('scatter.mp4')

# To depict the magnetization in a plot
# plt.plot(N_1, M)
# plt.xlabel("T")
# plt.ylabel("m")
# plt.show()
