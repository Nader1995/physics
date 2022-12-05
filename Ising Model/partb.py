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
import matplotlib.pyplot as plt
import constants
import numpy as np
import math
import random
from eflip import e_flip
from datetime import datetime


# constants
h = constants.h
J = constants.J
T = float(10)      # Temperature; b) T = 10
L = constants.L
N = constants.N
epsilon = 0.05
runTime = constants.runTime

# variables
N_1 = np.arange(0, runTime, dtype=int)           # Monte Carlo step per spin for part (a),
M_1 = np.zeros(runTime)
M_2 = np.zeros(runTime)
s = np.zeros((L, L))

# To randomize the elements of (s)
for i in range(L):
    for j in range(L):
        x = random.random()
        if x >= 0.5:
            s[i][j] = 1
        else:
            s[i][j] = -1


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

ki = [0]
magnet = [0]
temp = [T]
while T > 0.5:

    counter = 0
    while counter < runTime:

#        num = 0
#        c = np.zeros(N)
#        for n in range(L):
#            for m in range(L):
#                c[num] = s[n][m]
#                num += 1

#        c += 1
#        colors = cm.rainbow(c)
#        plt.scatter(*points, c=colors, s=10)
#        camera.snap()

        for n in range(L):
            for m in range(L):
                M_1[counter] += s[n][m]
        M_2[counter] = (M_1[counter] / N) ** 2
        M_1[counter] = M_1[counter] / N

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
                if r <= math.exp(-e / T):
                    s[i][j] = - s[i][j]
                else:
                    pass

        counter += 1

    T = T * (1 - epsilon)
    ki.append(((sum(M_2)/runTime) - (sum(M_1)/runTime) ** 2) / T)
    magnet.append(sum(M_1)/runTime)
    temp.append(T)
    print("here!")
    print(T)

# anim = camera.animate(blit=True)
# anim.save('scatter1.mp4')

# plt.plot(temp, magnet)
# plt.xlabel("T")
# plt.ylabel("m")
# plt.show()

plt.plot(temp, ki)
plt.xlabel("T")
plt.ylabel("χ")
plt.show()
