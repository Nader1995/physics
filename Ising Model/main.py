# ********************************************************************
# *                                                                  *
# * Ising Model, v 1.10                                              *
# * Date: 13990603                                                   *
# * Nader Safarinia, Student ID: 981104                              *
# * Copyleft 🄯 2020 N. Safari. All left reserved!                    *
# *                                                                  *
# ********************************************************************

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from celluloid import Camera
import constants

# constants
L = constants.L
N = constants.N

a = np.arange(0, N)
for i in range(N):
    a[i] = a[i] / L

b = np.arange(0, N)
for i in range(N):
    b[i] = b[i] % L

numpoints = N
points = [a, b]
c = np.zeros(numpoints)
camera = Camera(plt.figure())
for _ in range(100):
    colors = cm.rainbow(c)
    c += 0.01
    plt.scatter(*points, c=colors, s=10)
    camera.snap()
anim = camera.animate(blit=True)
anim.save('scatter.mp4')

