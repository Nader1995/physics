# ********************************************************************
# *                                                                  *
# * Random Walk, v 1.00                                              *
# * Date: 13990408                                                   *
# * Nader Safarinia, Student ID: 981104                              *
# * Copyleft 🄯 2020 N. Safari. All left reserved!                    *
# *                                                                  *
# ********************************************************************

# modules
import numpy as np
import random
from DoPlot import DoPlot
import matplotlib.pyplot as plt
import constants

# constants
color = constants.color
m = constants.m
t_max = constants.t_max
h = constants.h
counter = constants.counter


# Init Plot
plt.ion()

# Define the object to do plot
obj = DoPlot(" Balanced Random Walk")

s2 = np.zeros(counter, dtype=float)
time = np.linspace(0, t_max, counter, endpoint=True)

for q in range(m):
    x = 0
    t = 0
    for i in range(counter):
        # to plot the result for two values of q
        if q is 1:
            # simply plot the result
            obj.plot_func1(t, x, color[0])
        if q is 2:
            obj.plot_func1(t, x, color[1])
        r = random.random()
        if r < 0.5:
            # Including two parts of exercise in one block
            # x += 0.2
            x += 0.1
        else:
            x -= 0.1

        s2[i] += x ** 2
        t += h

# we need to separate plotting s2 from plotting x because ....
# we 're gonna use " poly.polyfit" function to plot a fit line
for i in range(counter):
    s2[i] /= m

# plot x² and the fit line
obj.plot_func2(time, s2)

plt.show(block=True)
