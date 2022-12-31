# ********************************************************************
# *                                                                  *
# * Projectile Motion, v 1.01                                        *
# * Date: 13990322                                                   *
# * Nader Safarinia, Student ID: 981104                              *
# * Copyleft 🄯 2020 N. Safari. All left reserved!                    *
# *                                                                  *
# ********************************************************************

# modules
from time import process_time
import matplotlib.pyplot as plt
import numpy as np
import constants as const
from DoPlot import DoPlot
import AirDensity as Air

# Constants
h = const.h
theta = const.theta
g = const.g
C = const.C
color = const.color
A = const.A
m = const.m

# Header
print("Projectile Motion, v 1.01")
print("Date: 13990323")
print("Parameters:")
print("x(0) = x₀ = %.1f\t\t y(0) = y₀ = %.1f" % (const.x_0, const.y_0))
print("Δt = %.e s \t\t v(0) = v₀ = %.1f m/s" % (h, const.v_0))
print("C = %.2f  \t\t\t m = %.1f kg" % (C[1], m))
print("A = %.4f m²\t\t g = %.1f kg/m^3" % (A, g))
print("h = %.4f " % h)


# Main part of the code
start_time = process_time()  # Holds the start time of the program.

# Init Plot
plt.ion()

# Projectile Motion, Type I: Linear Approximation for Air Density
# ==============================================================

# Define the Object named 'linear'
linear = DoPlot('Projectile Motion, Type I: Linear Approximation in Air Density')
linear.xy.set_title('Y vs. X')
linear.R.set_title('Range vs. θ')
linear.xy.text(900, -100, "$Red$" + ' ' + '$ Dot:$' + ' ' + '$No$' + ' ' + '$Air$\n' +
               "$Blue$" + ' ' + '$ Dot:$' + ' ' + '$With$' + ' ' + '$Air$\n' +
               '$Green$' + ' ' + '$ Line:$' + ' ' + '$Exact$' + ' ' + '$ Solution$', style='oblique',
               bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})
linear.R.text(20, 70, "$Red$" + ' ' + '$ Dot:$' + ' ' + '$No$' + ' ' + '$Air$\n' +
              "$Blue$" + ' ' + '$ Dot:$' + ' ' + '$With$' + ' ' + '$Air$\n' +
              '$Green$' + ' ' + '$ Line:$' + ' ' + '$Exact$' + ' ' + '$ Solution$', style='oblique',
              bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})

# First things first, let's plot the analytic graph for y vs. x.
linear.analytic_y_vs_x()

for c in C:
    x = const.x_0
    y = const.y_0
    vx = const.vx_0
    vy = const.vy_0
    while y >= 0:
        x += h * vx
        y += h * vy
        v = np.sqrt(vx ** 2 + vy ** 2)
        # Here we can choose the specific approximation
        rho = Air.linear(y)
        psi = (0.5 * c * rho * A * v) / m
        vx -= h * psi * vx
        vy -= h * (g + psi * vy)
        linear.y_vs_x(x, y, color[C.index(c)])

# To plot the analytic graph for Range vs. θ
linear.analytic_range()

for c in C:
    for j in theta:  # j goes from 0 to 90 degree
        x = const.x_0
        y = const.y_0
        vx = float(const.v_0 * np.cos(j))
        vy = float(const.v_0 * np.sin(j))
        while y >= 0:
            x += h * vx
            y += h * vy
            v = np.sqrt(vx ** 2 + vy ** 2)
            rho = Air.linear(y)
            psi = (0.5 * c * rho * A * v) / m
            vx -= h * psi * vx
            vy -= h * (g + psi * vy)
        # To plot Range vs. θ, we turn j into degree
        linear.range(x, j * 180 / np.pi, color[C.index(c)])

# Projectile Motion, Type II: Isothermal Approximation for Air Density
# ==============================================================

isothermal = DoPlot('Projectile Motion, Type II: Isothermal Approximation for Air Density')
isothermal.xy.set_title('Y vs. X')
isothermal.R.set_title('Range vs. θ')
isothermal.xy.text(900, -100, "$Red$" + ' ' + '$ Dot:$' + ' ' + '$No$' + ' ' + '$Air$\n' +
                   "$Blue$" + ' ' + '$ Dot:$' + ' ' + '$With$' + ' ' + '$Air$\n' +
                   '$Green$' + ' ' + '$ Line:$' + ' ' + '$Exact$' + ' ' + '$ Solution$', style='oblique',
                   bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})
isothermal.R.text(20, 70, "$Red$" + ' ' + '$ Dot:$' + ' ' + '$No$' + ' ' + '$Air$\n' +
                  "$Blue$" + ' ' + '$ Dot:$' + ' ' + '$With$' + ' ' + '$Air$\n' +
                  '$Green$' + ' ' + '$ Line:$' + ' ' + '$Exact$' + ' ' + '$ Solution$', style='oblique',
                  bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})

isothermal.analytic_y_vs_x()

for c in C:
    x = const.x_0
    y = const.y_0
    vx = const.vx_0
    vy = const.vy_0
    while y >= 0:
        x += h * vx
        y += h * vy
        v = np.sqrt(vx ** 2 + vy ** 2)
        rho = Air.isothermal(y)
        psi = (0.5 * c * rho * A * v) / m
        vx -= h * psi * vx
        vy -= h * (const.g + psi * vy)
        isothermal.y_vs_x(x, y, color[C.index(c)])

isothermal.analytic_range()

for c in C:
    for j in theta:
        x = const.x_0
        y = const.y_0
        vx = float(const.v_0 * np.cos(j))
        vy = float(const.v_0 * np.sin(j))
        while y >= 0:
            x += h * vx
            y += h * vy
            v = np.sqrt(vx ** 2 + vy ** 2)
            rho = Air.isothermal(y)
            psi = (0.5 * c * rho * A * v) / m
            vx -= h * psi * vx
            vy -= h * (const.g + psi * vy)
        isothermal.range(x, j * 180 / np.pi, color[C.index(c)])

# Projectile Motion, Type III: adiabatic Approximation for Air Density
# Here we also include Type I and Type II, just for the sake of comparision
# ==============================================================

adiabatic = DoPlot('Projectile Motion, Type I, TypeII and Type III Comparision')
adiabatic.xy.set_title('Y vs. X')
adiabatic.R.set_title('Range vs. θ')
adiabatic.xy.text(900, 0, "$Red$" + ' ' + '$ Dot:$' + ' ' + '$Linear$' + ' ' + '$Approximation$\n' +
                  "$Blue$" + ' ' + '$ Dot:$' + ' ' + '$Isothermal$' + ' ' + '$Approximation$\n' +
                  '$Black$' + ' ' + '$ Dot:$' + ' ' + '$Adiabatic$' + ' ' + '$Approximation$', style='oblique',
                  bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})
adiabatic.R.text(10, 10, "$Red$" + ' ' + '$ Dot:$' + ' ' + '$Linear$' + ' ' + '$Approximation$\n' +
                 "$Blue$" + ' ' + '$ Dot:$' + ' ' + '$Isothermal$' + ' ' + '$Approximation$\n' +
                 '$Black$' + ' ' + '$ Dot:$' + ' ' + '$Adiabatic$' + ' ' + '$Approximation$', style='oblique',
                 bbox={'facecolor': 'red', 'alpha': 0.3, 'pad': 5})

# Here we are not interested in analytical answers or no_air case anymore
# adiabatic.analytic_y_vs_x()

for i in range(3):
    x = const.x_0
    y = const.y_0
    vx = const.vx_0
    vy = const.vy_0
    while y >= 0:
        x += h * vx
        y += h * vy
        v = np.sqrt(vx ** 2 + vy ** 2)
        # In this way, we accumulate all the approximations in ρ
        rho = [Air.linear(y), Air.isothermal(y), Air.adiabatic(y)]
        # We are not interested in C = 0 ( no_air case ) anymore. Hence , we write C[1]
        psi = (0.5 * C[1] * rho[i] * A * v) / m
        vx -= h * psi * vx
        vy -= h * (g + psi * vy)
        adiabatic.y_vs_x(x, y, color[i])

# adiabatic.analytic_range()

for i in range(3):
    for j in theta:
        x = const.x_0
        y = const.y_0
        vx = float(const.v_0 * np.cos(j))
        vy = float(const.v_0 * np.sin(j))
        while y >= 0:
            x += h * vx
            y += h * vy
            v = np.sqrt(vx ** 2 + vy ** 2)
            rho = [Air.linear(y), Air.isothermal(y), Air.adiabatic(y)]
            psi = (0.5 * C[1] * rho[i] * A * v) / m
            vx -= h * psi * vx
            vy -= h * (g + psi * vy)
        adiabatic.range(x, j * 180 / np.pi, color[i])

executing_time = 1000 * (process_time() - start_time)
print("\nFinish in %.1f ms!" % executing_time)

plt.show(block=True)
