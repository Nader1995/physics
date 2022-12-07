# Modules
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import constants


class DoPlot:
    def __init__(self, caption):

        # Figure
        self.size = [10, 10]
        self.fig = plt.figure(figsize=self.size)
        self.fig.canvas.set_window_title(caption)
        # plt.show()        # to see plots during the run_time

        # Sub figures
        self.x1 = self.fig.add_subplot(1, 2, 1)  # x(t)
        self.x2 = self.fig.add_subplot(1, 2, 2)  # x²(t)
        self.fig.tight_layout(pad=4.0)

        # Axis label
        # Here I change the axis label automatically
        self.x1.set_xlabel(f'$t (1 = {int(1/constants.h)} steps)$')
        self.x1.set_ylabel('$x(t) [m]$')
        self.x2.set_xlabel(f'$t (1 = {int(1/constants.h)} steps)$')
        self.x2.set_ylabel('$<x²(t)> [m²]$')

    # Plot x vs. t
    def plot_func1(self, t, x, color):
        self.x1.scatter([t], [x], s=2, c=color)

    # Plot x² vs. t, and the fit line
    def plot_func2(self, t, x):
        self.x2.scatter(t, x, s=2, c="blue")
        # To draw the fit line
        # The number determines the polynomial order
        coefs = poly.polyfit(t, x, 1)
        ffit = poly.polyval(t, coefs)
        self.x2.plot(t, ffit, c="red")
