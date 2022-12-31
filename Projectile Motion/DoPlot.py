# Modules
import matplotlib.pyplot as plt
import numpy as np
import constants as const


class DoPlot:
    def __init__(self, caption):

        # Figure
        self.size = [10, 10]
        self.fig = plt.figure(figsize=self.size)
        self.fig.canvas.set_window_title(caption)
        # plt.show()        # to see plots during the run_time

        # Sub figures
        self.xy = self.fig.add_subplot(1, 2, 1)  # x(t)
        self.R = self.fig.add_subplot(1, 2, 2)  # v(t)
        self.fig.tight_layout(pad=4.0)

        # Axis label
        self.xy.set_xlabel('$x (m)$')
        self.xy.set_ylabel('$y (m)$')
        self.R.set_xlabel('$θ₀ (degree)$')
        self.R.set_ylabel('$R (m)$')

    def y_vs_x(self, x, y, color):
        self.xy.scatter([x], [y], s=2, c=color)

    def range(self, r, theta, color):
        self.R.scatter([theta], [r], s=2, c=color)

    # To draw the Analytical plot for y vs. x
    def analytic_y_vs_x(self):
        t = 0
        y_exact = [0]
        x_exact = [0]
        while y_exact[-1] >= 0:
            x_exact.append(const.x_0 + const.vx_0 * t)
            y_exact.append(const.y_0 + const.vy_0 * t - 0.5 * const.g * (t ** 2))
            t += const.h
        self.xy.plot(x_exact, y_exact, c='g', alpha=0.5)

    # To draw the Analytical plot for Range vs. θ
    def analytic_range(self):
        t = 0
        Range = []
        for j in const.theta:
            y_exact = [0]
            x_exact = [0]
            vx = float(const.v_0 * np.cos(j))
            vy = float(const.v_0 * np.sin(j))
            # y_exact is an array, and we have to check the last element
            # to see whether it's positive or not
            while y_exact[-1] >= 0:
                x_exact.append(const.x_0 + vx * t)
                y_exact.append(const.y_0 + vy * t - 0.5 * const.g * (t ** 2))
                t += const.h
            # range is the biggest value of x_exact
            Range.append(x_exact[-1])
        self.R.plot(const.theta_deg, Range, c='g', alpha=0.5)
