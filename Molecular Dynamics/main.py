# *****************************************************************************
# *                                                                           *
# * Molecular Dynamics - part I, v 1.20                                       *
# * Date: 13990406                                                            *
# * Nader Safarinia, Student ID: 981104                                       *
# * Copyleft 🄯 2020 N. Safari. All left reserved!                             *
# * Here I got help from the code written in the coming link:                 *
# * https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot *
# *                                                                           *
# *****************************************************************************

# Modules
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import constants
from initialization import init
from Update import update


# To run the MD code without time limit
class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""

    def __init__(self, numpoints=constants.N):
        self.numpoints = numpoints
        self.stream = self.data_stream()

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.upgrade, interval=5,
                                           init_func=self.setup_plot, blit=True)

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        x, y = next(self.stream).T
        self.scat = self.ax.scatter(x, y, vmin=0, vmax=1,
                                    cmap="jet", edgecolor="k")
        self.ax.axis([0, constants.L, 0, constants.L])
        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

    def data_stream(self):
        """Generate the coordinates, velocities and accelerations."""
        x = np.zeros(constants.N, dtype=float)  # The position,
        y = np.zeros(constants.N, dtype=float)
        vx = np.zeros(constants.N, dtype=float)  # velocity and
        vy = np.zeros(constants.N, dtype=float)
        ax = np.zeros(constants.N, dtype=float)  # acceleration of particles.
        ay = np.zeros(constants.N, dtype=float)
        bx = np.zeros(constants.N, dtype=float)  # The temporary acceleration of particles in (n+1) time step.
        by = np.zeros(constants.N, dtype=float)
        x, y, vx, vy, ax, ay = init(x, y, vx, vy, ax, ay)
        # We feed the x and y coordinates to xy variable
        xy = np.zeros((constants.N, 2))

        for i in range(constants.N):
            xy[i][0] = x[i]
            xy[i][1] = y[i]
        t = float(0)
        # There is no time limit here, none.
        while True:
            t, x, y, vx, vy, ax, ay = update(t, x, y, vx, vy, ax, ay, bx, by)
            for i in range(constants.N):
                xy[i][0] = x[i]
                xy[i][1] = y[i]

            yield np.c_[xy[:, 0], xy[:, 1]]

    def upgrade(self, i):
        """Update the scatter plot."""
        data = next(self.stream)

        # Set x and y data...
        self.scat.set_offsets(data[:, :2])

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,


if __name__ == '__main__':
    a = AnimatedScatter()
    plt.show()
