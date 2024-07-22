import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import EllipseCollection, LineCollection

from simul import Simul


class Animate:
    def __init__(self, simul):
        self.simul = simul
        self.fig, self.ax = plt.subplots(figsize=(5, 5))  # initialise  graphics
        self.circles = EllipseCollection(widths=2*simul.sigma, heights=2*simul.sigma, angles=0, units='x',
                                         offsets=simul.position, transOffset=self.ax.transData)  # circles at position
        self.ax.add_collection(self.circles)

        self.segment = [[[0, 0], [0, simul.L], [simul.L, simul.L], [simul.L, 0], [0, 0]]]  # simulation cell
        self.line = LineCollection(self.segment, colors='#000000')  # draw square
        self.ax.add_collection(self.line)

        self.ax.set_xlim(left=-0.5, right=self.simul.L+0.5)  # plotting limits on screen
        self.ax.set_ylim(bottom=-0.5, top=self.simul.L+0.5)
        self._ani = 0

    def init(self):  # this is the first thing drawn to the screen
        self.circles.set_offsets(self.simul.position)

    def anim_step(self, m):  # m is the number of calls that have occurred to this function
        print('anim_step m = ', m)
        if m == 0:
            time.sleep(0.5)

        self.simul.md_step()  # perform simulation step
        self.circles.set_offsets(self.simul.position)  # update positions on screen

    def go(self, nframes):
        self._ani = animation.FuncAnimation(self.fig, func=self.anim_step, frames=nframes,
                                            repeat=False, interval=10, init_func=self.init)  # run animation
        plt.show()
