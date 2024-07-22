from simul import Simul
from animate import Animate
import numpy as np


def main():
    np.random.seed(1)  # set random numbers to be always the same
    simulation = Simul(simul_time=0.008, sigma=1., L=7)  # sigma particle radius # L box size
    print(simulation.__doc__)  # print the documentation from the class

    animate = Animate(simulation)
    animate.go(nframes=75)  # number of animation steps
    print(simulation)  # print last configuration to screen


if __name__ == '__main__':
    main()
