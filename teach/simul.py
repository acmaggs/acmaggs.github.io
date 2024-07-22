import numpy as np
import math


class Simul:
    """ 
    This is the prototype of the simulation code
    It moves the particles with at velocity, using a vector notation: numpy should be used.
    """
    def __init__(self, simul_time, sigma, L):
        np.seterr(all='ignore')  # remove errors in where statements
        self.position = np.array([[2., 2.], [5., 2.], [2., 5.], [5., 5.]])  # starting positions
        self.velocity = 3*np.random.normal(size=self.position.shape)  # random velocities
        self.l, self.m = np.triu_indices(self.position.shape[0], k=1)  # all pairs of indices between particles
        self.sigma = sigma  # particle radius
        self.simul_time = simul_time
        self.L = L

    def wall_time(self):
        first_collision_time = np.inf
        particle = np.inf
        direction = np.inf
        return first_collision_time, particle, direction
        # calculate time of first collision, particle involved and direction

    def pair_time(self):
        pass

    def md_step(self):
        print('Simul::md_step')
        ke_start = (self.velocity**2).sum()/2.   # starting kinetic energy

        pressure = -1
        current_time = 0
        condition_on_time_variables = False

        w_time, particle, direction = self.wall_time()

        while condition_on_time_variables:   # think about this
            # do something
            w_time, particle, direction = self.wall_time()  # update collisions times
            # update current_time

        # adapt the position update  as a function of your logic
        self.position += self.simul_time * self.velocity

        assert math.isclose(ke_start,  (self.velocity**2).sum()/2.)  # check that we conserve energy after all the collisions

        return pressure

    def __str__(self):   # this is used to print the position and velocity of the particles
        p = np.array2string(self.position)
        v = np.array2string(self.velocity)
        return 'pos= '+p+'\n'+'vel= '+v+'\n'
