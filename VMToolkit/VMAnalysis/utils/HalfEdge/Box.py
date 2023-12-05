#
# \file Box.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 04-Dec-2024
# Handles simulation box data
#

import numpy as np


class Box:

    def __init__(self, a, b=None):
        """
        Construct simulation box.
          Parameters
          ----------
          a : list 
             Components of the a vector of the simulation box a = [ax, ay]
          b : list
             Components of the a vector of the simulation box a = [ax, ay]
          Note
          ----
            Simulation box is centred at (0,0).
        """
        has_b = None
        if type(a) is not list:
            raise ValueError('a vector has to be a list.')
        if len(a) != 2:
            raise ValueError('a vector has to be give with two components.')
        if b is not None:
            if (type(b) is list) or (len(b) == 2):
                has_b = True

        self.a = np.array(a)
        if has_b:
            self.b = np.array(b)
        else:
            self.b = np.array([a[1], a[0]])

        self.h = np.vstack((self.a, self.b)).T
        self.invh = np.linalg.inv(self.h)
