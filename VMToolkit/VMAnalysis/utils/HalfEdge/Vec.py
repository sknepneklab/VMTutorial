#
# \file Vec.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 04-Dec-2024
# Handles operations on vectors in 2d


import numpy as np


class Vec:
    """Class that handles 2D vector algebra in periodic boundaries."""

    def __init__(self, r, box=None):
        """
          Construct a vector in 2d.
          Parameters
          ----------
            r : list
              list of two coordinates of a vector
            box : Box (optional)
              simulation box (used for periodic systems)
        """
        self.r = np.array(r, dtype=np.float64)
        self.box = box
        if box is not None:
            s = np.dot(self.box.invh, self.r)
            s -= np.rint(s)
            self.r = np.dot(self.box.h, s)

    def __add__(self, v):
        """
          Define vector addition.
          Parameters
          ----------
            v : Vec
              Vector to add to self
        """
        rr = self.r + v.r
        box = self.box if self.box is not None else v.box
        if self.box is not None or v.box is not None:
            si = np.dot(box.invh, self.r)
            sj = np.dot(box.invh, v.r)
            s = si + sj
            s -= np.rint(s)
            rr = np.dot(box.h, s)
        return Vec(rr, box)

    def __sub__(self, v):
        """
          Define vector subtraction.
          Parameters
          ----------
            v : Vec
              Vector to subtract from self
        """
        rr = self.r - v.r
        box = self.box if self.box is not None else v.box
        if self.box is not None or v.box is not None:
            si = np.dot(box.invh, self.r)
            sj = np.dot(box.invh, v.r)
            s = si - sj
            s -= np.rint(s)
            rr = np.dot(box.h, s)
        return Vec(rr, box)

    def __mul__(self, scale):
        """
          Define vector scaling (from left).
          Parameters
          ----------
            scale : float
              Scaling factor
        """
        return Vec(scale*self.r, self.box)

    def __rmul__(self, scale):
        """
          Define vector scaling (from right).
          Parameters
          ----------
            s : float
              Scaling factor
        """
        return Vec(scale*self.r, self.box)

    def __iadd__(self, v):
        """
          Increment the vector by another vector.
          Parameters
          ----------
            v : Vec
              Vector to increment by
        """
        rr = self + v
        self.r = rr.r
        return self

    def __isub__(self, v):
        """
          Decrement the vector by another vector.
          Parameters
          ----------
            v : Vec
              Vector to decrement by
        """
        rr = self - v
        box = self.box if self.box is not None else v.box
        if self.box is not None or v.box is not None:
            s = np.dot(box.invh, rr.r)
            s -= np.rint(s)
            rr.r = np.dot(self.box.h, s)
        self.r = rr.r
        return self
    
    def __neg__(self):
        """
            Unary minus operator.
        """
        return Vec(-self.r)

    def __repr__(self):
        """
          Return vector components as a tuple
        """
        return "({:.6f}, {:.6f})".format(self.r[0], self.r[1])

    def __str__(self):
        """
          Return vector components as a string for printing.
        """
        return "({:.6f}, {:.6f})".format(self.r[0], self.r[1])

    def dot(self, v):
        """
          Compute dot product between two vectors.
          Parameters
          ----------
            v : Vec
              Vector to dot product with
        """
        return np.dot(self.r, v.r)

    def outer(self, v):
        """
           Compute outer product between two vectors.
           Parameters
           ----------
             v : Vec
               Vector to compute outer product with   
        """
        return np.outer(self.r, v.r)

    def length(self):
        """ 
          Compute length of the vector.
        """
        return np.sqrt(np.sum(self.r**2))

    def unit(self):
        """ 
          Return a unit-length vector in the direction of this vector.
        """
        l = self.length()
        if l > 0.0:
            return (1.0/l)*self

    def rotate(self, phi):
        """
          Rotate the vector in plane.
          Parameter
          ---------
            phi : rotaton angle
        """
        R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        self.r = np.dot(R, self.r)
        if self.box is not None:
            s = np.dot(self.box.invh, self.r)
            s -= np.rint(s)
            self.r = np.dot(self.box.h, s)

    def to_list(self):
        """
          Return a list with x an y components of the vector.
        """
        return self.r.tolist()

    def ez_cross(self):
        """
          Return cross product of e_z with this vector.
        """
        return Vec([-self.r[1], self.r[0]], self.box)

    def cross(self, b):
        """
            Cross product of this vector with vector b.
            Since these are 2d vectors, this is a scalar.
        """
        return self.r[0]*b.r[1] - self.r[1]*b.r[0]
