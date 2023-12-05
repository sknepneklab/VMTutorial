# \file Face.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 04-Dec-2024
# Handles vertices in the half-edge implementation of the mesh

from .Vec import Vec


class Vertex:

    def __init__(self, idx, r):
        self.idx = idx
        self.r = r
        self.type = None
        self.boundary = False
        self.constraint = None
        self.erased = False
        self.he = None
