#
# \file HalfEdge.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 04-Dec-2024
# Handles half-edges in the half-edge implementation of the mesh


class HalfEdge:

    def __init__(self, idx):
        self.idx = idx
        self.vfrom = None
        self.vto = None
        self.prev = None
        self.next = None
        self.pair = None
        self.face = None
        self.boundary = False
        self.myo = None
        self.tension = None
        # one can add properties here (e.g., myosin)
