#
# \file TensionTensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 09-Dec-2023
# \brief Computes a tensor from tensions on edges and produces a VTP file for it
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *
from .Tensor import Tensor


class TensionTensor:

    def __init__(self, infile):
        ext = os.path.splitext(infile)[1]
        if ext != '.json' and ext != '.bz2':
            raise Exception('A JSON file with the mesh has to be provided.')
        self.mesh = Mesh()
        self.mesh.read(infile)
        self.Tension = Tensor(self.mesh.num_inner_faces)
        idx = 0
        for f in self.mesh.faces:
            if not f.outer:
                ev = f.edge_vectors()
                i = 0
                he, first = f.he, f.he
                while True:
                    v = ev[i]
                    tension = he.tension
                    if tension == None:
                        raise Exception('Tension on half-edges not defined.')
                    self.Tension.T[idx, :, :] += tension * \
                        np.array(v.outer(v))/v.length()
                    i += 1
                    he = he.next
                    if he.idx == first.idx:
                        break
                self.Tension.T[idx, :, :] /= f.area()
                idx += 1

    def plot_vtk_tensor(self, filename):
        self.Tension.plot_vtk_tensor(filename, self.mesh, 'Tension_tensor')
