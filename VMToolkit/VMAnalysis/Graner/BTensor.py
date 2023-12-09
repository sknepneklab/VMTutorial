# \file BTensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 09-Dec-2023
# Computes the B tensor (Graner statistics) and produces a VTP file for it
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *
from .Tensor import Tensor


class BTensor:

    def __init__(self, frame_1, frame_2, dt=1.0):
        ext_1 = os.path.splitext(frame_1)[1]
        ext_2 = os.path.splitext(frame_2)[1]
        if (ext_1 != '.json' and ext_1 != '.bz2') or (ext_2 != '.json' and ext_2 != '.bz2'):
            raise Exception(
                'JSON files with two consecutive frames of the mesh state have to be provided.')
        self.mesh_1 = Mesh()
        self.mesh_2 = Mesh()
        self.mesh_1.read(frame_1)
        self.mesh_2.read(frame_2)
        self.B = Tensor(self.mesh_2.num_inner_faces)
        self.C = Tensor(self.mesh_2.num_inner_faces)
        if self.mesh_1.system['time'] and self.mesh_2.system['time']:
            dt = self.mesh_2.system['time'] - self.mesh_1.system['time']

        self.faces_1 = []
        self.faces_2 = []
        for f_1 in self.mesh_1.faces:
            if not f_1.outer:
                if 'original_id' not in f_1.params:
                    self.faces_1 = self.mesh_1.faces 
                    self.faces_2 = self.mesh_2.faces 
                    break
                f_2 = self.mesh_2.get_original_face(f_1.params['original_id'])
                if f_2 != None:
                    if not (f_1.erased or f_2.erased):
                        self.faces_1.append(f_1)
                        self.faces_2.append(f_2)
        # Implicit assumption that faces are not relabeled here
        idx = 0
        for (f_1, f_2) in zip(self.faces_1, self.faces_2):
            if not (f_1.outer or f_2.outer):
                # Number of links of either cells
                try:
                    neigh_1, neigh_2 = f_1.neighbours(), f_2.neighbours()
                except:
                    print(f_1.idx, f_1.erased, f_2.idx, f_2.erased)
                # Conserved links
                neigh = list(set(neigh_1) & set(neigh_2))
                nv_1 = f_1.neigh_vectors(self.mesh_1.box)
                nv_2 = f_2.neigh_vectors(self.mesh_2.box)
                for n in neigh:
                    idx_1, idx_2 = neigh_1.index(n), neigh_2.index(n)
                    # l vectors at two consecutive time steps
                    l_1, l_2 = nv_1[idx_1], nv_2[idx_2]
                    l = 0.5*(l_1 + l_2)
                    dl = l_2 - l_1
                    # Adding up the conserved links
                    self.C.T[idx, :, :] += np.array(l.outer(dl))
                # Need to normalise by the total number of links. Per cell, that comes out to (symmetrised) 0.5(len(neigh_1) + len(neigh_2))
                zav = 0.5*(len(neigh_1) + len(neigh_2))
                self.C.T[idx, :, :] = (1/dt)*self.C.T[idx, :, :]/zav
                self.B.T[idx, :, :] = self.C.T[idx, :, :] + \
                    self.C.T[idx, :, :].T
                idx += 1

    def plot_vtk_tensor(self, filename):
        self.B.plot_vtk_tensor(filename, self.mesh_2, 'B_tensor')
    
    def plot_vtk_ellipse(self, filename, N = 20, scale = 1.0):
        self.B.plot_vtk_ellipse(filename, self.mesh_2, N = N, scale = scale)

    def plot_vtk_lines(self, filename, scale = 1.0):
        self.B.plot_vtk_lines(filename, self.mesh_2, scale = scale)
