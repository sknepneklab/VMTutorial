#
# \file TTensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 09-Dec-2023
# Computes the T tensor (Graner statistics) and produces a VTP file for it
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *
from .Tensor import Tensor


class TTensor:

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
        self.T = Tensor(self.mesh_2.num_inner_faces)
        if ('time' in self.mesh_1.system) and ('time' in self.mesh_2.system):
            dt = self.mesh_2.system['time'] - self.mesh_1.system['time']
        else:
            dt = 1.0
        fidx = 0
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
        for (f_1, f_2) in zip(self.faces_1, self.faces_2):
            if not (f_1.outer or f_2.outer):
                neigh_1, neigh_2 = f_1.neighbours(), f_2.neighbours()
                # compute l(t) and l(t+dt)
                nv_1 = f_1.neigh_vectors(self.mesh_1.box)
                nv_2 = f_2.neigh_vectors(self.mesh_2.box)
                # links that appear after a T1 transition
                appear = list(set(neigh_2) - set(neigh_1))
                # links that dissapear in T1 transition
                disappear = list(set(neigh_1) - set(neigh_2))
                ma = np.zeros((2, 2))
                md = np.zeros((2, 2))
                if len(appear) > 0:
                    for a in appear:
                        idx = neigh_2.index(a)
                        l = nv_2[idx]
                        ma += np.array(l.outer(l))
                if len(disappear) > 0:
                    for d in disappear:
                        idx = neigh_1.index(d)
                        l = nv_1[idx]
                        md += np.array(l.outer(l))
                # Do not keep the time here for now
                # Use central average of contact numbers again
                zav = 0.5*(len(neigh_1) + len(neigh_2))
                self.T.T[fidx, :, :] = (1/dt)*(ma - md)/zav
                fidx += 1

    def plot_vtk_tensor(self, filename):
        self.T.plot_vtk_tensor(filename, self.mesh_2, 'T_tensor')

    def plot_vtk_ellipse(self, filename, N = 20, scale = 1.0):
        self.T.plot_vtk_ellipse(filename, self.mesh_2, N = N, scale = scale)

    def plot_vtk_lines(self, filename, scale = 1.0):
        self.T.plot_vtk_lines(filename, self.mesh_2, scale = scale)