###########################################################################
#
#  Copyright (C) 2017, 2018, 2019, 2020 University of Dundee
#  All rights reserved.
#
#  This file is part of AJM (Active Junction Model) program.
#
#  AJM is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  AJM is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ############################################################################

#
# \file VTensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 21-Jul-2020
# Computes the V tensor (Graner statistics) and produces a VTP file for it
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *
from .Tensor import Tensor
from .MTensor import MTensor
from .BTensor import BTensor


class VTensor:

    def __init__(self, frame_1, frame_2, dt=1.0):
        ext_1 = os.path.splitext(frame_1)[1]
        ext_2 = os.path.splitext(frame_2)[1]
        if (ext_1 != '.json' and ext_1 != '.bz2') or (ext_2 != '.json' and ext_2 != '.bz2'):
            raise Exception(
                'JSON files with two consecutive frames of the mesh state have to be provided.')
        self.M = MTensor(frame_1)
        self.B = BTensor(frame_1, frame_2)
        self.V = Tensor(self.B.B.N)
        self.Omega = Tensor(self.B.B.N)
        for i in range(self.B.B.N):
            if self.M.M.N != self.B.B.N:
                midx = self.M.mesh.faces[self.B.mesh_2.faces[i].params["original_id"]].idx
            else:
                midx = i
            invM = np.linalg.inv(self.M.M.T[midx, :, :])
            C = self.B.C.T[i, :, :]
            term_1, term_2 = np.dot(invM, C), np.dot(C.T, invM)
            self.V.T[i, :, :] = 0.5*(term_1 + term_2)
            self.Omega.T[i, :, :] = 0.5*(term_1 - term_2)

    def plot_vtk_tensor(self, filename, plot_omega=False):
        self.V.plot_vtk_tensor(filename, self.B.mesh_2, 'V_tensor')
        if plot_omega:
            self.Omega.plot_vtk_tensor(
                'omega_'+filename, self.B.mesh_2, 'Omega_tensor')
