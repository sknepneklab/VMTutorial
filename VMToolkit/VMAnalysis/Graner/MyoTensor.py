###########################################################################
#
#  Copyright (C) 2017, 2018 University of Dundee
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
# \file MyoTensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 25-May-2020
# \brief Computes a tensor from myosin on edges and produces a VTP file for it
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *
from .Tensor import Tensor


class MyoTensor:

    def __init__(self, infile):
        ext = os.path.splitext(infile)[1]
        if ext != '.json' and ext != '.bz2':
            raise Exception('A JSON file with the mesh has to be provided.')
        self.mesh = Mesh()
        self.mesh.read(infile)
        self.Myo = Tensor(self.mesh.num_inner_faces)
        idx = 0
        for f in self.mesh.faces:
            if not f.outer:
                ev = f.edge_vectors()
                i = 0
                he, first = f.he, f.he
                while True:
                    v = ev[i]
                    myo = he.myo
                    if myo == None:
                        raise Exception('Myosin on half-edges not defined.')
                    self.Myo.T[idx, :, :] += myo * \
                        np.array(v.outer(v))/v.length()
                    i += 1
                    he = he.next
                    if he.idx == first.idx:
                        break
                self.Myo.T[idx, :, :] /= f.area()
                idx += 1

    def plot_vtk_tensor(self, filename):
        self.Myo.plot_vtk_tensor(filename, self.mesh, 'Myo_tensor')
