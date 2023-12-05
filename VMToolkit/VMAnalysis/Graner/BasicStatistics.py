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
# \file BasicStatistics.py
# \author Silke Henkes, silke.henkes@bristol.ac.uk
# \revision Rastko Sknepnek, sknepnek@gmail.com
# \date 31-May-2020
# \brief Returns basic information for a mesh
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *


class BasicStatistics:

    def __init__(self, infile):
        ext = os.path.splitext(infile)[1]
        if ext != '.json' and ext != '.bz2':
            raise Exception('A JSON file with the mesh has to be provided.')
        self.mesh = Mesh()
        self.mesh.read(infile)
        # Set the total number of cells, ignoring the outer edge
        self.N = self.mesh.num_inner_faces

    def getAreas(self):
        areas = np.zeros(self.N)
        i = 0
        for f in self.mesh.faces:
            if not f.outer:
                areas[i] = f.area()
                i += 1
        return areas

    def getPerimeters(self):
        perimeters = np.zeros(self.N)
        i = 0
        for f in self.mesh.faces:
            if not f.outer:
                perimeters[i] = f.perimeter()
                i += 1
        return perimeters

    def getPos(self):
        positions = np.zeros((self.N, 2))
        i = 0
        for f in self.mesh.faces:
            if not f.outer:
                positions[i, :] = f.rc(self.mesh.box).to_list()
                i += 1
        return positions

    # tension on individual edges for set of faces
    def getEdgeTension(self, facelist):
        tensions = []
        for f in self.mesh.faces:
            if f.idx in facelist:
                he, first = f.he, f.he
                while True:
                    if he.tension == None:
                        raise Exception('Tension on half-edges not defined.')
                    else:
                        tensions.append(he.tension)
                    he = he.next
                    if he.idx == first.idx:
                        break
        return tensions

    def getTension(self):
        tensions = np.zeros(self.N)
        i = 0
        for f in self.mesh.faces:
            if not f.outer:
                he, first = f.he, f.he
                while True:
                    if he.tension == None:
                        raise Exception('Tension on half-edges not defined.')
                    else:
                        tensions[i] += he.tension
                    he = he.next
                    if he.idx == first.idx:
                        break
                i += 1
        return tensions

    # myosin on individual edges for set of faces
    def getEdgeMyosin(self, facelist):
        myosins = []
        for f in self.mesh.faces:
            if f.idx in facelist:
                he, first = f.he, f.he
                while True:
                    if he.myo == None:
                        raise Exception('Tension on half-edges not defined.')
                    else:
                        myosins.append(he.myo)
                    he = he.next
                    if he.idx == first.idx:
                        break
        return myosins

    def getMyosin(self,):
        myosins = np.zeros(self.N)
        i = 0
        for f in self.mesh.faces:
            if not f.outer:
                he, first = f.he, f.he
                while True:
                    if he.myo == None:
                        raise Exception('Myosin on half-edges not defined.')
                    else:
                        myosins[i] += he.myo
                    he = he.next
                    if he.idx == first.idx:
                        break
                i += 1
        return myosins
