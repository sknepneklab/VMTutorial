#
# \file texture.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 09-Dec-2023
# \brief Computes eigenvalues of the texture matrix and produces a VTP file for it
#

import json
import os
from copy import copy

import numpy as np
try:
    import shapely
    import shapely.geometry
except ImportError:
    raise Exception('Shapely needs to be installed. \nYou may try installing it with:\n\t pip install shapely')

from ..utils.HalfEdge import *
from .Tensor import Tensor


class MTensor:

    def __init__(self, infile, region=1, pair_sum=False):
        """
          Compute M tensor in Graner statistics (Graner, et al., Eur. Phys. J. E 25, 349–369 (2008))

          Parameters:
          ----------

             infile - JSON file with the current configuration 
             region - number or list (default: 1)
                      if number, averaging is topological
                      region has to be a positive number >= 1. It counts how many next-nearest neighbour shells to include
                      if list, averaging is over a region
                      In this case, there are two options
                        Option 1: list is a list of cells to average over
                        Option 2: the region is list of points that define a polygon (preferably a rectangle) over which the averaging is performed. 
             pair_sum - True or False (default: False)
                      If region > 1 and True, compute average by directly summing over pairs of edges rather than averaging M per cell
                      Otherwise, ignored
        """
        ext = os.path.splitext(infile)[1]
        if ext != '.json' and ext != '.bz2':
            raise Exception('A JSON file with the mesh has to be provided.')
        if type(region) == int:
            topological = True
        elif type(region) == list:
            topological = False
            if any(isinstance(el, list) for el in region):
                poly = shapely.geometry.Polygon(region)
                poly_region = True
            else:
                poly_region = False
        else:
            raise Exception(
                'MTenor: "region" argument has to be either a number or a list.')
        self.mesh = Mesh()
        self.mesh.read(infile)
        if topological:
            self.M = Tensor(self.mesh.num_inner_faces)
            self.Nneigh = np.zeros(self.mesh.num_inner_faces)
        else:
            self.M = Tensor(1)
        self.original_face_id = []
        if topological:             # Topological averaging
            idx = 0
            for f in self.mesh.faces:
                if not (f.outer or f.erased):
                    if "original_id" in f.params:
                        self.original_face_id.append(f.params["original_id"])
                    nv = f.neigh_vectors(self.mesh.box)
                    for v in nv:
                        self.M.T[idx, :, :] += np.array(v.outer(v))
                    # Eq. (A4) in Graner, et al. Eur. Phys. J. E 25, 349–369 (2008)
                    self.M.T[idx, :, :] /= len(nv)
                    self.Nneigh[idx] = len(nv)
                    idx += 1
            if region > 1:           # If region is larger than 1, we do averaging over the number of nearest neighbours set by the region parameter
                idx = 0
                T = Tensor(self.mesh.num_inner_faces)
                for f in self.mesh.faces:
                    if not (f.outer or f.erased):
                        nb = self.mesh.nneighbours(idx, region)
                        Ntot = 0
                        for n in nb:
                            if not pair_sum:
                                ni = self.Nneigh[n]
                                T.T[idx] += ni*self.M.T[n, :, :]
                                Ntot += ni
                            else:
                                nv = self.mesh.faces[n].neigh_vectors(
                                    self.mesh.box)
                                for v in nv:
                                    T.T[idx] += np.array(v.outer(v))
                                    Ntot += 1
                        # Eq. (A5) in Graner, et al. Eur. Phys. J. E 25, 349–369 (2008)
                        T.T[idx] /= Ntot
                        idx += 1
                self.M.T = np.copy(T.T)
        else:
            if poly_region:    # average over a polygonal region
                # Identify all cells with centres inside the interrogation region
                avg_region = []
                for f in self.mesh.faces:
                    if not (f.outer or f.erased):
                        fc = f.rc(self.mesh.box)
                        p = shapely.geometry.Point(fc.r)
                        if poly.contains(p):
                            avg_region.append(f.idx)
                # Loop over all cells inside the interrogation region
                visited = []  # Keep track of those bonds that have been included
                Ntot = 0
                for idx in avg_region:
                    nb = self.mesh.faces[idx].neighbours()
                    nv = self.mesh.faces[idx].neigh_vectors(self.mesh.box)
                    # Find centre of the current face
                    rc1 = self.mesh.faces[idx].rc(self.mesh.box).r
                    for i in range(len(nb)):
                        sp = sorted([idx, nb[i]])
                        if not sp in visited:          # Include only if visited
                            # Find centre of the neighbour
                            rc2 = self.mesh.faces[nb[i]].rc(self.mesh.box).r
                            ln = shapely.geometry.LineString(
                                [rc1, rc2])   # Connect the two with a line
                            v = nv[i]
                            # Weight is the fraction inside the interrogation region
                            w = poly.intersection(ln).length/ln.length
                            self.M.T[0, :, :] += w*np.array(v.outer(v))
                            Ntot += w
                            visited.append(sp)
                self.M.T /= Ntot
            else:   # Average over list of cells
                Ntot = 0
                for idx in region:
                    nv = self.mesh.faces[idx].neigh_vectors(self.mesh.box)
                    for v in nv:
                        self.M.T[0, :, :] += np.array(v.outer(v))
                        Ntot += 1
                self.M.T /= Ntot

    def plot_vtk_tensor(self, filename):
        self.M.plot_vtk_tensor(filename, self.mesh, 'M_tensor')
    
    def plot_vtk_ellipse(self, filename, N = 20, scale = 1.0):
        self.M.plot_vtk_ellipse(filename, self.mesh, N = N, scale = scale)

    def plot_vtk_lines(self, filename, scale = 1.0):
        self.M.plot_vtk_lines(filename, self.mesh, scale = scale)
