import numpy as np
from cell import Cell
from vertex import Vertex
from scipy.spatial import Delaunay
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from CellList2D import *
from itertools import permutations
import json
import copy


class MakeMesh:

    def __init__(self, lattice):
        self.lattice = lattice
        self.points = np.zeros((len(self.lattice.points), 2))
        for i in range(len(self.lattice.points)):
            self.points[i, :] = self.lattice.points[i].r
        self.has_triangulation = False
        self.periodic_neighbours_matched = False

    def polygon_area(self, verts, points):
        r0 = points[verts[0], :]
        dr = np.zeros((verts.size, 2))
        for i in range(verts.size):
            ri = points[verts[i], :]
            dr[i, :] = self.apply_periodic(ri - r0)
        i = np.arange(verts.size)
        j = np.roll(i, -1)
        return 0.5 * np.sum(dr[i, 0] * dr[j, 1] - dr[j, 0] * dr[i, 1])

    def polygon_perim(self, verts, points):
        r0 = points[verts[0], :]
        dr = np.zeros((verts.size, 2))
        for i in range(verts.size):
            ri = points[verts[i], :]
            dr[i, :] = self.apply_periodic(ri - r0)
        return np.sum(np.sqrt((dr[1:, 0]-dr[:-1, 0])**2 + (dr[1:, 1]-dr[:-1, 1])**2)) + np.sqrt((dr[-1, 0]-dr[0, 0])**2 + (dr[-1, 1]-dr[0, 1])**2)

    def polygon_centre(self, verts, points):
        r0 = points[verts[0], :]
        dr = np.zeros((verts.size, 2))
        for i in range(verts.size):
            ri = points[verts[i], :]
            dr[i, :] = self.apply_periodic(ri - r0)
        return r0 + np.mean(dr, axis=0)

    def compute_l0(self, verts, points):
        r0 = points[verts[0], :]
        dr = np.zeros((verts.size, 2))
        for i in range(verts.size):
            ri = points[verts[i], :]
            dr[i, :] = self.apply_periodic(ri - r0)
        l0 = np.sqrt((dr[1:, 0]-dr[:-1, 0])**2 +
                     (dr[1:, 1]-dr[:-1, 1])**2).tolist()
        l0.append(np.sqrt((dr[-1, 0]-dr[0, 0])**2 + (dr[-1, 1]-dr[0, 1])**2))
        return l0

    def apply_periodic(self, r):
        Lx, Ly = self.lattice.Lx, self.lattice.Ly
        x, y = r
        x = x + Lx if x < -0.5 * Lx else x
        x = x - Lx if x > 0.5 * Lx else x
        y = y + Ly if y < -0.5 * Ly else y
        y = y - Ly if y > 0.5 * Ly else y
        return np.array([x, y])

    def triangle_circumcentre(self, verts):
        if verts.size != 3:
            raise Exception('A triangle is required to find circumcentre.')
        Lx, Ly = self.lattice.Lx, self.lattice.Ly
        r0 = self.points[verts[0]]
        r1 = self.points[verts[1]]
        r2 = self.points[verts[2]]
        r10 = self.apply_periodic(r1 - r0)
        r20 = self.apply_periodic(r2 - r0)
        Bx, By = r10
        Cx, Cy = r20
        D = 2 * (Bx * Cy - By * Cx)
        if D == 0:
            raise Exception('Encountered three collinear points.')
        Ux = (Cy * (Bx ** 2 + By ** 2) - By * (Cx ** 2 + Cy ** 2)) / D + r0[0]
        Uy = (Bx * (Cx ** 2 + Cy ** 2) - Cx * (Bx ** 2 + By ** 2)) / D + r0[1]
        return self.apply_periodic(np.array([Ux, Uy]))

    def build_triangulation(self):
        if not self.lattice.built:
            raise Exception(
                'Primal lattice needs to be built first. Please use build_lattice() function.')
        self.triang = Delaunay(self.points)
        self.inner_idx = []
        if self.lattice.is_periodic:
            for p in self.lattice.original_points:
                self.inner_idx.append(p.id)
        else:
            for p in self.lattice.points:
                self.inner_idx.append(p.id)
        self.has_triangulation = True

    def match_periodic_neighbours(self):
        if not self.has_triangulation:
            raise Exception(
                'Triangulation needs to be built before one can match neighours.')
        if not self.lattice.is_periodic:
            raise Exception(
                'Primal lattice has to be periodic in order to find its periodic neighbours.')

        self.simplices = []
        for s in self.triang.simplices:
            if any(map(lambda i: i in self.inner_idx, s)):
                simplex = []
                for v in s:
                    if v in self.inner_idx:
                        simplex.append(v)
                    else:
                        simplex.append(self.lattice.points[v].original_idx)
                if not any(map(lambda x: list(x) in self.simplices, permutations(simplex))):
                    self.simplices.append(simplex)
        self.simplices = np.array(self.simplices)
        for f in range(len(self.simplices)):
            for p in self.lattice.original_points[self.simplices[f]]:
                if not f in p.faces:
                    p.faces.append(f)
        for p in self.lattice.original_points:
            r0 = p.r
            angles = np.zeros(len(p.faces))
            i = 0
            for f in p.faces:
                rn = self.triangle_circumcentre(self.simplices[f])
                p.face_centres.append(rn)
                dr = self.apply_periodic(rn - r0)
                angles[i] = np.arctan2(dr[1], dr[0])
                i += 1
            p.faces = np.array(p.faces)[np.argsort(angles)]
            p.face_centres = np.array(p.face_centres)[np.argsort(angles), :]
        self.periodic_neighbours_matched = True

    def find_neighbours(self):
        if not self.has_triangulation:
            raise Exception('Triangulation has not been generated.')
        for i in range(len(self.triang.simplices)):
            verts = self.triang.simplices[i]
            for s in self.triang.neighbors[i]:
                nverts = self.triang.simplices[s]
                intsec = set(verts).intersection(set(nverts))
                if len(intsec) == 2:
                    vi, vj = intsec
                    if vi in self.inner_idx or vj in self.inner_idx:
                        if not vi in self.inner_idx:
                            vi = self.lattice.points[vi].original_idx
                        if not vj in self.inner_idx:
                            vj = self.lattice.points[vj].original_idx
                        if not vi in self.lattice.points[vj].neigh:
                            self.lattice.points[vj].neigh.append(vi)
                        if not vj in self.lattice.points[vi].neigh:
                            self.lattice.points[vi].neigh.append(vj)
        if self.lattice.is_periodic:
            points = self.lattice.original_points
        else:
            points = self.lattice.points
        for p in points:
            r0 = p.r
            angles = np.zeros(len(p.neigh))
            i = 0
            for pn in points[p.neigh]:
                rn = pn.r
                dr = self.apply_periodic(rn - r0)
                angles[i] = np.arctan2(dr[1], dr[0])
                i += 1
            p.neigh = np.array(p.neigh)[np.argsort(angles)]
            angles = np.sort(angles)
            da = np.abs(angles - angles[np.roll(np.arange(angles.size), -1)])
            da = np.where(da > np.pi, da - 2*np.pi, da)
            if np.sum(np.abs(da)) < 2 * np.pi:
                p.boundary = True

    def make_mesh(self):
        if not self.periodic_neighbours_matched:
            raise Exception(
                'One needs to run match_periodic_neighbours() before the mesh can be made.')
        self.mesh_vertices = []
        self.mesh_faces = []
        for simplex in self.simplices:
            c = self.triangle_circumcentre(simplex)
            self.mesh_vertices.append(c)
        for p in self.lattice.original_points:
            self.mesh_faces.append(p.faces)
        self.mesh_vertices = np.array(self.mesh_vertices)
        #self.mesh_faces = np.array(self.mesh_faces)

    def json_out(self, fname, params=None, stress_free=False):
        jsonData = {}
        jsonData["mesh"] = {}
        jsonData["mesh"]["vertices"] = []
        if hasattr(self.lattice, 'seed'):
            jsonData["mesh"]["seed"] = self.lattice.seed
        if not stress_free:
            jsonData["mesh"]["l0"] = 1.0
        jsonData["mesh"]["box"] = {"periodic": True,
                                   "lx": self.lattice.Lx, "ly": self.lattice.Ly}
        for i in range(self.mesh_vertices.shape[0]):
            vd = {}
            vd["id"] = i
            vd["r"] = self.mesh_vertices[i].tolist()
            vd["type"] = "regular"
            vd["erased"] = False
            vd["boundary"] = False
            vd["constraint"] = 'none'
            jsonData["mesh"]["vertices"].append(vd)
        jsonData["mesh"]["faces"] = []
        for c in range(len(self.mesh_faces)):  # .shape[0]):
            cd = {}
            cd["id"] = c
            cd["outer"] = False
            cd["nsides"] = len(self.mesh_faces[c])
            cd["type"] = 'passive'
            cd["vertices"] = self.mesh_faces[c].tolist()
            cd["P0"] = self.polygon_perim(
                self.mesh_faces[c], self.mesh_vertices)
            if params != None:
                if type(params) == dict:
                    if "A0" in params:
                        cd["A0"] = params['A0']
                    else:
                        cd["A0"] = self.polygon_area(
                            self.mesh_faces[c], self.mesh_vertices)
                    if "kappa" in params:
                        cd["kappa"] = params['kappa']
                    if 'gamma' in params:
                        cd["gamma"] = params['gamma']
                    if 'lambda' in params:
                        cd['lambda'] = params['lambda']
                    if stress_free:
                        if 'gamma' in cd:
                            cd['lambda'] = cd["gamma"]*cd["P0"]
                        else:
                            raise(
                                'Gamma needs to be set before we can find lambda for the stress free configuration.')
                        cd['l0'] = self.compute_l0(
                            self.mesh_faces[c], self.mesh_vertices)
                    if 'maxA0' in params:
                        cd['maxA0'] = params['maxA0']
                else:
                    raise Exception(
                        'params has to be a dictionary with parameters')
            else:
                cd["A0"] = self.polygon_area(
                    self.mesh_faces[c], self.mesh_vertices)
            jsonData["mesh"]["faces"].append(cd)
        with open(fname, 'w') as out:
            json.dump(jsonData, out, sort_keys=True, indent=4)

    def make_initial_configuration(self):
        self.build_triangulation()
        self.match_periodic_neighbours()
        self.make_mesh()
