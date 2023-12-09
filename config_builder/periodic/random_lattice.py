import numpy as np
from vertex import Vertex
from CellList2D import CellList2D


class RandomLattice:

    def __init__(self, N, Lx, Ly=None, min_dist=1.0, seed=None):
        self.N = N
        self.Lx = Lx
        self.Ly = Ly if Ly != None else Lx
        self.min_dist = min_dist
        self.cl = CellList2D(Lx, Ly, 2.5*min_dist)
        self.built = False
        self.is_periodic = False
        self.max_attempts = 5*N
        self.seed = seed
        np.random.seed(self.seed)

    def apply_periodic(self, r):
        Lx, Ly = self.Lx, self.Ly
        x, y = r
        x = x + Lx if x < -0.5 * Lx else x
        x = x - Lx if x > 0.5 * Lx else x
        y = y + Ly if y < -0.5 * Ly else y
        y = y - Ly if y > 0.5 * Ly else y
        return np.array([x, y])

    def build_lattice(self, points=None):
        if points is None:
            i = 0
            points = []
            attempt = 0
            while i < self.N:
                x, y = np.random.uniform(-0.5*self.Lx, 0.5 *
                                         self.Lx), np.random.uniform(-0.5*self.Ly, 0.5*self.Ly)
                r = np.array([x, y])
                neighs = self.cl.get_neighbours(r)
                can_add = True
                for n in neighs:
                    vn = points[n]
                    rn = vn.r
                    dr = rn - r
                    dr = self.apply_periodic(dr)
                    if np.sqrt(np.dot(dr, dr)) < self.min_dist:
                        can_add = False
                        break
                if can_add:
                    points.append(Vertex(i, r))
                    self.cl.add_particle(r, i)
                    i += 1
                if attempt > self.max_attempts:
                    raise Exception('Failed to place points after '+str(attempt) +
                                    ' attempts. Please reduce the cutoff distance.')
                else:
                    attempt += 1
            self.built = True
            self.points = np.array(points)
        else:
            pts = []
            for i in range(points.shape[0]):
                pts.append(Vertex(i, points[i, :]))
            self.points = np.array(pts)

    def copy_periodic(self):
        ix = [-1, 0, 1]
        iy = [-1, 0, 1]
        copy = []
        idx = len(self.points)
        for iix in ix:
            for iiy in iy:
                if iix != 0 or iiy != 0:
                    for p in self.points:
                        r = p.r + np.array([iix*self.Lx, iiy*self.Ly])
                        copy.append(Vertex(idx, r))
                        copy[-1].original_idx = p.id
                        idx += 1

        self.original_points = np.copy(self.points)
        self.points = np.append(self.points, np.array(copy), axis=0)
        self.is_periodic = True

    def build_periodic_lattice(self, points=None):
        self.build_lattice(points)
        self.copy_periodic()
