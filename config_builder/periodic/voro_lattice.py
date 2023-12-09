import numpy as np
from scipy.spatial import Voronoi
from CellList2D import CellList2D
from vertex import Vertex


class Box:

    def __init__(self, lx=10.0, ly=None):
        self.lx = lx
        if ly == None:
            self.ly = self.lx
        else:
            self.ly = ly


class VoroLattice:

    def __init__(self, N, box, min_dist=1.0, seed=None):
        self.N = N
        self.box = box
        self.Lx = box.lx
        self.Ly = box.ly
        self.min_dist = 1.0
        self.cl = CellList2D(box.lx, box.ly, 2.5*min_dist)
        self.seed = seed
        self.max_attempts = 5*N
        self.is_periodic = False
        np.random.seed(seed)

    def __apply_periodic(self, r):
        lx, ly = self.box.lx, self.box.ly
        x, y = r
        x = x + lx if x < -0.5 * lx else x
        x = x - lx if x > 0.5 * lx else x
        y = y + ly if y < -0.5 * ly else y
        y = y - ly if y > 0.5 * ly else y
        return np.array([x, y])

    def __build_lattice(self):
        i = 0
        attempt = 0
        self.r = np.zeros((self.N, 2))
        while i < self.N:
            x, y = np.random.uniform(-0.5*self.box.lx, 0.5 *
                                     self.box.lx), np.random.uniform(-0.5*self.box.ly, 0.5*self.box.ly)
            r = np.array([x, y])
            neighs = self.cl.get_neighbours(r)
            can_add = True
            for n in neighs:
                rn = self.r[n, :]
                dr = rn - r
                dr = self.__apply_periodic(dr)
                if np.sqrt(np.dot(dr, dr)) < self.min_dist:
                    can_add = False
                    break
            if can_add:
                self.r[i, :] = r
                self.cl.add_particle(r, i)
                i += 1
            if attempt > self.max_attempts:
                raise Exception('Failed to place points after '+str(attempt) +
                                ' attempts. Please reduce the cutoff distance.')
            else:
                attempt += 1
        self.built = True
        self.r_original = np.copy(self.r)

    def __pad(self, f=0.5):
        ax = f*self.box.lx
        ay = f*self.box.ly
        xmin, xmax = -0.5*self.box.lx, 0.5*self.box.lx
        ymin, ymax = -0.5*self.box.ly, 0.5*self.box.ly
        self.rp = np.copy(self.r)
        rl = self.r[self.r[:, 0] - xmin <= ax] + np.array([self.box.lx, 0])
        rr = self.r[xmax - self.r[:, 0] <= ax] - np.array([self.box.lx, 0])
        rd = self.r[self.r[:, 1] - ymin <= ay] + np.array([0, self.box.ly])
        ru = self.r[ymax - self.r[:, 1] <= ay] - np.array([0, self.box.ly])
        rc1 = self.r[(self.r[:, 0] - xmin <= ax) & (self.r[:, 1] -
                                                    ymin <= ay)] + np.array([self.box.lx, self.box.ly])
        rc2 = self.r[(self.r[:, 0] - xmin <= ax) & (ymax - self.r[:, 1]
                                                    <= ay)] + np.array([self.box.lx, -self.box.ly])
        rc3 = self.r[(xmax - self.r[:, 0] <= ax) & (self.r[:, 1] -
                                                    ymin <= ay)] + np.array([-self.box.lx, self.box.ly])
        rc4 = self.r[(xmax - self.r[:, 0] <= ax) & (ymax - self.r[:, 1]
                                                    <= ay)] + np.array([-self.box.lx, -self.box.ly])
        self.Npad = rl.shape[0] + rr.shape[0] + rd.shape[0] + ru.shape[0]
        self.rp = np.vstack((self.rp, rl, rr, rd, ru, rc1, rc2, rc3, rc4))

    def __vorobuild(self, f=0.5):
        self.__pad(f)
        self.vor = Voronoi(self.rp)

    def __centroid(seld, verts):
        x = verts[:, 0]
        y = verts[:, 1]
        a = x*np.roll(y, 1) - np.roll(x, 1)*y
        area = 0.5*np.sum(a)
        cx = np.sum((x + np.roll(x, 1))*a)/(6.0*area)
        cy = np.sum((y + np.roll(y, 1))*a)/(6.0*area)
        return [cx, cy]

    def __relax(self, f=0.1, eps=1e-3, max_iter=1000):
        new_r = np.zeros_like(self.r)
        dlmin = 1e10
        iteration = 0
        while iteration < max_iter:
            self.__vorobuild(f)
            i = 0
            for r in self.vor.point_region[:self.N]:
                region = self.vor.regions[r]
                new_r[i, :] = self.__apply_periodic(
                    self.__centroid(self.vor.vertices[region, :]))
                i += 1
            dl = []
            for i in range(self.N):
                dr = self.__apply_periodic(new_r[i, :] - self.r[i, :])
                dl.append(np.sqrt(dr[0]**2 + dr[1]**2))
            smallest_dl = np.max(dl)
            if smallest_dl < dlmin:
                best_r = np.copy(new_r)
                dlmin = smallest_dl
            if all(map(lambda x: x < eps, dl)):
                break
            else:
                iteration += 1
                print('Iteration {:d} with max eps value {:.6f}.'.format(
                    iteration, np.max(dl)))
            self.r = np.copy(new_r)
        else:
            self.r = np.copy(best_r)
            print('Failed to converge after {:d} iterations.'.format(max_iter))

    def __copy_periodic(self):
        ix = [-1, 0, 1]
        iy = [-1, 0, 1]
        copy = []
        idx = len(self.points)
        for iix in ix:
            for iiy in iy:
                if iix != 0 or iiy != 0:
                    for p in self.points:
                        r = p.r + np.array([iix*self.box.lx, iiy*self.box.ly])
                        copy.append(Vertex(idx, r))
                        copy[-1].original_idx = p.id
                        idx += 1

        self.original_points = np.copy(self.points)
        self.points = np.append(self.points, np.array(copy), axis=0)
        self.is_periodic = True

    def build_lattice(self, f=0.5, eps=1e-3, max_iter=1000):
        self.__build_lattice()
        self.points = []
        if eps > 0:
            self.__relax(f, eps, max_iter)
        for i in range(self.N):
            self.points.append(Vertex(i, self.r[i]))
        self.__copy_periodic()
