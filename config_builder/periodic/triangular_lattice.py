import numpy as np
from vertex import Vertex


class TriangularLattice:

    def __init__(self, Lx, Ly=None, a=1):
        self.a = a
        self.Lx = Lx
        self.Ly = Ly if Ly != None else Lx
        self.a1 = np.array([1, 0])
        self.a2 = np.array([0.5, 0.5 * np.sqrt(3)])
        self.built = False
        self.is_periodic = False

    def __inside(self, r):
        x, y = r
        return (-0.5*self.Lx <= x < 0.5*self.Lx and -0.5*self.Ly <= y < 0.5*self.Ly)

    def build_lattice(self):
        points = []
        Nx = int(self.Lx / self.a + 1)
        Ny = int(self.Ly / self.a + 1)
        if Ny % 2 != 0:
            print(
                'Error! An odd number (Ny = {:d}) rows of points.'.format(Ny))
            raise Exception(
                'Ly has to be such that there is an even number of row of points.')
        self.Ly = 0.5*Ny*np.sqrt(3)
        idx = 0
        for i in range(-3*Nx // 2, 3*Nx // 2 + 1):
            for j in range(-3*Ny // 2, 3*Ny // 2 + 1):
                # - 0.25*self.a*np.array([1,0]) - 0.5*self.a*np.array([0.5,0.5*np.sqrt(3)])
                r = i * self.a * self.a1 + j * self.a * self.a2
                if self.__inside(r):
                    points.append(Vertex(idx, r))
                    idx += 1
        self.points = np.array(points)
        self.built = True

    def copy_periodic(self):
        # for i in range(self.points.shape[0]):
        #   self.points[i,:] -= np.array([0.0,0.25*self.a*np.sqrt(3)])
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

    def build_periodic_lattice(self):
        self.build_lattice()
        self.copy_periodic()
