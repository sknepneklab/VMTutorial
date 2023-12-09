import numpy as np
from cell import Cell
from vertex import Vertex
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from CellList2D import *
import json
import copy


def area(r):
    i = np.arange(r.shape[0])
    j = np.roll(i, 1)
    return 0.5*np.abs(np.sum(r[i, 0]*r[j, 1]-r[i, 1]*r[j, 0]))


def perim(r):
    i = np.arange(r.shape[0])
    j = np.roll(i, 1)
    dx, dy = r[j, 0] - r[i, 0], r[j, 1] - r[i, 1]
    return np.sum(np.sqrt(dx**2 + dy**2))


class HoneycombLattice:

    def __init__(self, lx, ly, a=1.0):
        if lx <= 0 or ly <= 0:
            raise ValueError("Box has to have positive side length.")
        self.lx = lx
        self.ly = ly
        self.a = a
        self.points = []
        self.cpoints = []
        self.faces = []
        self.vertices = []
        self.cells = []
        self.cell_list = CellList2D([1.5*lx, 1.25*ly], 1.5*a)
        self.boundary = []

    def __build_grid(self):
        vecA, vecB = self.a*np.array([1.0, 0]), self.a*np.array([-1.0, 0.0])
        a1, a2 = self.a*np.array([1.5, 0.5*np.sqrt(3.0)]
                                 ), self.a*np.array([1.5, -0.5*np.sqrt(3.0)])
        l = (a1 + a2)[0]
        N = int(self.lx/l)
        idx = 0
        for i in range(-3*N, 3*N):
            for j in range(-3*N, 3*N):
                x, y = vecA + i*a1 + j*a2
                if -0.6*self.lx <= x <= 0.6*self.lx and -0.5*self.ly <= y <= 0.5*self.ly:
                    self.points.append(np.array([x, y]))
                    self.cell_list.add_particle([x, y], idx)
                    idx += 1
                x, y = vecB + i*a1 + j*a2
                if -0.6*self.lx <= x <= 0.6*self.lx and -0.5*self.ly <= y <= 0.5*self.ly:
                    self.points.append(np.array([x, y]))
                    self.cell_list.add_particle([x, y], idx)
                    idx += 1
                x, y = i*a1 + j*a2
                if -0.6*self.lx <= x <= 0.6*self.lx and -0.5*self.ly <= y <= 0.5*self.ly:
                    self.cpoints.append([x, y])
        self.points = np.asarray(self.points)

    def __build_circle(self, R):
        vecA, vecB = self.a*np.array([1.0, 0]), self.a*np.array([-1.0, 0.0])
        a1, a2 = self.a*np.array([1.5, 0.5*np.sqrt(3.0)]
                                 ), self.a*np.array([1.5, -0.5*np.sqrt(3.0)])
        l = (a1 + a2)[0]
        N = int(2*R/l)
        idx = 0
        for i in range(-5*N, 5*N):
            for j in range(-5*N, 5*N):
                x, y = vecA + i * a1 + j * a2
                r = np.sqrt(x*x + y*y)
                if r <= (R + 0.5*self.a):
                    self.points.append(np.array([x, y]))
                    self.cell_list.add_particle([x, y], idx)
                    idx += 1
                x, y = vecB + i * a1 + j * a2
                r = np.sqrt(x*x + y*y)
                if r <= (R + 0.5*self.a):
                    self.points.append(np.array([x, y]))
                    self.cell_list.add_particle([x, y], idx)
                    idx += 1
                x, y = i * a1 + j * a2
                r = np.sqrt(x*x + y*y)
                if r <= (R + 0.5*self.a):
                    self.cpoints.append([x, y])
        self.points = np.asarray(self.points)

    def __build_polygon(self, region):
        vecA, vecB = self.a*np.array([1.0, 0]), self.a*np.array([-1.0, 0.0])
        a1, a2 = self.a * np.array([1.5, 0.5 * np.sqrt(3.0)]
                                   ), self.a * np.array([1.5, -0.5 * np.sqrt(3.0)])
        l = (a1 + a2)[0]
        poly = Polygon(region)
        xmin, ymin, xmax, ymax = poly.bounds
        L = np.max([xmax - xmin, ymax - ymin])
        N = int(L / l)
        idx = 0
        for i in range(-5*N, 5*N):
            for j in range(-5*N, 5*N):
                x, y = vecA + i * a1 + j * a2
                p = Point([x, y])
                if poly.contains(p):
                    self.points.append(np.array([x, y]))
                    self.cell_list.add_particle([x, y], idx)
                    idx += 1
                x, y = vecB + i * a1 + j * a2
                p = Point([x, y])
                if poly.contains(p):
                    self.points.append(np.array([x, y]))
                    self.cell_list.add_particle([x, y], idx)
                    idx += 1
                x, y = i * a1 + j * a2
                p = Point([x, y])
                if poly.contains(p):
                    self.cpoints.append([x, y])
        self.points = np.asarray(self.points)

    def order_vertices(self, l, rc, clockwise=False):
        angles = []
        for i in l:
            dx, dy = self.vertices[i].r - rc
            angles.append(np.arctan2(dy, dx))
        sidx = np.argsort(angles)
        if clockwise:
            sidx = list(reversed(sidx))
        return l[sidx]

    def order_vertices_boundary(self, l, rc, clockwise=False):
        remain = copy.deepcopy(l)
        i0 = l[0]
        l = [l[0]]  # arbitrary start
        remain = np.setdiff1d(remain, i0)
        while len(remain) > 1:
            v = self.vertices[l[-1]]
            nn = self.nearest_neighbour(remain, v)
            l.append(nn)
            # remove nearest neighbour from remaining list
            remain = np.setdiff1d(remain, nn)
        l.append(remain[0])

        # compute signed area (cf https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order/ for instance ) :
        As = 0
        n = len(l)
        for i in range(n):
            v1, v2 = self.vertices[l[i]], self.vertices[l[(i+1) % n]]
            As += v1.r[0] * v2.r[1] - v1.r[1]*v2.r[0]
        is_clockwise = True if As < 0 else False
        if not is_clockwise and clockwise:
            print("Reversing boundary orientation")
            l.reverse()
        elif is_clockwise and not clockwise:
            print("Reversing boundary orientation")
            l.reverse()

        return np.array(l)

    def nearest_neighbour(self, l, v0):
        n = len(l)
        vs = [self.vertices[i] for i in l]
        dmin = 2*self.lx
        for v in vs:
            d = np.linalg.norm(v.r - v0.r)
            if d < dmin:
                nn = v.id
                dmin = d
        return nn

    def __build_faces(self):
        a2 = self.a**2
        cnt = []
        for c in self.cpoints:
            xc, yc = c
            pts = []
            neighs = self.cell_list.get_neighbours([xc, yc])
            for i in neighs:
                x, y = self.points[i]
                dx, dy = x - xc, y - yc
                dr2 = dx*dx + dy*dy
                if dr2 <= 1.1*a2:
                    pts.append(i)
            if len(pts) == 6:
                self.faces.append(pts)
                cnt.append([xc, yc])
        # remove dangling vertices
        index_map = [-1 for __ in range(self.points.shape[0])]
        vidx = 0
        for i in range(self.points.shape[0]):
            for face in self.faces:
                if i in face:
                    self.vertices.append(Vertex(vidx, self.points[i]))
                    index_map[i] = vidx
                    vidx += 1
                    break
        for fid in range(len(self.faces)):
            c = Cell(fid, cnt[fid])
            for f in self.faces[fid]:
                c.verts.append(index_map[f])
                c.verts = self.order_vertices(np.array(c.verts), c.rc).tolist()

            self.cells.append(c)
        for c in self.cells:
            for i in range(len(c.verts)):
                ip = (i+1) % len(c.verts)
                f, fp = c.verts[i], c.verts[ip]
                if not f in self.vertices[fp].neigh:
                    self.vertices[fp].neigh.append(f)
                if not fp in self.vertices[f].neigh:
                    self.vertices[f].neigh.append(fp)
                if not c.id in self.vertices[c.verts[i]].faces:
                    self.vertices[c.verts[i]].faces.append(c.id)

    def __build_boundary(self, R=None):
        for v in self.vertices:
            if len(v.faces) < 3:
                self.boundary.append(v.id)

        self.boundary = np.array(self.boundary)
        self.boundary = self.order_vertices_boundary(
            self.boundary, np.array([0, 0]), True)

        for b in self.boundary:
            self.vertices[b].boundary = True

    def build(self, R=None, poly=None):
        if R == None and poly == None:
            print("Building grid...")
            self.__build_grid()
        elif R != None:
            print("Building disk...")
            self.__build_circle(R)
        else:
            print('Building polygon...')
            self.__build_polygon(poly)
        print("Building cells...")
        self.__build_faces()
        print("Building boundary...")
        self.__build_boundary(R)

    def set_cell_type(self, region, tp):
        if all(map(lambda x: isinstance(x, list), region)):
            poly = Polygon(region)
            for c in self.cells:
                p = Point(c.rc)
                if poly.contains(p):
                    c.type = tp
        elif all(map(lambda x: isinstance(x, int), region)):
            for c in region:
                self.cells[c].type = tp
        else:
            raise RuntimeError('Unknown region type.')

    def set_vertex_type(self, region, tp, boundary_only=False):
        if len(region) > 2:
            poly = Polygon(region)
            for v in self.vertices:
                can_add = False
                if (not boundary_only) or v.boundary:
                    can_add = True
                p = Point(v.r)
                if can_add and poly.contains(p):
                    v.type = tp
        else:
            if len(region) == 2:
                rmin, rmax = region
                for v in self.vertices:
                    can_add = False
                    if (not boundary_only) or v.boundary:
                        can_add = True
                    x, y = v.r
                    r = np.sqrt(x*x + y*y)
                    if can_add and rmin <= r <= rmax:
                        v.type = tp

    def set_constraint(self, tp, constraint):
        for v in self.vertices:
            if v.type == tp:
                v.constraint = constraint

    def set_l0(self):
        for c in self.cells:
            for i in range(len(c.verts)):
                ip = (i + 1) % len(c.verts)
                x1, y1 = self.vertices[c.verts[i]].r
                x2, y2 = self.vertices[c.verts[ip]].r
                c.l0.append(np.sqrt((x2-x1)**2 + (y2-y1)**2))

    def minmax(self):
        self.minx, self.maxx = self.vertices[0].r[0], self.vertices[0].r[0]
        self.miny, self.maxy = self.vertices[0].r[1], self.vertices[0].r[1]
        for v in self.vertices:
            if v.r[0] < self.minx:
                self.minx = v.r[0]
            if v.r[0] > self.maxx:
                self.maxx = v.r[0]
            if v.r[1] < self.miny:
                self.miny = v.r[1]
            if v.r[1] > self.maxy:
                self.maxy = v.r[1]

    def get_cells_with_vert_type(self, vert_type):
        cells = []
        for c in self.cells:
            for v in c.verts:
                if self.vertices[v].type == vert_type and not c.id in cells:
                    cells.append(c.id)
        return cells

    def json_out(self, fname):
        jsonData = {}
        jsonData["mesh"] = {}
        jsonData["mesh"]["vertices"] = []
        jsonData["mesh"]["box"] = {
            "periodic": False, "lx": self.lx, "ly": self.ly}
        for v in self.vertices:
            vd = {}
            vd["id"] = v.id
            vd["r"] = v.r.tolist()
            vd["type"] = v.type
            vd["erased"] = False
            vd["boundary"] = v.boundary
            vd["coordination"] = len(v.neigh)
            vd["neighbours"] = v.neigh
            vd["constraint"] = v.constraint
            jsonData["mesh"]["vertices"].append(vd)
        jsonData["mesh"]["faces"] = []
        for c in self.cells:
            cd = {}
            cd["id"] = c.id
            cd["outer"] = False
            cd["nsides"] = len(c.verts)
            cd["type"] = c.type
            cd["vertices"] = c.verts
            if len(c.l0) > 0:
                cd["l0"] = c.l0
            r = []
            for v in c.verts:
                r.append(self.vertices[v].r)
            r = np.asarray(r)
            cd["A0"] = area(r)
            cd["P0"] = perim(r)
            jsonData["mesh"]["faces"].append(cd)
        bd = {}
        bd["id"] = self.cells[-1].id + 1
        bd["outer"] = True
        bd["nsides"] = self.boundary.size
        bd["type"] = "passive"
        bd["vertices"] = self.boundary.tolist()
        bd["A0"] = 0.0
        bd["P0"] = 0.0
        jsonData["mesh"]["faces"].append(bd)
        with open(fname, 'w') as out:
            json.dump(jsonData, out, sort_keys=True, indent=4)
