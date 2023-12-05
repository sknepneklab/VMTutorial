#
# \file Face.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 04-Dec-2024
# Handles faces in the half-edge implementation of the mesh

from .Vec import Vec
from .Tensor import *

class Face:

    def __init__(self, idx):
        self.idx = idx
        self.he = None
        self.outer = False
        self.type = None
        self.erased = False
        self.n = None
        self.A0 = None
        self.property = None
        self.params = {}

    def area(self):
        # Use shoe-lace formula to compute area
        # By defintion, outer face will have negative area
        first = self.he
        he = self.he
        area = 0.0
        r0 = he.vfrom.r
        while True:
            vi = he.vfrom
            vj = he.vto
            ri = vi.r - r0
            rj = vj.r - r0
            area += ri.r[0]*rj.r[1] - ri.r[1]*rj.r[0]
            he = he.next
            if he.idx == first.idx:
                break
        return 0.5*abs(area)

    def perimeter(self):
        first = self.he
        he = self.he
        perim = 0.0
        while True:
            ri = he.vfrom.r
            rj = he.vto.r
            dr = rj - ri
            perim += dr.length()
            he = he.next
            if he.idx == first.idx:
                break
        return perim

    def rc(self, box=None):
        first = self.he
        he = self.he
        Rc = Vec([0.0, 0.0])
        N = 0
        while True:
            dr = he.vfrom.r - first.vfrom.r
            Rc += Vec(dr.r)
            N += 1
            he = he.next
            if he.idx == first.idx:
                break
        Rc = (1.0/N)*Rc + first.vfrom.r
        return Vec(Rc.r, box)

    def neighbours(self):
        neigh = []
        first = self.he
        he = self.he
        while True:
            if not he.pair.face.outer:
                neigh.append(he.pair.face.idx)
            he = he.next
            if he.idx == first.idx:
                break
        return neigh

    def neigh_vectors(self, box=None):
        nv = []
        first = self.he
        he = self.he
        Rc = self.rc(box)
        while True:
            if not he.pair.face.outer:
                Rci = he.pair.face.rc(box)
                dR = Rci - Rc
                nv.append(dR)
            he = he.next
            if he.idx == first.idx:
                break
        return nv

    def edge_vectors(self):
        ev = []
        first = self.he
        he = self.he
        while True:
            dR = (he.vto.r - he.vfrom.r).unit()
            ev.append(dR)
            he = he.next
            if he.idx == first.idx:
                break
        return ev

    def num_sides(self):
        first = self.he
        he = self.he
        N = 0
        while True:
            N += 1
            he = he.next
            if he.idx == first.idx:
                break
        return N

    def centroid(self):
        first = self.he
        he = self.he
        r0 = he.vfrom.r
        Rc = Vec([0.0, 0.0])
        while True:
            ri = he.vfrom.r - r0
            rj = he.vto.r - r0
            fact = ri.r[0]*rj.r[1] - ri.r[1]*rj.r[0]
            Rc.r[0] += (ri.r[0] + rj.r[0])*fact
            Rc.r[1] += (ri.r[1] + rj.r[1])*fact
            he = he.next
            if he.idx == first.idx:
                break
        return (1.0/(6*self.area()))*Rc + r0
    
    def get_shape(self):
        first = self.he
        he = self.he
        G = Tensor([0,0,0,0])
        while True:
            r = he.vto.r - he.vfrom.r
            T = Tensor(r).traceless()
            G += T
            he = he.next
            if he.idx == first.idx:
                break
        return (1.0/self.area())*G 

    def get_nematic(self):
        return self.get_shape().principal_evec()

    def crosses_boundary(self):
        first = self.he
        he = self.he
        if first.vfrom.r.box is None:
            return False
        box = first.vfrom.r.box
        l = box.h.diagonal()

        while True:
            dr = he.vto.r.r - he.vfrom.r.r
            if np.any(np.abs(dr) > 0.5*l):
                return True 
            he = he.next
            if he.idx == first.idx:
                break   
        return False