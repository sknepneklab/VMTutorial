from triangular_lattice import *
from make_mesh import *

L = 21

t = TriangularLattice(L)
t.build_periodic_lattice()
m = MakeMesh(t)
m.make_initial_configuration()
m.json_out("hexagonal.json")
