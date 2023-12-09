#from random_lattice import *
from voro_lattice import *
from make_mesh import *
import numpy as np
import os

N = 1000 
Lx, Ly = 60, 60
a = 0.5*np.sqrt(Lx*Ly/N)
eps = -1
max_iter = 5000

seed = int(os.urandom(4).hex(), base=16)  # generate some large number as seed, it will we stored in the JSON file

t = VoroLattice(N, Box(Lx,Ly), a, seed = seed)   # change second and third arguments to make a larger system
print('Generating Voronoi lattice...')
t.build_lattice(0.5, eps, max_iter)
print('Building mesh...')
m = MakeMesh(t)
m.make_initial_configuration()
m.json_out('random_test.json', {'kappa': 1.0, 'gamma': 0.25, 'A0': Lx*Ly/N}, stress_free = True)
