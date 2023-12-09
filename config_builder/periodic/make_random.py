from random_lattice import *
from make_mesh import *
import os 

seed = int(os.urandom(4).hex(), base=16)  # generate some large number as seed, it will we stored in the JSON file

# Parameters:
# 1st : number of particles
# 2nd : x-direction box size
# 3rd : y-direction box size
# 4th : minimum distance between particles
# 5th : random number generator seed
t = RandomLattice(400, 20, 20, 0.6204032394013997, seed = seed)   

# Buld lattice of points
t.build_periodic_lattice()

# Create Voronoi mesh
m = MakeMesh(t)

# Make actual mesh
m.make_initial_configuration()

# Save it into ouput file
m.json_out('random_conf.json')
