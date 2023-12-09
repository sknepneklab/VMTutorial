from honeycomb_lattice import *

# Genrate hexagonal lattice in a rectangular box [-10,10] x [-10,10]
h = HoneycombLattice(20.0,20.0,1.0)

# build lattice 
h.build()
h.minmax()

# Set a rectangle around the left side of the system
left = [[1.2*h.minx,0.95*h.miny], 
        [0.9*h.minx,0.95*h.miny],
        [0.9*h.minx,0.95*h.maxy],
        [1.2*h.minx,0.95*h.maxy]]

# Set a rectangle around the right side of the system
right = [[0.9*h.maxx,0.95*h.miny], 
         [1.2*h.maxx,0.95*h.miny],
         [1.2*h.maxx,0.95*h.maxy],
         [0.9*h.maxx,0.95*h.maxy]]

top =  [[0.95*h.minx,0.9*h.maxy], 
         [0.95*h.minx,1.2*h.maxy],
         [0.95*h.maxx,1.2*h.maxy],
         [0.95*h.maxx,0.9*h.maxy]]

bottom  =  [[0.95*h.minx,0.9*h.miny], 
         [0.95*h.minx,1.2*h.miny],
         [0.95*h.maxx,1.2*h.miny],
         [0.95*h.maxx,0.9*h.miny]]



# Set boundary vertices on the left as type 1
h.set_vertex_type(left, "left", True)
# Set boundary vertices on the right as type 2
h.set_vertex_type(right, "right", True)

# Print everything out as JSON
h.json_out("honeycomb.json")
