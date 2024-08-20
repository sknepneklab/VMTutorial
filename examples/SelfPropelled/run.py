
# ###############################################################
#
#  Performs simulations of a single self-propelled cell moving 
#  in a tissue
#
#  Author: Rastko Sknepnek, (c) 2024
#  Date: 20-Aug-2024
#
# ###############################################################


# ###############################################################
#
# Load standard Python modules
#
# ###############################################################

import os
import argparse
import numpy as np

# ###############################################################
#
# Load VMTutorial modules
#
# ###############################################################
from VMToolkit.VM import *

# ###############################################################
#
# Read command line arguments 
#
# ###############################################################

argparser = argparse.ArgumentParser()
argparser.add_argument("--input", type=str, default="random_normal_sp.json", help="Input file with the initial configuration")
argparser.add_argument("--dt", type=float, default=0.05, help="Time step for the simulation")
argparser.add_argument("--p0", type=float, default=3.85, help="shape parameter p0 = P0/sqrt(A0)")
argparser.add_argument("--v0", type=float, default=0.1, help="self-propulsion velocity")
argparser.add_argument("--nrun", type=int, default=500, help="number of steps to run")
argparser.add_argument("--nrelax", type=int, default=50, help="number of relaxation steps to run")
argparser.add_argument("--dumpfreq", type=int, default=1, help="how often to dump data")
argparser.add_argument("--cellout", type=str, default="cells", help="output file name prefix for cells")
argparser.add_argument("--dirout", type=str, default="dir", help="output file name prefix for directors")
argparser.add_argument("--jsonout", type=str, default="cells", help="output file name prefix for JSON files")
argparser.add_argument("--seed", type=int, default=0, help="random seed")
argparser.add_argument("--Dr", type=float, default=0.1, help="rotation diffusion constant")
args = argparser.parse_args()

# ###############################################################
#
# Specify random number generator seed
#
# ###############################################################

if args.seed < 0:
    seed = int(os.urandom(4).hex(), base=16) % 2**30  # generate some large number as seed
else:
    seed = args.seed

# ################################################################
#
# Set up simulation objects
#
# ################################################################
tissue  = Tissue()                                               # initialise mesh
sim_sys = System(tissue)                                         # base object for the system
forces = Force(sim_sys)                                          # handles all types of forces
integrators = Integrate(sim_sys, forces, args.seed)              # handles all integrators
topology = Topology(sim_sys, forces)                             # handles all topology changes (T1, division, ingression)
dumps = Dump(sim_sys, forces)                                    # handles all data output 
simulation = Simulation(sim_sys, integrators, forces, topology)  # simulation object

# ################################################################
#
# Read the initial configuration
#
# ################################################################
sim_sys.read_input(args.input)

# ###############################################################
#
# Calculate average cell area and use it as A0
#
# ###############################################################
A0 = 0.0
for c in tissue.cells():
    A0 += tissue.cell_area(c.id)
A0 /= len(tissue.cells())

# ###############################################################
#
# Set simulation parameters
#
# ###############################################################
kappa = 1.0
gamma = 1.0
lam = args.p0 * np.sqrt(A0) * gamma

# ---------- Forces ----------------------
forces.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
forces.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 - lambda*P 

#----------- Force parameters for both cell types -----------------------
forces.set_params('area', 'passive', {'kappa' : kappa})
forces.set_params('area', 'active', {'kappa' : kappa})

forces.set_params('perimeter', 'passive',  {'gamma': gamma, 'lambda': lam})   
forces.set_params('perimeter', 'active',  {'gamma': gamma, 'lambda': lam})   

# ##################################################################
#
# unit of time
#
# #################################################################
eta = 1.0 # friction coefficient
tau = eta/(kappa*A0) # characteristic time scale

# ##################################################################
#
# Set time step to dt*tau
#
# #################################################################
dt = args.dt*tau

# ##################################################################
#
# Define the number of steps to run
#
# #################################################################
nsteps = int(round(0.5/dt))  # number of steps to run is set to 0.5*tau. Please change if you need more frequent data output

# #################################################################
#
# Set conditions for the T1 transition
#
# #################################################################
topology.set_params({'min_edge_len': 0.01, 'new_edge_len': 0.011}) # T1 transition parameters. Please ignore myosin.


# ##################################################################
#
# Define integrators
#
# #################################################################
integrators.add('brownian')         # add Brownian integrator that handles all vertex movements


# #################################################################
#
# Simulation starts here
#
# #################################################################
integrators.set_dt(dt)       # all integrators have the same time step size

# ##############################################################
#
# We first do the relaxation part
#
# ##############################################################

print('Running relaxation...')
for i in range(0, args.nrelax):
   if i % args.dumpfreq == 0:
      dumps.dump_cells(f'relax_{args.cellout}_p0_{args.p0:.3f}_v0_{args.v0:.2f}_{i:06d}.vtp', draw_periodic=True)
   simulation.run(nsteps)

dumps.dump_mesh(f'relaxed_{args.jsonout}_p0_{args.p0:.3f}_v0_{args.v0:.2f}.json')

print('Relaxation done...')

# ################################################################
#
# Now we proceed to the actual simulation with self-propulsion
#
# ################################################################

# We need to add self-propulsion force to the system and apply it to active cells only

#----------- Forces -------------------------
forces.add('self-propulsion') # add self-propulsion force term

#----------- Force parameters --------------------------------------
forces.set_params('self-propulsion', 'active', {'v0': args.v0})
forces.set_vec_params('self-propulsion', 'active', {'n': Vec(1.0, 0.0)})  # set the direction of self-propulsion


# ################################################################
#
# In case we want to update the direction of self-propulsion we 
# need to uncomment the following two lines
#
# ################################################################
#integrators.set_flag('brownian', 'update_n')
#integrators.set_params('brownian', {'Dr': args.Dr}) # set the parameters for the brownian integrator


# Run the simulation for nrun steps 
print('Running self-propulsion simulation...')
for i in range(0, args.nrun):
   if i % args.dumpfreq == 0:
      dumps.dump_cells(f'{args.cellout}_p0_{args.p0:.3f}_v0_{args.v0:.2f}_{i:06d}.vtp', draw_periodic=True)
      dumps.dump_cell_directors(f'dir_{args.cellout}_p0_{args.p0:.3f}_v0_{args.v0:.2f}_{i:06d}.vtp', draw_periodic=True)
   simulation.run(nsteps)


