# ###############################################################
#
#  Performs simulations of a random patch of particles
#
#  Author: Rastko Sknepnek, (c) 2023
#  Date: 03-Dec-2023
#
# ###############################################################


# ###############################################################
#
# Load standard Python modules
#
# ###############################################################
import sys as s
import argparse 
import numpy as np

# ###############################################################
#
# Load VMYutorial modules
#
# ###############################################################
from vmtutorial.vmtutorial import *

# ###############################################################
#
# Read command line arguments 
#
# ###############################################################
parser = argparse.ArgumentParser()
parser.add_argument('--input', dest = 'input', type = str, default = 'random_conf.json' , help = 'input fule')
parser.add_argument('--lambda', dest = 'lam', type = float, default = 1.7 , help = 'lambda for the vertex model')
parser.add_argument('--dt', dest = 'dt', type = float, default = 0.1, help = 'timestep')
parser.add_argument('--seed', dest = 'seed', type = int, default = None, help = 'random number generator seed')
parser.add_argument('--dumpfreq', dest = 'dumpfreq', type = int, default = 100, help = 'how often to produce output')
parser.add_argument('--nrun', dest = 'nrun', type = int, default = 100, help = 'number of run steps')
args = parser.parse_args()

# ###############################################################
#
# Set parameters
#
# ###############################################################
freq = int(round(1.0/args.dt))   # This makes sure that we output data once per unit of time

# Cell mechanics parameters
kappa = 1.0     # area stiffness
gamma = 0.25     # perimeter stiffness
lam = args.lam


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

# #################################################################
#
# Create the initial configuration and read it
#
# #################################################################

sim_sys.read_input(args.input)           # read input configuration


# #################################################################
#
# Add forces to the system
#
# #################################################################

forces.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
forces.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)

# Set parameters for each cell type
forces.set_params('area', 'passive', {'kappa' : kappa})
forces.set_params('perimeter', 'passive',  {'gamma': gamma, 'lambda': lam})


# #################################################################
#
# Set conditions for the T1 transition
#
# #################################################################

topology.set_params({'min_edge_len': 0.05, 'new_edge_len': 0.055}) 


# #################################################################
#
# Add Brownian integrator that will handle mechanical part
#
# #################################################################

integrators.add('brownian')    


# #################################################################
#
# Simulation starts here
#
# #################################################################

integrators.set_dt(args.dt) # set time step

step = 0       # Step counter in terms of time units

# Pulling on the passive system
for i in range(args.nrun):
    if i % args.dumpfreq == 0:
        dumps.dump_cells(f'cells_{i:08d}.vtp', draw_periodic=True)
        dumps.dump_mesh(f'mesh_{i:08d}.json')
    simulation.run(int(round(freq)))
    step += 1 















