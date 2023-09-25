from __future__ import division, print_function
from sys import stdout, argv
# OpenMM Imports
from openmm.unit import *
from openmm.app import *
from openmm import *
import numpy as np


pdb_file = argv[1]
pdb = PDBFile(pdb_file)

modeller = Modeller(pdb.topology, pdb.positions)
ff = ForceField('charmm36.xml')
modeller.addHydrogens(ff, pH=7.0)

system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('all.pdb', 'w'))
print('Done!')
