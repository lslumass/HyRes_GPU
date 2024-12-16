from __future__ import division, print_function
import argparse
from HyresBuilder import HyresFF
# OpenMM Imports
from openmm.unit import *
from openmm.app import *
from openmm import *
import numpy as np

# 0) set variables in the simulation
gpu_id = "0"
top_inp = './hyres_forcefield/param_hyres_GPU.inp'
param_inp = './hyres_forcefield/top_hyres_GPU.inp'

# input parameters
parser = argparse.ArgumentParser()
parser.add_argument('-c', "--pdb", default='conf.pdb', help="pdb file")
parser.add_argument('-p', "--psf", default='conf.psf', help="psf file")
parser.add_argument('-t', "--temp", default=303, type=float, help="system temperature")
parser.add_argument('-b', "--box", nargs='+', type=float, help="box dimensions in nanometer, e.g., '50 50 50' ")
parser.add_argument('-s', "--salt", default=150.0, type=float, help="salt concentration in mM")

args = parser.parse_args()
pdb_file = args.pdb
psf_file = args.psf
T = args.temp
c_ion = args.salt/1000.0                                   # concentration of ions in M
# pbc box length
if len(args.box) == 1:
    lx, ly, lz = args.box[0], args.box[0], args.box[0]
elif len(args.box) == 3:
    lx = args.box[0]*unit.nanometer
    ly = args.box[1]*unit.nanometer
    lz = args.box[2]*unit.nanometer
else:
    print("Error: You must provide either one or three values for box.")
    exit(1)
a = Vec3(lx, 0.0, 0.0)
b = Vec3(0.0, ly, 0.0)
c = Vec3(0.0, 0.0, lz)

# simulation parameters
dt = 0.0001*unit.picoseconds		                               # equilibration time step, production time step is 0.004
total_step = 250000000                                             # total step
equil_step = 10000
temperture = T*unit.kelvin                                      # temperature
log_freq = 1000                                                 # frequency of log file
dcd_freq = 5000                                                 # frequency of dcd file
pdb_freq = 12500000                                             # frequence of dpd_traj file
pressure = 1*unit.atmosphere                                    # pressure in NPT
friction = 0.1/unit.picosecond                                  # friction coefficient in Langevin

# force field parameters                
dh = 0.304/(np.sqrt(c_ion))
print('Debye-Huckel screening length: ', dh)
ffs = {
    'temp': T,                                                  # Temperature
    'dh': dh*unit.nanometer,                                    # Debye Huckel screening length
    'ke': 138.935456,                                           # Coulomb constant, ONE_4PI_EPS0
    'er': 20.0,                                         # relative dielectric constant
    'eps_hb': 2.0*unit.kilocalorie_per_mole,                    # hydrogen bond strength
    'sigma_hb': 0.29*unit.nanometer,                            # sigma of hydrogen bond
}

# 1) import coordinates and topology form charmm pdb and psf
print('\n################## load coordinates, topology and parameters ###################')
pdb = PDBFile(pdb_file)
psf = CharmmPsfFile(psf_file)
psf.setBox(lx, lx, lx)
top = psf.topology
params = CharmmParameterSet(top_inp, param_inp)
system = psf.createSystem(params, nonbondedMethod=CutoffPeriodic, constraints=HBonds)
system = HyresFF.HyresProteinSystem(psf, system, ffs)
system.setDefaultPeriodicBoxVectors(a, b, c)

with open('system.xml', 'w') as output:
    output.write(XmlSerializer.serialize(system))

print('\n# Now, the system has:')
for force in system.getForces():
    print('      ', force.getName())

# simulation
print('\n################### prepare simulation system with NVT ####################')
integrator = LangevinMiddleIntegrator(temperture, friction, dt)
plat = Platform.getPlatformByName('CUDA')
prop = {'Precision': 'mixed', 'DeviceIndex': gpu_id}
simulation = Simulation(top, system, integrator, plat, prop)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperture)
print('Langevin, CUDA, {temperature} K')

print('\n################### Minimization, Equilibriation, Production simulation ####################')
print('# minimizeEnergy:')
print('before: ', simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy(maxIterations=50000, tolerance=0.1)
print('after: ', simulation.context.getState(getEnergy=True).getPotentialEnergy())

print('\n# Equilibriation running:')
simulation.step(equil_step)
print('\nset the time step as 0.004 fs')
simulation.integrator.setStepSize(0.004*unit.picoseconds)
simulation.step(equil_step)

## save a pdb traj using large step, dcd traj using small step, and log file
simulation.reporters.append(PDBReporter('system.pdb', pdb_freq))
simulation.reporters.append(DCDReporter('system.dcd', dcd_freq))
simulation.reporters.append(StateDataReporter('system.log', log_freq, progress=True, remainingTime=True, speed=True, totalSteps=total_step, temperature=True))
#simulation.reporters.append(CheckpointReporter('system.chk', dcd_freq*10))

print('\n# NVT simulation running:')
simulation.step(total_step)

simulation.saveCheckpoint('system.chk')
print('\n# Finished!')