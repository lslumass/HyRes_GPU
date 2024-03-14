from __future__ import division, print_function
from sys import argv
import numpy as np
from HyresBuilder import HyresFF
# OpenMM Imports
from openmm.unit import *
from openmm.app import *
from openmm import *

# 0) set variables in the simulation
pdb_file = argv[1]
psf_file = argv[2]
top_inp = './top_hyres_GPU.inp'
param_inp = './param_hyres_GPU.inp'

# simulation parameters
dt = 0.004*unit.picoseconds		                            # time step
total_step = 250000000                                          # total step
equil_step = 10000
log_freq = 1000                                                 # frequency of log file
dcd_freq = 5000                                                 # frequency of fine dcd file
T = 300
temperture = T*unit.kelvin                                      # temperature
c_ion = float(argv[3])/1000.0                                   # concentration of ions 
lx = 100.0*unit.nanometer                                       # pbc box length
a = Vec3(lx, 0.0, 0.0)
b = Vec3(0.0, lx, 0.0)
c = Vec3(0.0, 0.0, lx)
pressure = 1*unit.atmosphere                                    # pressure in NPT
friction = 0.1/unit.picosecond                                  # friction coefficient in Langevin
freq = 25  

# force field parameters
ffs = {
    'temp': T,
    'kf': np.sqrt(c_ion/9.480)*AngstromsPerNm,
    'er': 20.0,                                                 # effective dielectric constant
    'eps_hb': 2.0*unit.kilocalorie_per_mole,                    # hydrogen bond strength
    'sigma_hb': 0.29*unit.nanometer,                            # sigma of hydrogen bond
    'eps_base': 2.05*unit.kilocalorie_per_mole,                 # base stacking strength
    'eps_gen': 0.3*unit.kilocalorie_per_mole,                   # general pairing strength
}

# 1) import coordinates and topology form charmm pdb and psf
print('\n################## load coordinates, topology and parameters ###################')
pdb = PDBFile(pdb_file)
psf = CharmmPsfFile(psf_file)
psf.setBox(lx, lx, lx)
params = CharmmParameterSet(top_inp, param_inp)
system = psf.createSystem(params, nonbondedMethod=CutoffPeriodic, constraints=HBonds)
system = HyresFF.HyresProteinSystem(psf, system, ffs)
system.setDefaultPeriodicBoxVectors(a, b, c)

with open('system.xml', 'w') as output:
    output.write(XmlSerializer.serialize(system))

print('\n# Now, the system has:')
for force in system.getForces():
    print('      ', force.getName())

# prepare simulation
print('\n################### prepare simulation system with NVT/NPT ####################')
#system.addForce(MonteCarloBarostat(pressure, temperture, freq))     # for NPT
integrator = LangevinMiddleIntegrator(temperture, friction, dt)
plat = Platform.getPlatformByName('CUDA')
prop = {'Precision': 'mixed', 'DeviceIndex': '0'}
simulation = Simulation(psf.topology, system, integrator, plat, prop)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperture)

print('\n################### Minimization, Equilibriation, Production simulation ####################')
print('\n# minimizeEnergy:')
simulation.minimizeEnergy()
simulation.step(equil_step)

simulation.reporters.append(PDBReporter('system.pdb', total_step))
simulation.reporters.append(DCDReporter('system.dcd', dcd_freq))
simulation.reporters.append(StateDataReporter('system.log', log_freq, step=True, time=True, progress=True, totalSteps=total_step, temperature=True))
#simulation.reporters.append(CheckpointReporter('system.chk', dcd_freq*10))

print('\n# NVT simulation:')
simulation.step(total_step)

simulation.saveCheckpoint('system.chk')
print('\n# Finished!')
