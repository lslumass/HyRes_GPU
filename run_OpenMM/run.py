from __future__ import division, print_function
import argparse
import importlib.resources as pkg_resources
from HyresBuilder import HyresFF, utils
# OpenMM Imports
from openmm.unit import *
from openmm.app import *
from openmm import *
import numpy as np


# 0) set variables in the simulation
gpu_id = "0"
top_inp, param_inp = utils.load_ff('protein')

# input parameters
parser = argparse.ArgumentParser()
parser.add_argument('-c', "--pdb", default='conf.pdb', help="pdb file, default is conf.pdb")
parser.add_argument('-p', "--psf", default='conf.psf', help="psf file, default is conf.psf")
parser.add_argument('-t', "--temp", default=303, type=float, help="system temperature, default is 303 K")
parser.add_argument('-b', "--box", nargs='+', type=float, help="box dimensions in nanometer, e.g., '50 50 50' ")
parser.add_argument('-s', "--salt", default=0.0, type=float, help="salt concentration in mM, default is 0 mM")
parser.add_argument('-e', "--ens", default='NVT', type=str, help="simulation ensemble, NPT, NVT, or non, non is for non-periodic system")

args = parser.parse_args()
pdb_file = args.pdb
psf_file = args.psf
T = args.temp
c_ion = args.salt/1000.0                                   # concentration of ions in M
ensemble = args.ens

## set pbc and box vector
if ensemble in ['NPT', 'NVT']:
    # pbc box length
    if len(args.box) == 1:
        lx, ly, lz = args.box[0], args.box[0], args.box[0]
    elif len(args.box) == 3:
        lx = args.box[0]
        ly = args.box[1]
        lz = args.box[2]
    else:
        print("Error: You must provide either one or three values for box.")
        exit(1)
    a = Vec3(lx, 0.0, 0.0)
    b = Vec3(0.0, ly, 0.0)
    c = Vec3(0.0, 0.0, lz)
elif ensemble not in ['NPT', 'NVT', 'non']:
    print("Error: The ensemble must be NPT, NVT or non. The input value is {}.".format(ensemble))
    exit(1)

# simulation parameters
dt = 0.0001*unit.picoseconds		                               # equilibration time step, production time step is 0.004
total_step = 250000000                                             # total step
equil_step = 10000
temperture = T*unit.kelvin                                      # temperature
log_freq = 1000                                                 # frequency of log file
dcd_freq = 10000                                               # frequency of dcd file
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
top = psf.topology
params = CharmmParameterSet(top_inp, param_inp)
if ensemble == 'non':
    system = psf.createSystem(params, nonbondedMethod=CutoffNonPeriodic, constraints=HBonds)
else:
    psf.setBox(lx, ly, lz)
    top.setPeriodicBoxVectors((a, b, c))
    top.setUnitCellDimensions((lx, ly,lz))
    system = psf.createSystem(params, nonbondedMethod=CutoffPeriodic, constraints=HBonds)
    system.setDefaultPeriodicBoxVectors(a, b, c)

system = HyresFF.HyresSystem(psf, system, ffs)

with open('system.xml', 'w') as output:
    output.write(XmlSerializer.serialize(system))

print('\n# Now, the system has:')
for force in system.getForces():
    print('      ', force.getName())

# simulation
print('\n################### prepare simulation system####################')
if ensemble == 'NPT':
    print('This is a NPT system')
    system.addForce(MonteCarloBarostat(pressure, temperture, 25))
elif ensemble == 'NVT':
    print('This is a NVT system')
elif ensemble == 'non':
    print('This is a non-periodic system')
integrator = LangevinMiddleIntegrator(temperture, friction, dt)
plat = Platform.getPlatformByName('CUDA')
prop = {'Precision': 'mixed', 'DeviceIndex': gpu_id}
simulation = Simulation(top, system, integrator, plat, prop)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperture)
print(f'Langevin, CUDA, {temperture}')

print('\n################### Minimization, Equilibriation, Production simulation ####################')
print('# minimizeEnergy:')
print('before: ', simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy(maxIterations=500000, tolerance=0.01)
print('after: ', simulation.context.getState(getEnergy=True).getPotentialEnergy())

print('\n# Equilibriation running:')
simulation.step(equil_step)

## save a pdb traj using large step, dcd traj using small step, and log file
simulation.reporters.append(PDBReporter('system.pdb', pdb_freq))
simulation.reporters.append(XTCReporter('system.xtc', dcd_freq))
simulation.reporters.append(StateDataReporter('system.log', log_freq, temperature=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True))
#simulation.reporters.append(CheckpointReporter('system.chk', dcd_freq*10))

simulation.integrator.setStepSize(0.004*unit.picoseconds)
print('\n# Production simulation running:')
simulation.step(total_step)

simulation.saveCheckpoint('system.chk')
print('\n# Finished!')
