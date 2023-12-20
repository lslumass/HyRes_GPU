from __future__ import division, print_function
from sys import stdout, argv
# OpenMM Imports
from openmm.unit import *
from openmm.app import *
from openmm import *
import numpy as np

# 0) set variables in the simulation
pdb_file = argv[1]
psf_file = argv[2]
gpu_id = argv[3]
top_inp = 'top_hyres_GPU.inp'
param_inp = 'param_hyres_GPU.inp'

dt = 0.004*unit.picoseconds		                                # time step
total_step = 500000000                                             # total step
equil_step = 100000
temperture = 300*unit.kelvin                                             # temperature
log_freq = 1000                                                 # frequency of log file
dcd_freq = 50000                                                 # frequency of fine dcd file
c_ion = 0.15                                                     # concentration of ions
lx = 40.0*unit.nanometer                                          # pbc box length
a = Vec3(lx, 0.0, 0.0)
b = Vec3(0.0, lx, 0.0)
c = Vec3(0.0, 0.0, lx)

er = 20.0                                                       # effective dielectric constant
eps_hb = 2.0*unit.kilocalorie_per_mole                          # hydrogen bond strength
sigma_hb = 0.29*unit.nanometer                                  # sigma of hydrogen bond
r_cut = 1.8*unit.nanometer                                      # cutoff distance of nonbondedforce
pressure = 1*unit.atmosphere                                    # pressure in NPT
friction = 0.1/unit.picosecond                                    # friction coefficient in Langevin
freq = 25  

# 1) import coordinates and topology form charmm pdb and psf
print('load coordinates, topology and parameters')
pdb = PDBFile(pdb_file)

psf = CharmmPsfFile(psf_file)
top = psf.topology
params = CharmmParameterSet(top_inp, param_inp)
system = psf.createSystem(params, nonbondedMethod=CutoffNonPeriodic, constraints=HBonds)

# 2) constructe the force field
print('constructe the Hyres2 force field')
# get nonbonded force
for force in system.getForces():
    if force.getName() == "NonbondedForce":
        nbforce = force
print('get the NonBondedForce:', nbforce.getName())

# add custom nonbondedforce: CNBForce
formula = '(step(Ron - r) + '+ \
'(step(r - Ron) * step(Roff - r) - '+ \
'step(r - Ron) * step(Ron - r)) * '+ \
'(((Roff2 - r2)^2 * (Roff2 + 2.0 * r2 - 3.0 * Ron2)) / '+ \
'(Roff2 - Ron2)^3)) * '+ \
'(4.0 * epsilon * six * (six - 1.0) + (138.935456 / eps * charge1 * charge2) / r * exp(-kf * r));'+ \
'six = (sigma / r)^6; '+ \
'sigma = 0.5 * (sigma1 + sigma2); '+ \
'epsilon = sqrt(epsilon1 * epsilon2);'+ \
'Ron2 = Ron * Ron; Roff2 = Roff * Roff; r2 = r * r; '
CNBForce = CustomNonbondedForce(formula)
CNBForce.setNonbondedMethod(nbforce.getNonbondedMethod())
CNBForce.addGlobalParameter('eps', er)
CNBForce.addGlobalParameter('Ron', r_cut - 0.2*unit.nanometer)
CNBForce.addGlobalParameter('Roff', r_cut)
CNBForce.addGlobalParameter('kf', np.sqrt(c_ion/9.480)*AngstromsPerNm)
CNBForce.setCutoffDistance(r_cut)

# perparticle variables: sigma, epsilon, charge,
CNBForce.addPerParticleParameter('charge')
CNBForce.addPerParticleParameter('sigma')
CNBForce.addPerParticleParameter('epsilon')
for idx in range(nbforce.getNumParticles()):
    particle = nbforce.getParticleParameters(idx)
    perP = [particle[0], particle[1], particle[2]]
    CNBForce.addParticle(perP)
# add exclusion
for idx in range(nbforce.getNumExceptions()):
    ex = nbforce.getExceptionParameters(idx)
    CNBForce.addExclusion(ex[0], ex[1])
system.addForce(CNBForce)

# add nonbondedforce of 1-4 interaction through custombondforece
formula = '(step(Ron - r) + '+ \
'(step(r - Ron) * step(Roff - r) - '+ \
'step(r - Ron) * step(Ron - r)) * '+ \
'(((Roff2 - r2)^2 * (Roff2 + 2.0 * r2 - 3.0 * Ron2)) / '+ \
'(Roff2 - Ron2)^3)) * '+ \
'(4.0 * epsilon * six * (six - 1.0) + (138.935456 / eps * charge) / r * exp(-kf * r));'+ \
'six = (sigma / r)^6; '+ \
'Ron2 = Ron * Ron; Roff2 = Roff * Roff; r2 = r * r; '
Force14 = CustomBondForce(formula)
Force14.addGlobalParameter('eps', er)
Force14.addGlobalParameter('Ron', r_cut - 0.2*unit.nanometers)
Force14.addGlobalParameter('Roff', r_cut)
Force14.addGlobalParameter('kf', np.sqrt(c_ion/9.480)*10.00)

Force14.addPerBondParameter('charge')
Force14.addPerBondParameter('sigma')
Force14.addPerBondParameter('epsilon')
for idx in range(nbforce.getNumExceptions()):
    ex = nbforce.getExceptionParameters(idx)
    Force14.addBond(ex[0], ex[1], [ex[2], ex[3], ex[4]])
system.addForce(Force14)

# Add the Custom hydrogen bond force
formula  = 'epsilon*(5.0*(sigma/r)^12-6.0*(sigma/r)^10)*swrad*cosd^4*swang; '+ \
'swrad = step(rcuton-r)+step(r-rcuton)*(step(rcutoff-r)-step(rcuton-r))*'+ \
'roff2*roff2*(roff2-3.0*ron2)/roffon2^3; '+ \
'roff2 = rcutoff*rcutoff-r*r; '+ \
'ron2 = rcuton*rcuton-r*r; '+ \
'roffon2 = rcutoff*rcutoff-rcuton*rcuton; '+ \
'rcutoff = CTOFHB; rcuton = CTONHB; r = distance(a1, d2); '+ \
'swang = step(cosdcuton-cosd)+step(cosd-cosdcuton)*(step(cosdcutoff-cosd)-step(cosdcuton-cosd))*'+ \
'cosdoff2*cosdoff2*(cosdoff2-3.0*cosdon2)/cosdoffon2^3; '+ \
'cosdoff2 = cosdcutoff*cosdcutoff-cosd*cosd; '+ \
'cosdon2 = cosdcuton*cosdcuton-cosd*cosd; '+ \
'cosdoffon2 = cosdcutoff*cosdcutoff-cosdcuton*cosdcuton; '+ \
'cosdcutoff = -cos(CTOFHA); cosdcuton = -cos(CTONHA); cosd = cos(angle(a1,d1,d2));'

Hforce = CustomHbondForce(formula)
Hforce.setNonbondedMethod(nbforce.getNonbondedMethod())
Hforce.addGlobalParameter('CTOFHB', 0.5*unit.nanometers)
Hforce.addGlobalParameter('CTONHB', 0.4*unit.nanometers)
Hforce.addGlobalParameter('CTOFHA', 91*unit.degree)
Hforce.addGlobalParameter('CTONHA', 90*unit.degree)
Hforce.addGlobalParameter('sigma', sigma_hb)
Hforce.addGlobalParameter('epsilon', eps_hb)
Hforce.setCutoffDistance(0.6*unit.nanometers)

Ns, Hs, Os, Cs = [], [], [], []
for atom in psf.topology.atoms():
    if atom.name == "N" and atom.residue.name != 'PRO':
        Ns.append(int(atom.index))
    if atom.name == "H":
        Hs.append(int(atom.index))
    if atom.name == "O":
        Os.append(int(atom.index))
    if atom.name == "C":
        Cs.append(int(atom.index))

for idx in range(len(Hs)):
    Hforce.addDonor(Hs[idx], Ns[idx], -1)
    Hforce.addAcceptor(Os[idx], -1, -1)
system.addForce(Hforce)

# delete the NonbondedForce
system.removeForce(nbforce.getForceGroup())

# prepare simulation
print('prepare simulation system with NVT:')
integrator = LangevinMiddleIntegrator(temperture, friction, dt)
plat = Platform.getPlatformByName('CUDA')
prop = {'Precision': 'mixed', 'DeviceIndex': gpu_id}
simulation = Simulation(top, system, integrator, plat, prop)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperture)
print('minimizeEnergy:')
simulation.minimizeEnergy(maxIterations=50000, tolerance=0.1)

simulation.reporters.append(PDBReporter('system.pdb', total_step))
simulation.reporters.append(DCDReporter('system.dcd', dcd_freq))
#simulation.reporters.append(DCDReporter('system_fine.dcd', dcd_freq))
simulation.reporters.append(StateDataReporter('system.log', log_freq, step=True, time=True, progress=True, totalSteps=total_step, temperature=True))
#simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, progress=True, remainingTime=True, speed=True, totalSteps=total_step))

print('# NVT simulation:')
simulation.step(total_step)

simulation.saveCheckpoint('system.chk')
