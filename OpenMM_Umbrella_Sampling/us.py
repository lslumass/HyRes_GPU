import openmm as mm
import openmm.app as app
import openmm.unit as unit
from sys import stdout
import numpy as np


### system construction
T = 300*unit.kelvin
pdb = app.PDBFile('deca-ala.pdb')
forcefield = app.ForceField('amber14-all.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds, hydrogenMass=1.5*unit.amu)
integrator = mm.LangevinMiddleIntegrator(T, 1/unit.picosecond, 0.004*unit.picoseconds)
simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.reporters.append(app.DCDReporter('smd_traj.dcd', 10000))
simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, time=True, potentialEnergy=True, temperature=True, speed=True))

### equilibrate
simulation.context.setVelocitiesToTemperature(T)
simulation.step(1000)

### define the CV as the distance between the CAs of the two end residues
index1 = 8
index2 = 98
cv = mm.CustomBondForce('r')
cv.addBond(index1, index2)

### now setup SMD
# starting value
r0 = 1.3*unit.nanometers
# force constant
fc_pull = 1000.0*unit.kilojoules_per_mole/unit.nanometers**2
# pulling speed
v_pulling = 0.02*unit.nanometers/unit.picosecond # nm/ps
# simulation time step
dt = simulation.integrator.getStepSize()
# total number of steps
total_steps = 30000 # 120ps
# number of steps to run between incrementing r0 (1 makes the simulation slow)
increment_steps = 10

### define a harmonic restraint on the CV
# the location of the restrain will be moved as we run the simulation
# this is constant velocity steered MD
pullingForce = mm.CustomCVForce('0.5 * fc_pull * (cv-r0)^2')
pullingForce.addGlobalParameter('fc_pull', fc_pull)
pullingForce.addGlobalParameter('r0', r0)
pullingForce.addCollectiveVariable("cv", cv)
system.addForce(pullingForce)
simulation.context.reinitialize(preserveState=True)

### define the windows
# during the pulling loop we will save specific configurations corresponding to the windows
window_number = 24
windows = np.linspace(1.3, 3.3, window_number)
window_coords = []
window_index = 0

### SMD pulling loop
for i in range(total_steps//increment_steps):
    simulation.step(increment_steps)
    current_cv_value = pullingForce.getCollectiveVariableValues(simulation.context)
    if (i*increment_steps)%5000 == 0:
        print("r0 = ", r0, "r = ", current_cv_value)
    # increment the location of the CV based on the pulling velocity
    r0 += v_pulling * dt * increment_steps
    simulation.context.setParameter('r0',r0)
    # check if we should save this config as a window starting structure
    if (window_index < len(windows) and current_cv_value >= windows[window_index]):
        window_coords.append(simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions())
        window_index += 1

### save the window structures
for i, coords in enumerate(window_coords):
    outfile = open(f'window_{i}.pdb', 'w')
    app.PDBFile.writeFile(simulation.topology,coords, outfile)
    outfile.close()

### define running for each window
def run_window(window_index):
    print('running window', window_index)
    # load in the starting configuration for this window
    pdb = app.PDBFile(f'window_{window_index}.pdb')
    simulation.context.setPositions(pdb.positions)
    # set the fixed location of the harmonic restraint for this window
    r0 = windows[window_index]
    simulation.context.setParameter('r0', r0)
    # run short equilibration with new positions and r0
    simulation.context.setVelocitiesToTemperature(T)
    simulation.step(1000)

    # run the data collection
    # total number of steps
    total_steps = 100000 # 400 ps
    # frequency to record the current CV value
    record_steps = 1000
    # run the simulation and record the value of the CV.
    cv_values=[]
    for i in range(total_steps//record_steps):
        simulation.step(record_steps)
        # get the current value of the cv
        current_cv_value = pullingForce.getCollectiveVariableValues(simulation.context)
        cv_values.append([i, current_cv_value[0]])
    # save the CV timeseries to a file so we can postprocess
    np.savetxt(f'cv_values_window_{window_index}.txt', np.array(cv_values))
    
    print('Completed window', window_index)

### loop for running all windows
for n in range(window_number):
    run_window(n)

### generate metafile for wham analysis
metafilelines = []
for i in range(window_number):
    metafileline = f'cv_values_window_{i}.txt {windows[i]} {fc_pull}\n'
    metafilelines.append(metafileline)

with open("metafile.txt", "w") as f:
    f.writelines(metafilelines)
