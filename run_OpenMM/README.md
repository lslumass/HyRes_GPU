# This is a tutorial for running Hyres-GPU simulation using OpenMM

## 1. Installation of OpenMM
   Please follow the official instructions for installing OpenMM.
   http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm
   
## 2. Preparation
   This simulation would start from CHARMM-based input files and force-filed files, which are as follows:  
   a. coordinate of protein: PDB file  
   b. topology of protein: PSF file  
   c. HyRes-GPU top file: top_hyres_GPU.inp  
   d. HyRes-GPU parameter file: param_hyres_GPU.inp   

   HyRes-style PDB and PSF files can be created from atomistic model following the [instructions](../README.md)   
  

## 3. Running script
   run.py is the simplified script for running simulation, but [HyresBuilder] is needed.  
   **Usage**:
   ```
   python run.py -h

   usage: run.py [-h] [-c PDB] [-p PSF] [-t TEMP] [-b BOX [BOX ...]] [-s SALT]

optional arguments:
  -h, --help            show this help message and exit
  -c PDB, --pdb PDB     pdb file
  -p PSF, --psf PSF     psf file
  -t TEMP, --temp TEMP  system temperature
  -b BOX [BOX ...], --box BOX [BOX ...]
                        box dimensions in nanometer, e.g., '50 50 50'
  -s SALT, --salt SALT  salt concentration in mM
   ```
the default values are ```-c conf.pdb -p conf.psf -t 303 -s 150```   

   for details:  
   a. run.nvt.py is for the NVT simulation  
   b. run.npt.py is for the NPT simulation  
 
   **Note**: The running script is based on Python 3 and numpy is needed. Parameters should be changed inside the scripts and submit using ```python run.nvt.py conf.pdb conf.psf 0```   
   
