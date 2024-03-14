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

   HyRes-style PDB and PSF files can be created from atomistic model following the [instructions](https://github.com/wayuer19/HyRes_GPU/blob/main/at2hyres/README.md)
  

## 3. Running script
   run.py is the simplified script for running simulation  
   **Usage**: python run.py **.pdb *.psf c_ion  
   where c_ion is the salt concentration  


   for details:  
   a. run.nvt.py is for the NVT simulation  
   b. run.npt.py is for the NPT simulation  
 
   **Note**: The running script is based on Python 3 and numpy is needed.  
   
