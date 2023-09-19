This is a tutorial for running Hyres-GPU simulation using OpenMM

1. Installation of OpenMM
   Please follow the official instructions for installing OpenMM.
   http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm
   
2. Preparation
   This simulation would start from CHARMM-based input files and force-filed files, which are as follows:
   a. coordinate of protein: PDB file
   b. topology of protein: PSF file
   c. HyRes-GPU top file: top-hyres-gpu.inp
   d. HyRes-GPU parameter file: param-hyres-gpu.inp

3. Running script
   a. run.nvt.py is for the NVT simulation
   b. run.npt.py is for the NPT simulation
   
   Usage: python run.*.py pdb_file psf_file GPU_id
            where GPU_id is the index of GPU device, such as 0, 1, 2
   Note: The running script is based on Python 3 and numpy is needed.
   
