# HyRes_GPU model tutorial 
HyRes model is a **Hy**brid **Res**olution coarse-grained protein model developed since 2017, where the backbone is atomistic-style description but coarse-grained side chains. For a full and complete description of the HyRes model, please move to the related papers as follows:  
1. HyRes: [Phys. Chem. Chem. Phys., 2017,19, 32421-32432](https://doi.org/10.1039/C7CP06736D)
2. HyRes-II: [J. Chem. Inf. Model. 2022, 62, 18, 4523–4536](https://doi.org/10.1021/acs.jcim.2c00974)
3. HyRes-GPU: [J. Am. Chem. Soc. 2024, 146, 1, 342–357](https://doi.org/10.1021/jacs.3c09195)  

This is a tutorial for running Hyres simulations. The tutorial contains the following:
1. IDP (intrinsically disordered protein) system: prepare PDB and PSF for a single IDP or IDPs (LLPS, liquid-liquid phase separation);
2. Folded protein system: prepare PDB and PSF for single-domain proteins and protein complexes.

Below is a short description of all folders in this tutorial:

**hyres_forcefield**: Hyres force filed file (written in the CHARMM format)

**at2hyres:** convert atomistic model to HyRes model

**run_OpenMM:** setup and script for running simulation on OpenMM   

**examples**: examples for IDP and folded protein system

**generate_psf**: generate PSF for IDP and folded proteins

**backmap:** rebuild atomistic model from HyRes results (for developers)


## Build environment  
well-tested in: OpenMM 8.0+, python 3.8/3.9
1. build a new conda environment   
conda create -n hyres python=3.8   
2. install openmm 8.0+   
reference: [OpenMM Installation](http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm)   
3. dependencies:  
a. [CHARMM-GUI](https://www.charmm-gui.org/) (registration required)   
b. [psfgen-python](https://psfgen.robinbetz.com/) (install through "conda install -c conda-forge psfgen", python < 3.10)   
c. [PeptideBuilder](https://github.com/clauswilke/PeptideBuilder) and Biopython (Optional)
d. mdtraj (>=1.9.9)

## Prepare Hyres PDB

### IDP system
For IDPs or simple peptides without folded structures, please using [HyresBuilder](https://github.com/wayuer19/HyresBuilder). The **HyresBuild** module can generate a Hyres PDB for a primary amino acid sequence of an IDP.

### Folded protein system
This is for the cases where you already have an all-atom PDB file, usually describing a folded protein or a multi-domain/chain protein. All hydrogens are not necessary but no need to remove them if your have them.

To convert atomistic model to HyRes model, please follow [README.md](at2hyres/README.md) in [at2hyres](at2hyres) for the step-by-step instructions.  


### LLPS system   
In the simulation of phase separation, there usually has several copies of protein. In this case, use [packmol](https://m3g.github.io/packmol/) or any other softwares to pack these copies into one box. And use monomer pdb file to create the psf file of the whole system ([details here](./generate_psf/README.md)). 

## Prepare Hyres PSF

We created a `psfgen_hyres.py` in the **generate_psf** folder to generate Hyres PSF for all kinds of scenarios.  
Detailed usage can be found [here](./generate_psf/README.md).   
 
## Run simulation    
Here, we prepared well-tested python scripts which incorpate the HyRes force field using OpenMM. 

1. run.nvt.py is for the NVT simulation  
2. run.npt.py is for the NPT simulation

Some parameters needed to be modified according to users' specific needs, such as box size, time step, ion concentration, and so on   

HyRes simulation can be run through:
```
python run.nvt.py pdb_file psf_file GPU_id
python run.npt.py pdb_file psf_file GPU_id
```  

