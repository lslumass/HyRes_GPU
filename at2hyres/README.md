# at2hyres
In this section, following the instructions, you will see how to convert an atomistic model to HyRes model.

## Dependencies
1. python 3.x (< 3.10, test with 3.8.17), conda environment recommended
2. [CHARMM-GUI](https://www.charmm-gui.org/) (registration needed)
3. psfgen-python (install through conda install -c conda-forge psfgen, python < 3.10)
4. PeptideBuilder and Biopython (Optional)

## step 1: prepare the atomistic PDB file  
HyRes was developed based CHARMM force field and software, so first a CHARMM-style PDB file is needed. Here CHARMM-GUI is a good choice.
1. open CHARMM-GUI and sign in.
2. select Input Generator --> PDB Reader&Manipulator
3. here, input the PDB id or upload the PDB file, and then NEXT>NEXT>  
   **Note:** for IDPs/IDRs, one can generate the PDB file using the PeptideBuilder package, **peptide_build.py** is a simple example.
4. 
