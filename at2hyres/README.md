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
   **Note:** for IDPs/IDRs, one can generate the PDB file using the PeptideBuilder package.
   TDP-43-LCD is included in the folder of examples as a simple example. Typing the following command:   
   `python peptide_build.py tdp-43-lcd.seq`   
   one can get the PDB file of TDP-43-LCD named tdp-43.pdb
4. here, one can choose any chains needed and also model missing residues if needed
5. select ACE and CT3 in the "Terminal group patching" as the two terminus
   **Note:** only ACE/CT3 are supported now in HyRes
6. Download the PDB file labelled as "CHARMM PDB"
