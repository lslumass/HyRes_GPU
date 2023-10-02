# at2hyres
In this section, following the instructions, you will see how to convert an atomistic model to HyRes model.

## Dependencies
1. python 3.x (< 3.10, test with 3.8.17), conda environment recommended
2. [CHARMM-GUI](https://www.charmm-gui.org/) (registration needed)
3. psfgen-python (install through conda install -c conda-forge psfgen, python < 3.10)
4. [PeptideBuilder](https://github.com/clauswilke/PeptideBuilder) and Biopython (Optional)

## step 1: prepare the atomistic PDB file  
HyRes was developed based CHARMM force field and software, so first a CHARMM-style PDB file is needed. Here CHARMM-GUI is a good choice.
1. open CHARMM-GUI and sign in.
2. select Input Generator --> PDB Reader&Manipulator
3. here, input the PDB id or upload the PDB file, and then NEXT>NEXT>  
   **Note:** for IDPs/IDRs, one can generate the PDB file using the PeptideBuilder package.
   TDP-43-LCD is included in the folder of examples as a simple example. Typing the following command:   
   `python peptide_build.py tdp-43-lcd.seq`   
   one can get the PDB file of TDP-43-LCD named tdp-43.pdb
4. here, one can choose any chains needed and also model missing residues if needed, NEXT>
5. select ACE and CT3 in the "Terminal group patching" as the two terminus, NEXT>
   **Note:** only ACE/CT3 are supported now in HyRes
6. download the PDB file labeled as "CHARMM PDB", (tdp_43_charmm.pdb in examples)

## step 2: convert atomistic pdb to HyRes pdb  
1. to make sure the PDB file are in the correct format, please fix the potential problems using **pdbfix_res.py**  
   **Usage:** `python pdbfix_res.py input_file output_file 0/1`  
   **Note:** for the last parameter, 1 is for ACE/CT3 terminus, otherwise 0.  
   **For example:** `python tdp_43_charmm.pdb tdp_43_fix.pdb 1`  
2. use at2hyres_v2.py to convert the atomistic model to HyRes model  
   **Usage:** `python at2hyres_v2.py atomistic_pdb hyres_pdb`  
   **For example:** `python at2hyres_v2.py tdp_43_fix.pdb tdp_43_hyres.pdb`  

## step 3: create psf file  
Use psfgen_hyres.py to generate psf of hyres model. Here, 'top_hyres_gpu.in' is needed.   
**Usage:** `python psfgen_hyres.py input_pbd_file outout_psf_file`  
**for example:** `python psfge_hyres.py tdp_43_hyres.pdb tdp_43_hyres.psf`


# More   
Following the above instructions, one can obtain the pdb and psf files for a single protein/peptide. But sometimes, the simulation consists of multi-copies of one protein (like LLPS simulation) or several different proteins (like protein complex or binding of folded protein and IDP). In these two cases, additional steps are required.   
## Case 1: multi-copies of one protein (for example, the LLPS simulation)   
Again, take the TDP-43-LCD as an example. Following the previous steps, we have the pdb and psf files for HyRes model of TDP-43-LCD, which are tdp_43_hyres.pdb and tdp_43_hyres.psf.   
1. To create the pdb files containing 100 copies of tdp_43_hyres.pdb, lots of tools can be used. Here, we use [packmol](https://m3g.github.io/packmol/) to pack 100 proteins together, named as 100.pdb. Both the pdb file and packmol input file (pack.inp) can be found in examples.   
