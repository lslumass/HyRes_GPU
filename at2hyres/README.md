# at2hyres
In this section, following the instructions, you will see how to convert an atomistic model to HyRes model.

## Dependencies
1. python 3.x (< 3.10, test with 3.8.17), conda environment recommended
2. [CHARMM-GUI](https://www.charmm-gui.org/) (registration required)
3. psfgen-python (install through "conda install -c conda-forge psfgen", python < 3.10)
4. [PeptideBuilder](https://github.com/clauswilke/PeptideBuilder) and Biopython (Optional)

# Simply starting from single-chain protein   
## step 1: prepare the atomistic PDB file  
HyRes was developed based CHARMM force field and software, so first a CHARMM-style PDB file is needed. Here CHARMM-GUI is a good choice.
1. open CHARMM-GUI and sign in.
2. select Input Generator --> PDB Reader&Manipulator
3. here, input the PDB id or upload the PDB file, and then NEXT>NEXT>  
   >[!NOTE]
   >for IDPs/IDRs, one can generate the PDB file using the PeptideBuilder package.
   TDP-43-LCD is included in the folder of examples as a simple example. Typing the following command:   
   `python peptide_build.py tdp-43-lcd.seq tdp_43.pdb`   
   one can get the PDB file of TDP-43-LCD named tdp-43.pdb.

   >[!WARNING]   
   >Sometimes after adding ACE/CT3 in the next steps, you will find the C-terminus in the pdb file has an unreasonable position of C=O, you should check this pdb file and change the last "O" to "OXT", which can be recognized by CHARMM-GUI in the next steps. But if no terminal group is patched in the next steps, nothing needs to be done.   
4. here, one can choose any chains needed and also model missing residues if needed, NEXT>
5. in this step, one can choose the type of terminal groups ("Terminal group patching").    
   >[!NOTE]   
   >For normal terminus, select NONE and NONE; for ACE/CT3 terminus, select ACE and CT3, NEXT>   
   Only normal N/C-terminus (charged or neutral) or ACE/CT3 are supported now in HyRes
6. download the PDB file labeled as "CHARMM PDB", (tdp_43_charmm.pdb in examples)

## step 2: convert atomistic pdb to HyRes pdb  
1. to make sure the PDB file are in the correct format, please fix the potential problems using **pdbfix_res.py**  
   **Usage:** `python pdbfix_res.py input_file output_file 0/1`  
   **Note:** for the last parameter, 1 is for ACE/CT3 terminus, otherwise 0. And, the chainid will be changed to "X".  
   **For example:** `python pdbfix_res.py tdp_43_charmm.pdb tdp_43_fix.pdb 1`  
3. use at2hyres_v2.py to convert the atomistic model to HyRes model  
   **Usage:** `python at2hyres_v2.py atomistic_pdb hyres_pdb`  
   **For example:** `python at2hyres_v2.py tdp_43_fix.pdb tdp_43_hyres.pdb`  


