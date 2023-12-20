import sys
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB

"""
Simple example script demonstrating how to use the PeptideBuilder library.
Usage: python peptide_build.py seq_file output_file_name
"""

# read each line for each protein
seq_file = sys.argv[1]
seqs = {}
with open(seq_file, 'r') as f:
    for line in f.readlines():
        name = line.split()[0]
        seq = line.split()[1].strip()
        seqs[name] = seq

# generate 
for name, seq in seqs.items():
    geo = Geometry.geometry(seq[0])
    structure = PeptideBuilder.initialize_res(geo)
    for seq in seq[1:]:
        geo = Geometry.geometry(seq)
        PeptideBuilder.add_residue(structure, geo)
    # add terminal oxygen (OXT)
    #PeptideBuilder.add_terminal_OXT(structure)

    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save(name+'.pdb')

exit()
