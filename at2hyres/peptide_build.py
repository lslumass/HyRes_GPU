"""
Simple example script demonstrating how to use the PeptideBuilder library.
Usage: python peptide_build.py seq_file
"""
import sys
from PeptideBuilder import Geometry
import PeptideBuilder

seq_file = sys.argv[1]

geo = Geometry.geometry("R")
#geo.phi = 0
#geo.psi_im1 = 0
structure = PeptideBuilder.initialize_res(geo)

with open(seq_file, 'r') as f:
    sequence = f.readlines()[0]
for seq in sequence:
    geo = Geometry.geometry(seq)
    PeptideBuilder.add_residue(structure, geo)
# add terminal oxygen (OXT)
#eptideBuilder.add_terminal_OXT(structure)

import Bio.PDB

out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save("tdp-43.pdb")

