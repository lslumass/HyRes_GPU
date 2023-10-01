"""
Simple example script demonstrating how to use the PeptideBuilder library.
"""

from PeptideBuilder import Geometry
import PeptideBuilder

geo = Geometry.geometry("R")
#geo.phi = 0
#geo.psi_im1 = 0
structure = PeptideBuilder.initialize_res(geo)
sequence = 'NRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSWGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGWGSASNAGSGSGFNGGFGSSMDSKSSGWGM'
for seq in sequence:
    geo = Geometry.geometry(seq)
    PeptideBuilder.add_residue(structure, geo)
# add terminal oxygen (OXT)
#eptideBuilder.add_terminal_OXT(structure)

import Bio.PDB

out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save("tdp-43.pdb")

