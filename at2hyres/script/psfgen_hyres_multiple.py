import sys
from psfgen import PsfGen


## this script is used to generate psf file for llps system
## Author: Shanlong Li @ UMass
## Usage: python *.py psf_output_filename pdb_input_filename chain_number

out = sys.argv[1]
pdb = sys.argv[2]
num = int(sys.argv[3])

gen = PsfGen()
gen.read_topology('./top_hyres_GPU.inp')

for i in range(num):
    segid = 'P'+str(i)
    gen.add_segment(segid=segid, pdbfile=pdb)

gen.write_psf(filename=out)
