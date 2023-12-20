import sys
from psfgen import PsfGen


## this script is used to combined different protein and ten generate the psf file
## the max num of protein is 10

## Author: Shanlong Li @ UMass

out = sys.argv[1]

gen = PsfGen()
gen.read_topology('top_hyres_GPU.inp')

num = len(sys.argv)
for i in range(2, num):
    if sys.argv[i] != '':
        segid = 'P'+str(i)
        gen.add_segment(segid=segid, pdbfile=sys.argv[i])

gen.write_psf(filename=out)
