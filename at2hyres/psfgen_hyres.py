import sys
from psfgen import PsfGen


inp = sys.argv[1]
out = sys.argv[2]

gen = PsfGen()
gen.read_topology('top_hyres_GPU.inp')
gen.add_segment(segid='A', pdbfile=inp)
gen.write_psf(filename=out)
