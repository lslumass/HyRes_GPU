import sys
from psfgen import PsfGen

## usage: python psfgen_hyres.py input_file output_file terminal_type(Neutral/charged/ACE_CT3)
## for terminal_type, input 'charged' for normal terminus, otherwise leave it blank

pdb = sys.argv[1]
num = int(sys.argv[2])
out = sys.argv[3]
if len(sys.argv) == 5:
    ter = sys.argv[4]
elif len(sys.argv) == 4:
    ter = 'neutral'

gen = PsfGen()
gen.read_topology('./top_hyres_GPU.inp')

for i in range(num):
    segid = 'P'+str(i)
    gen.add_segment(segid=segid, pdbfile=pdb)
    if ter == 'charged':
        res_start, res_end = 1, len(gen.get_resids(segid))
        gen.set_charge(segid, res_start, "N", 1.00)
        gen.set_charge(segid, res_end, "O", -1.00)

gen.write_psf(filename=out)

exit()
