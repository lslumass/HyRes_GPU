import sys
from psfgen import PsfGen

## usage: python psfgen_hyres.py output_psf input_psf chain_number terminal_type(Neutral/charged/ACE_CT3)
## for terminal_type, input 'charged' for normal terminus, otherwise leave it blank

out = sys.argv[1]
pdb = sys.argv[2]
num = int(sys.argv[3])
ter = sys.argv[4] if len(sys.argv) > 4 else 'neutral'

gen = PsfGen()
gen.read_topology('./top_hyres_GPU.inp')

for i in range(num):
    segid = 'P'+str(i)
    gen.add_segment(segid=segid, pdbfile=pdb)
    if ter == 'charged':
        res_start, res_end = gen.get_resids(segid)[0], gen.get_resids(segid)[-1]
        gen.set_charge(segid, res_start, "N", 1.00)
        gen.set_charge(segid, res_end, "O", -1.00)

gen.write_psf(filename=out)

exit()
