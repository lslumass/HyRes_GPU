import sys
from psfgen import PsfGen

## usage: python psfgen_hyres.py input_file output_file terminal_type(Neutral/charged/ACE_CT3)
## for terminal_type, input 'charged' for normal terminus, otherwise leave it blank

inp = sys.argv[1]
out = sys.argv[2]
if len(sys.argv) == 4:
    ter = sys.argv[3]
elif len(sys.argv) == 3:
    ter = 'neutral'

gen = PsfGen()
gen.read_topology('top_hyres_GPU.inp')
gen.add_segment(segid='A', pdbfile=inp)
if ter == "charged":
    res_start, res_end = gen.get_resids(segid)[0], gen.get_resids(segid)[-1]
    gen.set_charge(segid, res_start, "N", 1.00)
    gen.set_charge(segid, res_end, "O", -1.00)
gen.write_psf(filename=out)

exit()
