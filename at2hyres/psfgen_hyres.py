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
gen.write_psf(filename=out)

if ter == 'charged':
    ## modify the terminus, add charge to N (+1) and O (-1)
    with open(out, 'r') as psf:
        line = psf.readlines()
        num = int(line[7].split()[0])
        line[8] = line[8][:36] + '1.000000' + line[8][44:]
        line[num+7] = line[num+7][:35] + '-1.000000' + line[num+7][44:]
        with open(out, 'w') as f:
            f.writelines(line)

exit()