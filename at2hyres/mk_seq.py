from modeller import *
from modeller.automodel import *

## using PDB code to obtain the sequence of protein
code = '6ssx'
e = Environ()
m = Model(e, file=code)
aln = Alignment(e)
aln.append_model(m, align_codes=code)
aln.write(file=code+'.seq')
