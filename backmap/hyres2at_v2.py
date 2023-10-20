import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.topology.guessers import guess_types
from rotamer import opt_side_chain

## this script is used to backmap HyRes model to atomistic model
## Athour: Shanlong Li

print('Hyres_rebuilder, version: 0.1.0')

inp = sys.argv[1]
out = sys.argv[2]
drt = '/home/lsl/Desktop/scripts/hyres_backmapping/map/'

hyres = mda.Universe(inp)
guessed_eles = guess_types(hyres.atoms.names)
hyres.add_TopologyAttr('elements', guessed_eles)

with mda.Writer(out, multiframe=False, reindex=False) as f:
    idx = 1
    for res in hyres.residues:
        name = res.resname
        mobile = mda.Universe(drt+name+'_ideal.pdb')
        segid = hyres.select_atoms("resid "+str(res.resid)).segids[0]
        chainID = hyres.select_atoms("resid "+str(res.resid)).chainIDs[0]
        for atom in mobile.atoms:
            atom.residue.resid = res.resid
            atom.segment.segid = segid
            atom.chainID = chainID
 
        #align.alignto(mobile.select_atoms("name N CA C O CB"), hyres.select_atoms("resid "+str(res.resid)+" and name N CA C O"), select='name N CA C O', match_atoms=False)
        align.alignto(mobile, hyres.select_atoms("resid "+str(res.resid)), select='name N CA C O', match_atoms=False)
 
        if name not in ['GLY', 'PRO', 'ALA']:
            refs = hyres.select_atoms("resid "+str(res.resid)+" and name CA CB CC CD CE CF")
            opt_side_chain(name, refs, mobile)
        
        if name != 'PRO':
            hyres_H = hyres.select_atoms("resid "+str(res.resid)+" and name H")
            for atom in hyres_H.atoms:
                atom.id = idx
                idx += 1
            f.write(hyres_H)
 
        mobile_sel = mobile.select_atoms("all")
        for atom in mobile_sel.atoms:
            atom.id = idx
            idx += 1
        #f.write(hyres_sel)
        f.write(mobile_sel)
